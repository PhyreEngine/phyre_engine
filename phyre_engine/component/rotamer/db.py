"""
Module containing classes for building a (1997) Dunbrack-style rotamer library.
"""
import sys

import scipy
from scipy.special import i0
from scipy.stats.mstats import gmean

import numpy as np
from phyre_engine.component import Component
from phyre_engine.tools.rotamer.data.generic import NUM_CHI_ANGLES
from phyre_engine.tools.rotamer.kdtree import PeriodicKDTree
from phyre_engine.tools.rotamer.lib import BackboneDependentRotamerLibrary


class GroupRotamers(Component):
    """
    Groups residues into a dictionary indexed by amino acid type and rotamer.
    """
    # TODO: Add a code example to docstring
    REQUIRED = ["residues"]
    ADDS = ["rotamer_dict"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Collect residues into a group of rotamers."""
        residues = self.get_vals(data)
        rotamer_dict = {}
        for aa, residue_list in collect_aas(residues).items():
            rotamer_dict[aa] = collect_rotamers(residue_list)
        data["rotamer_dict"] = rotamer_dict
        return data

class RotamerScaleFactors(Component):
    """
    Calculate the pilot density at each residue for all rotamers.

    The pilot density is calculated by non-adaptive kernel density estimation
    (KDE). This class calculates the pilot density for each data point by
    KDE over all data points with the same amino acid, *not* just those data
    points adopting the same rotamer. This allows meaningful predictions even
    for sparsely-populated rotamers.
    """

    REQUIRED = ["rotamer_dict"]
    ADDS = []
    REMOVES = []

    def __init__(self, concentration, alpha = 1/2):
        self.concentration = concentration
        self.alpha = alpha

    def run(self, data, config=None, pipeline=None):
        """
        Calculate pilot densities at each data point.

        The pilot density is added to each residue in the rotamer dictionary as
        the "pilot" key. The scale factor is calculated at the same time, and
        is added to the "scale_factor" key.

        This component also calculates pilot estimates on a per-rotamer basis
        for later use when calculating rotamer χ angles. These are added to each
        residue in the "rotamer_pilot" key.
        """
        aa_dict = self.get_vals(data)
        for aa, rotamer_dict in aa_dict.items():

            # Collect data for all rotamers
            all_data_points = self.join_rotamers(rotamer_dict)
            # Get all phi, psi angles as numpy array
            all_phi_psi = to_numpy_array(
                all_data_points,
                lambda res: np.deg2rad(np.array(res["torsion"])))

            # Calculate pilot estimates for all those residues
            pilots = self.pilot_density(
                all_phi_psi,
                self.concentration[aa])
            scales = self.scale_factor(pilots)
            # Add info to the rotamer dictionary
            for residue, pilot, scale in zip(all_data_points, pilots, scales):
                residue["pilot"] = pilot
                residue["scale_factor"] = scale

            # Now calculate pilot estimates for each rotamer.
            # This data is used when calculating angles, and needs to be
            # calculated for each set of chi angles.
            #
            # That is, we need to calculate the pilot estimate for (1,), (1,1),
            # (1,2),(1,3), (1,1,1), etc.
            for chi_index in range(0, NUM_CHI_ANGLES[aa]):
                for rotamer, residue_list in collect_rotamers(
                        all_data_points, chi_index + 1).items():
                    rotamer_phi_psi = to_numpy_array(
                        residue_list,
                        lambda res: np.deg2rad(np.array(res["torsion"])))

                    rotamer_pilots = self.pilot_density(
                        rotamer_phi_psi,
                        self.concentration[aa])

                    for residue, rotamer_pilot in zip(residue_list,
                                                      rotamer_pilots):
                        if "rotamer_pilot" not in residue:
                            residue["rotamer_pilot"] = {}
                        residue["rotamer_pilot"][rotamer] = rotamer_pilot
        return data


    def pilot_density(self, data_points, concentration):
        # First, build an array of all data points. This is done because the
        # scale factors are calculated over all phi, psi pairs rather than just
        # the phi, psi pairs for each rotamer. The reasoning behind this is
        # given in Dunbrack's paper (Supplementary Information)
        pilot_estimates = np.apply_along_axis(
            prob_density_kde,
            1,
            data_points,
            data_points,
            concentration,
            None)
        """
        # For speed, split all_data_points into chunks and process in parallel
        thread_pool = ThreadPool(20)
        chunks = np.array_split(data_points, 20, 0)
        def estimate_chunk(chunk):
            return np.apply_along_axis(
                prob_density_kde,
                1,
                chunk,
                data_points,
                concentration,
                None)
        pilot_estimates = np.concatenate(
            tuple(thread_pool.map(estimate_chunk, chunks)))
        """
        return pilot_estimates

    def scale_factor(self, pilot_estimates):
        # The actual scale factors, lambda_i, are given by the geometric mean of
        # the pilot estimates divided by pilot estimate i, all raised to the
        # power of alpha.
        geom_mean = gmean(pilot_estimates)
        scale_factors = np.power(geom_mean / pilot_estimates, self.alpha)
        return scale_factors

    def join_rotamers(self, rotamer_dict):
        # Join data points of all rotamers together
        all_data_points = []
        for residue_list in rotamer_dict.values():
            all_data_points.extend(residue_list)
        return all_data_points

class CreateRotamerLibrary(Component):
    """
    Create a backbone-dependent rotamer library.

    The technique used by Dunbrack in their 2012 library is conceptually simple.
    A function is calculated giving the probability density of finding a (φ, ψ)
    angle given a known rotamer. This can be inverted using Baye's rule to give
    the desired output: the probability of finding a particular rotamer type *r*
    given a (φ, ψ) angle.

    This is done using kernel density estimation (KDE). At each data point,
    place a Gaussian function (the *kernel*): this gives an expression for the
    probability density of (φ, ψ) as a sum over all data points of that
    particular rotamer type. In reality, we use the von Mises kernel rather than
    the more common Gaussian kernel; this accounts for the periodic nature of
    φ and ψ.

    If the width (scale factor) of each kernel is constant, then sparse areas
    of (φ, ψ) space will contain strange peaks. We use adaptive KDE (AKDE) to
    compensate for this by increasing the width of the kernels placed upon
    data points in sparse regions. The width of each kernel is determined by the
    scale factor.

    This component requires the key ``rotamer_dict`` to be present in the
    pipeline state. The dictionary ``rotamer_dict`` should be indexed by amino
    acid type. Each value should be another dictionary, indexed by rotamer name
    and pointing to a list of dictionaries containing information about
    residues. Each residue should have a scale factor assigned in order to
    determine the width of the kernel at each data point. Scale factors can be
    assigned using the
    :py:class:`~phyre_engine.component.rotamer.db.RotamerScaleFactors`
    component.
    """

    REQUIRED = ["rotamer_dict"]
    ADDS = ["rotamer_library"]
    REMOVES = []

    def __init__(self, concentrations):
        """
        :param concentrations: Dictionary mapping amino acid names to
            concentration factors. Concentration factors alter the width of each
            kernel.
        """
        self.concentrations = concentrations

    def adaptive_prob_density(self, phi, psi, residue_list, concentration):
        r"""
        Calculate the adaptive probability density at a given φ, ψ point.

        Each residue in ``residue_list`` should adopt the same rotamer. This method
        returns the probability density of seeing a residue in that rotamer with
        the given :math:`\phi` and :math:`\psi` angles.

        :param float phi: Backbone :math:`\phi` angle in radians.
        :param float psi: Backbone :math:`\psi` angle in radians.
        :param list residue_list: List of residues all adopting the same
            rotamer.
        :param float concentration: Concentration used to determine the
            bandwidth of the kernels used when estimating kernel density.

        :return: Probability density.
        :rtype: float
        """
        phi_psi_array = to_numpy_array(
            residue_list,
            lambda residue: np.deg2rad(np.array(residue["torsion"])))

        scale_factor_array = to_numpy_array(
            residue_list,
            lambda residue: residue["scale_factor"])

        return prob_density_kde(
                np.array([phi, psi]),
                phi_psi_array,
                concentration,
                scale_factor_array)

    def rotamer_prob_densities(self, phi, psi, rotamer_dict, concentration):
        r"""
        Calculate probability densities for all rotamers.

        Calculates the probability density of seeing residues with backbone
        angles given by ``phi`` and ``psi``. The probability density is
        calculated for each rotamer in ``rotamer_dict``.

        :param float phi: Backbone :math:`\phi` angle in radians.
        :param float psi: Backbone :math:`\psi` angle in radians.
        :param rotamer_dict: Dictionary containing lists of residues indexed by
            rotamer name.
        :param float concentration: Base concentration parameter describing the
            bandwidth of each kernel.

        :return: Dictionary of probability densities indexed by rotamer name.
        """
        prob_densities = {}
        for rotamer, residue_list in rotamer_dict.items():
            prob_densities[rotamer] = self.adaptive_prob_density(
                phi, psi, residue_list, concentration)
        return prob_densities

    def prob_of_rotamer(self, rotamer, rotamer_probs, prob_dens_dict):
        r"""
        Calculate probability of rotamers at a particular (:math:`\phi, \psi`)
        angle.

        The probability :math:`P(r|\phi, \psi)` of each rotamer :math:`r` in the
        dictionary ``rotamer_dict`` for a residue with backbone angles
        :math:`(\phi, \psi)` is given by the expression

        .. math ::

            P(r|\phi, \psi) = \frac{%
                P(r) \rho(\phi, \psi | r)
            }{%
                \sum_{r'} P(r') \rho(\phi, \psi | r')
            },

        where :math:`P(r)` gives the total probability of rotamer :math:`r`  (as
        calculated by :py:meth:`.rotamer_probs` and passed in the dictionary
        ``rotamer_probs``) and :math:`\rho(\phi, \psi | r)` gives the
        probability density of seeing a residue with rotamer :math:`r` with the
        backbone angles :math:`(\phi, \psi)` (calculated by
        :py:meth:`.adaptive_prob_density` and passed in the parameter
        ``prob_dens_dict``).

        :param rotamer: Name of the rotamer for which the probability is
            calculated.
        :param rotamer_probs: Dictionary indexed by rotamer name giving the
            total probability of seeing each rotamer. Calculated by
            :py:meth:`.rotamer_probs`.
        :param prob_dens_dict: Dictionary indexed by rotamer name giving the
            probability density of seeing a residue with particular backbone
            angles. Calculated by :py:meth:`.adaptive_prob_density`.

        :return: Probability of seeing ``rotamer``.
        """

        # Calculate the probability of "rotamer" at backbone angle (phi, psi)
        # using Baye's rule.
        rotamer_prob_density = prob_dens_dict[rotamer]
        numerator = rotamer_probs[rotamer] * rotamer_prob_density

        denominator = 0
        for r in rotamer_probs.keys():
            prob_dens = prob_dens_dict[r]
            denominator += prob_dens * rotamer_probs[r]
        return numerator / denominator

    def rotamer_probs(self, rotamer_dict, num_data_points):
        r"""
        Calculate probability of this rotamer occurring.

        The probability :math:`P(r)` of rotamer `r` (for a particular amino) is
        simply given by

        .. math::

            P(r) = \frac{N_r}{\sum_{r'} N_{r'}}.
        """
        # Probability of each rotamer
        return {r: len(s) / num_data_points for r, s in rotamer_dict.items()}

    def run(self, data, config=None, pipeline=None):
        r"""
        Calculate probability of a residue adopting a particular rotamer at a
        range of :math:`\phi` and :math:`\psi` angles.
        """
        residue_dict = self.get_vals(data)
        rotamer_library = BackboneDependentRotamerLibrary(residue_dict)

        for aa, rotamer_dict in residue_dict.items():
            # Calculate total number of data points for calculating rotamer
            # probability.
            data_points_per_aa = np.sum([len(r) for r in rotamer_dict.values()])

            # Get the probability of each rotamer, indexed by rotamer
            rotamer_probs = self.rotamer_probs(rotamer_dict, data_points_per_aa)

            for phi  in range(-180, 190, 10):
                for  psi in range(-180, 190, 10):
                    print("Calculating rotamer probabilities for {aa} at ({phi}, {psi})".format(
                        aa=aa,
                        phi=phi,
                        psi=psi
                    ), file=sys.stderr)

                    probability_densities = self.rotamer_prob_densities(
                        np.deg2rad(phi), np.deg2rad(psi),
                        rotamer_dict, self.concentrations[aa])

                    for rotamer in rotamer_dict.keys():
                        rotamer_prob = self.prob_of_rotamer(
                            rotamer, rotamer_probs, probability_densities)
                        rotamer_bin = rotamer_library[aa][(phi, psi)][rotamer]
                        rotamer_bin.probability = rotamer_prob
        data["rotamer_library"] = rotamer_library
        return data

class IndependentRotamerAngles(Component):
    REQUIRED = ["residues", "rotamer_library"]
    ADDS = []
    REMOVES = []

    def collect_aas(self, residues):
        """Create dictionary of residues, indexed by residue type."""
        aa_dict = {}
        for residue in residues:
            res_name = residue["residue"].get_resname()
            if res_name not in aa_dict:
                aa_dict[res_name] = []
            aa_dict[res_name].append(residue)
        return aa_dict

    def collect_rotamers(self, residues, chi_index):
        """Create dict of side chains indexed by rotamer"""
        rotamer_dict = {}
        for residue in residues:
            if residue["rotamer"] is None:
                continue

            rot = residue["rotamer"].name[0:chi_index + 1]
            if rot not in rotamer_dict:
                rotamer_dict[rot] = []
            rotamer_dict[rot].append(residue)
        return rotamer_dict

    def run(self, data, config=None, pipeline=None):
        residues = data["residues"]
        rotamer_library = data["rotamer_library"]

        # Collect bins by amino acid and rotamer
        bin_dict = {}
        for rot_bin in rotamer_library:
            if rot_bin.aa not in bin_dict:
                bin_dict[rot_bin.aa] = {}
            if rot_bin.rotamer not in bin_dict[rot_bin.aa]:
                bin_dict[rot_bin.aa][rot_bin.rotamer] = []
            bin_dict[rot_bin.aa][rot_bin.rotamer].append(rot_bin)

        for aa, res_list in self.collect_aas(residues).items():
            rot_angles = {name: [None] * NUM_CHI_ANGLES[aa]
                          for name in bin_dict[aa]}

            for chi_index in range(0, NUM_CHI_ANGLES[aa]):
                rotamer_dict = self.collect_rotamers(res_list, chi_index)
                for rotamer, rotamer_res_list in rotamer_dict.items():
                    angles = np.array([r["sidechain"].angles[chi_index]
                                       for r in rotamer_res_list])
                    mean = scipy.stats.circmean(angles, high=180, low=-180)

                    for angle_key in rot_angles.keys():
                        if angle_key[0:chi_index + 1] == rotamer:
                            rot_angles[angle_key][chi_index] = mean

            for rotamer, mean in rot_angles.items():
                for rot_bin in bin_dict[aa][rotamer]:
                    rot_bin.mean = tuple(mean)
        return data

class CalculateRotamerAngles(Component):
    REQUIRED = ["residues", "rotamer_dict", "rotamer_library"]
    ADDS = []
    REMOVES = []

    def __init__(self, concentrations):
        self.concentrations = concentrations

    def scale_parameter(
            self,
            phi, psi,
            rotamer, residues,
            concentration, alpha,
            quadtree):
        r"""
        Calculates query-dependent scale parameter at (``phi``, ``psi``).

        When calculating the mean χ angles of a rotamer at a particular (φ, ψ)
        point, each χ angle is calculated as a weighted mean:

        .. math::

            \mu(\chi|\phi, \psi, r) = \frac{%
                \sum_i^{N_r} K_m(\phi - \phi_i) K_m(\psi - \psi_i) \chi_i
            }{%
                \sum_i^{N_r} K_m(\phi - \phi_i) K_m(\psi - \psi_i)
            }.

        The scale factor for the kernels :py:math:`K_m` is evaluated once at
        each query (φ, ψ) point. We follow the method of Shapovalov and Dunbrack
        (2010) and ensure that the bandwidth of the kernel encompasses at least
        25 points.
        """

        pilot_estimates = to_numpy_array(
            residues,
            lambda res: res["rotamer_pilot"][rotamer])
        phi_psi_list = to_numpy_array(
            residues,
            lambda res: np.deg2rad(np.array(res["torsion"])))

        # Get the non-adaptive kde at the query point
        query_estimator = prob_density_kde(
            np.array([phi, psi]),
            phi_psi_list,
            concentration)
        scale_param = np.power(gmean(pilot_estimates) / query_estimator, alpha)

        # We want to expand the scale parameter in sparse regions.
        # We do this by following Shapovalov and Dunbrack.
        # First, we find the nearest 25 points. If the distance of the farthest
        # point is less than our required distance, we accept the current scale
        # parameter. Otherwise, we take the distance of the 25th point and
        # convert that into our scale factor.

        # The conversion between distances r and scale factors λ is given by
        # r = 1/√(k / λ). So λ = kr^2
        cutoff_distance = np.sqrt(scale_param / concentration)
        distances, _ = quadtree.query(
            np.array([phi, psi]),
            25)

        if distances[-1] >= cutoff_distance:
            scale_param = concentration * distances[-1] ** 2

        return scale_param

    def rotamer_mean_chi(
            self,
            phi, psi,
            chi_index,
            scale_parameter,
            residues, concentration,
            quadtree):
        """
        Calculates the mean χ angle for a given rotamer at a particular point
        in (φ, ψ) space.

        In particularly dense regions calculating the kernel width may cause an
        overflow or underflow. In these regions, we look up all points within
        10 degrees of the query point and simply take the mean χ angles of those
        points. If this occurs, a warning will be emitted by numpy.
        """
        phi_psi_list = to_numpy_array(
            residues,
            lambda res: np.deg2rad(np.array(res["torsion"])))

        chi_list = to_numpy_array(
            residues,
            lambda res: np.deg2rad(np.array(res["sidechain"].angles[chi_index])))

        # Expand scale so it's the same length as the number of residues
        scale_list = np.repeat(scale_parameter, len(residues))

        kde_list = kernel(
            np.array([phi, psi]),
            phi_psi_list,
            concentration,
            scale_list)
        mean = np.sum(kde_list * chi_list) / np.sum(kde_list)

        if not np.isfinite(mean):
            # In dense regions, there is a good chance to overflow our floats.
            # Because these regions are densely populated, we can get a
            # reasonable estimate of the mean chi angle by just taking the mean
            # chi angle of all data points around our query.
            _, close_indices = quadtree.query(
                np.array([phi, psi]),
                len(chi_list), distance_upper_bound=np.deg2rad(10))
            # The final element of "_" will be inf because the underlying kdtree
            # implementation represents all points outside of the allowed upper
            # bounds with inf. The corresponding index is len(_), which does not
            # exist, so we will only use up to the last index
            mean = np.mean(chi_list[close_indices[0:-1]])

        return mean

    def run(self, data, config=None, pipeline=None):
        all_residues, rotamer_dict, rotamer_library = self.get_vals(data)

        for aa, aa_residues in collect_aas(all_residues).items():
            # Indexed by ((phi, psi), (chi1, chi2...)).
            # Called "partial_means" because the value is a single chi value
            # for the final rotating unit of the (sub) rotamer.
            partial_means = {}

            for chi_index in range(0, NUM_CHI_ANGLES[aa]):
                for rotamer, residues in collect_rotamers(aa_residues, chi_index + 1).items():
                    # Build a kd-tree for each rotamer
                    rot_bb = to_numpy_array(
                        residues,
                        lambda r: np.deg2rad(np.array(r["torsion"])))
                    quadtree = PeriodicKDTree(rot_bb, 2 * np.pi, 2 * np.pi)

                    for phi in range(-180, 190, 10):
                        for  psi in range(-180, 190, 10):
                            scale_param = self.scale_parameter(
                                np.deg2rad(phi), np.deg2rad(psi),
                                rotamer, residues,
                                self.concentrations[aa], 0.5,
                                quadtree)
                            mean = self.rotamer_mean_chi(
                                np.deg2rad(phi), np.deg2rad(psi),
                                chi_index, scale_param,
                                residues, self.concentrations[aa],
                                quadtree)
                            # Assign mean chi angle to all matching chis
                            key = ((phi, psi), rotamer)
                            partial_means[key] = np.rad2deg(mean)

            # Go from "partial" means (i.e. means calculated on sub-rotamers) to
            # full sets of means (r1, r2, r3, r4 => chi1, chi2, chi3, chi4)
            for phi in range(-180, 190, 10):
                for psi in range(-180, 190, 10):
                    for rotamer in rotamer_dict[aa].keys():
                        self._set_bin_mean(
                            rotamer_library, aa, partial_means,
                            rotamer, phi, psi)
        return data

    def _set_bin_mean(
            self, rotamer_library, aa, partial_means,
            rotamer, phi, psi):
        """Set the mean of rotamer_library[aa][(phi, psi)][rotamer]."""

        chi_values = [None] * NUM_CHI_ANGLES[aa]
        for chi_index in range(0, NUM_CHI_ANGLES[aa]):
            sub_rotamer = rotamer[0:chi_index + 1]
            chi_values[chi_index] = partial_means[(phi, psi), sub_rotamer]
            rotamer_bin = rotamer_library[aa][phi, psi][rotamer]
            rotamer_bin.mean_chi = chi_values


class WriteRotamerLibrary(Component):
    """
    Write a rotamer library in the format used by Shapovalov and Dunbrack.

    The library consists of the following columns, separated by at least one
    space: T, Phi, Psi, Count, r1, r2, r3, r4, Probabil, chi1Val, chi2Val,
    chi3Val, chi4Val, chi1Sig, chi2Sig, chi3Sig and chi4Sig. These are defined
    as:

    T:
        Three-letter amino acid name.

    Phi, Psi:
        Backbone dihedral angles in degrees.

    Count:
        Number of residues of this amino acid type with backbone angles in the
        current bin.

    r1--r4:
        Numbers indicating rotamer type. Each number gives the rotamer about
        the corresponding χ angle.

    Probabil:
        The probability of residues falling in this (ϕ, ψ) bin adopting this
        rotamer.

    chi1Val--chi4Val:
        Mean χ angles of residues adopting this rotamer in this bin.

    chi1Sig--chi4Sig:
        Standard deviation of the χ angles of residues adopting this rotamer in
        this bin.

    Rotamers Within each (ϕ, ψ) bin will be sorted by probability.

    .. warning::
        We don't currently calculate the count or standard deviation fields, so
        these will be written with a dummy value of zero.
    """
    REQUIRED = ["rotamer_library"]
    ADDS = []
    REMOVES = []
    # ARG  -180 -180    10     1  2  2  1  0.249730    62.5   176.9   176.6    85.7       6.9    11.1    10.5     9.9

    OUTPUT_LINE = (
        "{AA:3s} {phi:4d} {psi:4d} {count:5d}    "
        "{rotamer[0]:2d} {rotamer[1]:2d} {rotamer[2]:2d} {rotamer[3]:3d} "
        "{probability: 8.6f} "
        "{mean[0]: 6.1f} {mean[1]: 6.1f} {mean[2]: 6.1f} {mean[3]: 6.1f} "
        "{sd[0]: 6.1f} {sd[1]: 6.1f} {sd[2]: 6.1f} {sd[3]: 6.1f} ")

    def __init__(self, output_file):
        """
        :param str output_file: File to write.
        """
        self.output_file = output_file

    def run(self, data, config=None, pipeline=None):
        rotamer_library = self.get_vals(data)
        with open(self.output_file, "w") as out_fh:
            for aa, aa_dict in rotamer_library.items():
                # Sort phi, psi by phi and then psi
                sorted_backbone = sorted(
                    aa_dict.items(),
                    key=lambda kv: kv[0])
                for (phi, psi), rotamer_dict in sorted_backbone:
                    # Sort rotamer bins by probability
                    sorted_bins = sorted(
                        rotamer_dict.items(),
                        key=lambda kv: kv[1].probability,
                        reverse=True)
                    for rotamer, rot_bin in sorted_bins:
                        # Extend mean and rotamer name to 4 digits
                        rotamer = rotamer + (0,) * (4 - len(rotamer))
                        mean = (rot_bin.mean_chi
                            + [0] * (4 - len(rot_bin.mean_chi)))

                        line = self.OUTPUT_LINE.format(
                            AA=aa, phi=phi, psi=psi, count=0, rotamer=rotamer,
                            probability=rot_bin.probability,
                            mean=mean, sd=[0] * 4)
                        print(line, file=out_fh)
        return data



def kernel(phi_psi, phi_psi_data, concentration, scale=None):
    # If scale is supplied, each element is multiplied with each element before
    # the sum.
    if scale is None:
        scale = np.ones(len(phi_psi_data))

    # Calculate cos(phi - phi_i) + cos(psi - psi_i) for each element
    arg = np.cos(phi_psi_data - phi_psi)
    # Do the residue-wise sum: cos(...) + cos(...)
    arg = np.sum(arg, 1)
    # Multiply by concentration constant (kappa)
    arg = arg * concentration / scale
    # Calculate exponential for each element
    arg = np.exp(arg)
    # Multiply by constant (1/i0(kappa))^2
    arg = arg / np.square(i0(concentration / scale))
    # Multiply by constant (1/(4 * pi^2))
    arg = arg / (4 * np.square(np.pi))
    return arg

def prob_density_kde(phi_psi, phi_psi_data, concentration, scale=None):
    arg = kernel(phi_psi, phi_psi_data, concentration, scale)
    arg = np.sum(arg, 0)
    # Multiply by constant (1/N_r))
    arg = arg / len(phi_psi_data)
    return arg

def collect_aas(residues):
    """
    Create dictionary of residues, indexed by residue type.

    Given a list of dictionaries from the ``residues`` key of the pipeline data,
    build a dictionary indexed by residue type. Each element of the dictionary
    will be a list of residues.

    :param list residues: List of residues to collate.
    """
    aa_dict = {}
    for residue in residues:
        res_name = residue["residue"].get_resname()
        if res_name not in aa_dict:
            aa_dict[res_name] = []
        aa_dict[res_name].append(residue)
    return aa_dict

def collect_rotamers(residues, max_chi=None):
    """
    Create dictionary of residues indexed by rotamer.

    Given a list of dictionaries from the ``residues`` key of the pipeline data,
    build a dictionary indexed by the rotamer name that contains lists of all
    residues adopting that rotamer.

    If the parameter ``max_chi`` is provided, only ``max_chi`` chi angles will
    be considered. That is, if ``max_chi`` is ``1``, residues beloinging to
    rotamers ``(1, 1)`` and ``(1, 2)`` will be grouped into a rotamer named
    ``(1,)``. If you are using a rotamer scheme that doesn't support indexing of
    the rotamer names, this method will fail.

    .. warning::

        Note that if your rotamers are named with strings and each chi angle is
        not represented by a single character, this will produce undesirable
        results. For example, if ``max_chi=2``, the rotamers ``m-25`` and
        ``m-120`` will be grouped together as ``m-``.

    :raises TypeError: If rotamer names are not indexable.
    :param list residues: List of residues to collate.
    :param int max_chi: Maximum number of angles to consider.
    :return: Dictionary of rotamer names mapped to lists of residues.
    """
    rotamer_dict = {}
    for residue in residues:
        if residue["rotamer"] is None:
            continue

        rot = residue["rotamer"].name
        if max_chi is not None:
            rot = rot[0:max_chi]
        if rot not in rotamer_dict:
            rotamer_dict[rot] = []
        rotamer_dict[rot].append(residue)
    return rotamer_dict

def to_numpy_array(array, selector):
    """
    Convert a list of arbitrary values to a numpy array.

    Elements of ``array``  will be mapped by calling the function ``selector``
    with the element as an argument. The selector should return a scalar or
    numpy array.

    :param list array: List of elements to map.
    :param function selector:  Function to map elements to scalar or numpy
        arrays.
    :return: Numpy array.
    """
    return np.array([selector(e) for e in array])
