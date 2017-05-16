"""
Module containing classes for building a (1997) Dunbrack-style rotamer library.
"""
import math
import sys
import pickle
from pathlib import Path
import numpy as np
from scipy.special import i0
from scipy.stats.mstats import gmean
from multiprocessing.pool import ThreadPool
from phyre_engine.component import Component

class Bin:
    def __init__(self, aa, dihedrals, rotamer, prob):
        self.aa = aa
        self.dihedrals = dihedrals
        self.rotamer = rotamer
        self.prob = prob

    def __repr__(self):
        return "<Bin({aa}, {dihedrals}, {rotamer}, {prob})>".format(
            aa=self.aa,
            dihedrals=self.dihedrals,
            rotamer=self.rotamer,
            prob=self.prob)

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

    This component calculates scale factors, probability densities and then
    uses Baye's rule to calculate the actual rotamer probabilities.
    """

    REQUIRED = ["residues"]
    ADDS = []
    REMOVES = []

    def __init__(self, concentrations):
        """
        :param concentrations: Dictionary mapping amino acid names to
            concentration factors. Concentration factors alter the width of each
            kernel.
        """
        self.concentrations = concentrations

    def collect_aas(self, residues):
        """Create dictionary of residues, indexed by residue type."""
        aa_dict = {}
        for residue in residues:
            res_name = residue["residue"].get_resname()
            if res_name not in aa_dict:
                aa_dict[res_name] = []
            aa_dict[res_name].append(residue)
        return aa_dict

    def collect_rotamers(self, residues):
        """Create dict of side chains indexed by rotamer"""
        rotamer_dict = {}
        for residue in residues:
            if residue["rotamer"] is None:
                continue

            rot = residue["rotamer"].name
            if rot not in rotamer_dict:
                rotamer_dict[rot] = []
            rotamer_dict[rot].append(residue)
        return rotamer_dict

    def rotamer_probs(self, rotamer_dict, num_data_points):
        """Calculate probability of this rotamer occuring."""
        # Probability of each rotamer
        return {r: len(s) / num_data_points for r, s in rotamer_dict.items()}

    def to_numpy_array(self, rotamer_dict):
        """Build a dictionary of rotamers containing (φ, ψ) as numpy arrays."""
        result_dict = {}
        for rotamer, res_list in rotamer_dict.items():
            # The numpy array only needs to contain phi and psi
            backbone = np.array([res["torsion"] for res in res_list])
            result_dict[rotamer] = np.deg2rad(backbone)
        return result_dict

    def adaptive_prob_density(self,
            phi, psi,
            rotamer, rotamer_dict,
            scale_factors, concentration):
        """Calculate the adaptive probability density at a given φ, ψ point."""

        return prob_density_kde(
                np.array([phi, psi]),
                rotamer_dict[rotamer],
                concentration,
                scale_factors[rotamer])

    def calculate_scale_parameters(self, rotamer_dict, concentration, alpha):
        """
        Calculate each scale parameter λ_i for every phi_psi data point.

        The pilot estimate is given by prob_density_kde with a scale of 1.
        """

        # First, build an array of all data points. This is done because the
        # scale factors are calculated over all phi, psi pairs rather than just
        # the phi, psi pairs for each rotamer. The reasoning behind this is
        # given in Dunbrack's paper (Supplementary Information)
        all_data_points = np.concatenate(
            tuple(v for v in rotamer_dict.values()),
            0)

        # For speed, split all_data_points into chunks and process in parallel
        thread_pool = ThreadPool(20)
        chunks = np.array_split(all_data_points, 20, 0)
        def estimate_chunk(chunk):
            return np.apply_along_axis(
                prob_density_kde,
                1,
                chunk,
                all_data_points,
                concentration,
                None)

        pilot_estimates = np.concatenate(
            tuple(thread_pool.map(estimate_chunk, chunks)))

        # The actual scale factors, lambda_i, are given by the geometric mean of the
        # pilot estimates divided by pilot estimate i, all raised to the power of
        # alpha.
        geom_mean = gmean(pilot_estimates)
        scale_factors = np.power(geom_mean / pilot_estimates, alpha)

        # Finally, split the list back into rotamers.
        i = 0
        pilot_estimate_dict = {}
        for rotamer, values in rotamer_dict.items():
            rotamer_scale_factors = scale_factors[i:i + len(values)]
            i += len(values)
            pilot_estimate_dict[rotamer] = rotamer_scale_factors
        return pilot_estimate_dict

    def prob_densities(
            self, rotamer_dict, phi, psi, scale_factors, concentration):
        # Calculate the probability density at phi, psi for all rotamer types
        # Do this in one step because all rotamers are included in the
        # denonominator when using Baye's rule to invert probabilities later,
        # so these should be cached.
        # Returns a dictionary mapping rotamers to probability density.
        prob_dens_dict = {}
        for r in rotamer_dict.keys():
            prob_dens = self.adaptive_prob_density(
                    phi, psi,
                    r, rotamer_dict,
                    scale_factors, concentration)
            prob_dens_dict[r] = prob_dens
        return prob_dens_dict

    def prob_of_rotamer(self,
            rotamer,
            rotamer_dict, rotamer_probs,
            prob_dens_dict,
            scale_factors, concentration):

        # Calculate the probability of "rotamer" at backbone angle (phi, psi)
        # using Baye's rule.
        rotamer_prob_density = prob_dens_dict[rotamer]
        numerator = rotamer_probs[rotamer] * rotamer_prob_density

        denominator = 0
        for r in rotamer_dict.keys():
            prob_dens = prob_dens_dict[r]
            denominator += prob_dens * rotamer_probs[r]
        return numerator / denominator

    def run(self, data):
        all_residues = self.get_vals(data)

        # Probabilities of each rotmer in each phi, psi bin
        bin_probs = []
        for aa, residues in self.collect_aas(all_residues).items():
            print("Calculating for {}".format(aa), file=sys.stderr)

            # Split into rotamers
            rotamer_dict = self.collect_rotamers(residues)
            # Get the probability of each rotamer, indexed by rotamer
            rotamer_probs = self.rotamer_probs(rotamer_dict, len(residues))
            # Convert rotamer_dict to numpy arrays of phi, psi angles
            rotamer_dict = self.to_numpy_array(rotamer_dict)
            # Calculate scale factors for each data point. This is expensive.
            scale_factors = self.calculate_scale_parameters(
                rotamer_dict,
                self.concentrations[aa],
                0.5)


            # For each rotamer, calculate prob density
            for rotamer  in rotamer_dict.keys():
                for phi  in range(-180, 190, 10):
                    for  psi in range(-180, 190, 10):
                        prob_densities = self.prob_densities(
                            rotamer_dict,
                            np.deg2rad(phi), np.deg2rad(psi),
                            scale_factors, self.concentrations[aa])

                        prob = self.prob_of_rotamer(
                                rotamer,
                                rotamer_dict, rotamer_probs,
                                prob_densities,
                                scale_factors,
                                self.concentrations[aa])
                        bin_probs.append(Bin(aa, (phi, psi), rotamer, prob))
        data["rotamer_library"] = bin_probs
        return data


def prob_density_kde(phi_psi, phi_psi_data, concentration, scale=None):
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
    # Sum arguments
    arg = np.sum(arg, 0)
    # Multiply by constant (1/(4 * pi * N_r))
    arg = arg / (4 * np.square(np.pi) * len(phi_psi_data))
    return arg
