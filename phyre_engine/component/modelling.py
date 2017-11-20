"""Module providing utilities for homology modelling."""

import json
import Bio.SeqUtils
import Bio.PDB.Structure
import Bio.PDB.Model
import Bio.PDB.Chain
import Bio.PDB.PDBIO
import Bio.PDB.Residue
from phyre_engine.component import Component
import phyre_engine.tools.pdb as pdb
from phyre_engine.tools.template import Template
import numpy as np
from phyre_engine.tools.external import ExternalTool
import tempfile
import re
import subprocess
import shutil
from pathlib import Path

BACKBONE_ATOMS = {"N", "C", "CA", "O"} # Set of backbone atoms

class HomologyModeller(Component):
    """
    Generate a model for each element of the ``templates`` field in the pipeline
    state.

    This modeller requires pre-generated PDB files for each chain in the
    ``templates`` list. Each PDB template *must* contain the ``REMARK`` fields
    written by :py:class:`phyre_engine.component.db.db.ChainPDBBuilder` so that
    the sequence alignment can be correctly mapped onto the PDB structure. Each
    template must contain an ``alignment`` key containing the mapping between
    the query sequence and the template sequence.

    This component will add the ``model`` key to each template in the pipeline
    state. This will contain a file path pointing to the generated model.

    :param str chain_dir: Top-level directory containing the PDB chains.
    :param str model_name: Python template string, formatted with all the
        elements of each template.

    .. seealso::

        :py:class:`phyre_engine.component.db.db.ChainPDBBuilder`
            For generating PDB files of the correct format for use as templates.
    """

    REQUIRED = ["templates", "sequence"]
    ADDS = ["model"]
    REMOVES = []

    def __init__(self, chain_dir, model_name="{rank:02d}-{PDB}_{chain}.pdb"):
        self.chain_dir = chain_dir
        self.model_name = model_name

    def run(self, data, config=None, pipeline=None):
        """Build a model."""
        templates, query_seq = self.get_vals(data)
        pdb_io = Bio.PDB.PDBIO()

        for template in templates:
            alignment = template["alignment"]
            pdb_id = template["PDB"]
            chain = template["chain"]


            model_file = self.model_name.format(**template)

            if not Path(model_file).exists():
                self.logger.debug("Creating model file %s", model_file)

                template_file = pdb.pdb_path(
                    pdb_id, ".pdb", chain, self.chain_dir)
                db_template = Template.load(template_file)

                model_name = "model from {}_{}".format(pdb_id, chain)
                model_structure = Bio.PDB.Structure.Structure(model_name)
                model_model = Bio.PDB.Model.Model(1)
                model_chain = Bio.PDB.Chain.Chain("A")
                model_model.add(model_chain)
                model_structure.add(model_model)

                for residue_pair in alignment:
                    # Residue indices
                    i, j = residue_pair[0:2]

                    # Residue ID from canonical sequence map
                    j_id = db_template.canonical_indices[j - 1]

                    # Template residue
                    template_res = db_template.chain[j_id]

                    query_res_type = Bio.SeqUtils.seq3(query_seq[i - 1]).upper()
                    query_res = Bio.PDB.Residue.Residue(
                        (" ", i, " "),
                        query_res_type, " ")
                    for atom in template_res:
                        if atom.get_name() in BACKBONE_ATOMS:
                            query_res.add(atom.copy())
                    model_chain.add(query_res)

                pdb_io.set_structure(model_structure)
                pdb_io.save(model_file)
            template["model"] = model_file
        return data

class SoedingSelect(Component):
    """
    Select templates for multi-template modelling according to the algorithm of
    `Meier and Söding <https://doi.org/10.1371/journal.pcbi.1004343>`_.

    The algorithm attempts to select a list of templates in such a way that
    optimises a combination of coverage and restraint quality. The estimated
    alignment quality is used as a proxy for restraint quality.

    .. warning::

        First, each template is sorted according to some quality metric. This
        component assumes that the list of templates has already been sorted.
        If you pass an unsorted list of templates to this component, you will
        get nonsense results.

    Each residue of each template is assigned a *local quality score*
    :math:`s(i,t)`, defined as the product of the probability
    :math:`P_\text{hom}(t)` that template :math:`t` is homologous to query
    :math:`q` and the probability :math:`p(i,j|q,t)` that residue :math:`i` from
    :math:`q` is correctly aligned to residue :math:`j` of :math:`t`:

    .. math::

        s(i, t) = P_\text{hom}(t) p(i, j|q, t).

    Then, a change in score :math:`\Delta s(i, t)` is calculated for each
    residue of :math:`t` by comparison to the highest-scoring residue pair at
    position :math:`i` in the list of accepted templates :math:`T_\text{acc}`:

    .. math::

        \Delta s(i, t) = s(i, t) - \max_{t' \in T_\text{acc}} s(i, t').

    The change in score :math:`\Delta s(i, t)` represents the increase or
    decrease in alignment quality at query residue :math:`i` resulting from the
    inclusion of template :math:`t` in the :math:`T_\text{acc}`. The list of
    accepted templates is initially populated with the top-ranked template.

    To weight the positive values of :math:`\Delta s` more strongly than
    negative values, the exponential function is applied. A global score
    :math:`S(t)` for template :math:`t` is then calculated by summing over all
    local scores.

    .. math::

        S(t) = \sum_{i} \left[
            \exp\left(\alpha \Delta s(i, t)\right) - \beta
        \right].

    The constant factors :math:`\alpha` and :math:`\beta` were found by Meier
    and Söding by a grid search. We use the values :math:`\alpha = 0.95` and
    :math:`\beta = 1` from their work. These can be adjusted via the instance
    constants :py:attr:`.alpha` and :py:attr:`.beta`.

    This component will add the field ``template_at_residue`` to the pipeline
    state. This is a list with the same length as the ``sequence`` field,
    each element of which is a tuple containing the maximum assigned residue
    score and the winning template for that position (or `None` if no template
    was chosen for that position).
    """
    REQUIRED = ["sequence", "templates"]
    ADDS = ["template_at_residue"]
    REMOVES = []

    def __init__(self, alpha=0.95, beta=1.00):
        self.alpha = alpha
        self.beta = beta

    def run(self, data, config=None, pipeline=None):
        """Picking models according to Meier and Söding's algorithm."""
        sequence, templates = self.get_vals(data)

        template_at_residue = [None] * len(sequence)
        accepted_templates = []
        top_confidences = np.zeros(len(sequence))

        for template in templates:
            alignment = template["alignment"]

            # Indices into the query of the alignment
            query_indices = np.array([pair.i - 1 for pair in alignment])

            # If there are no residue pairs in the template, skip it
            if len(query_indices) == 0:
                continue

            prob_homologous = template["prob"] / 100
            # Defined over each residue in the index
            local_scores = (np.array([pair.probab for pair in alignment])
                            * prob_homologous)
            Δ_local_scores = local_scores - top_confidences[query_indices]

            global_score = np.sum(np.exp(self.alpha * Δ_local_scores)
                                  - self.beta)
            if global_score > 0:
                accepted_templates.append(template)

            # Defined over the residues in the alignment
            update_required = local_scores > top_confidences[query_indices]
            # Subset of indices to update defined on the query
            update_indices = query_indices[update_required]
            top_confidences[update_indices] = local_scores[update_required]
            for i in update_indices:
                template_at_residue[i] = (top_confidences[i], template)

        data["templates"] = accepted_templates
        data["template_at_residue"] = template_at_residue
        return data

class LoopModel(Component):
    """
    Use Alex Herbert's loop modeler to fill in as many gaps as possible.

    :param str bin_dir: Location of the loop modelling executable.
    :param str config: Loop modeller configuration file.
    :param str executable: Name of the executable to run, under `bin_dir`.
    """

    # TODO: Write a component to transfer top-level keys to lower levels so
    # we can run this component with component.Map instead of explicitly
    # looping over ``templates.``

    REQUIRED = ["pssm", "templates"]
    ADDS = []
    REMOVES = []

    LOOP_MODELLER = ExternalTool({
        "config": "c",
        "query": "f",
        "out_dir": "d",
        "model_list": "l",
    }, long_prefix="-")

    def __init__(self, bin_dir, config, executable="assembler.loop"):
        self.bin_dir = bin_dir
        self.config = config
        self.executable = executable

    def convert_ascii_pssm(self, pssm, output):
        """
        Convert the PSSM generated by blastpgp into the loop modeller format.

        This means chopping off the header and the column with residue IDs.
        """
        with open(pssm, "r") as pssm_in:
            for line in pssm_in:
                if re.match(r"^\s+\d", line):
                    cols = line.split()
                    residue = cols[1]
                    counts = cols[2:22]
                    output.write(residue + " ")
                    for count in counts:
                        output.write("{:4d}".format(int(count)))
                    output.write("\n")
        output.flush()

    def run(self, data, config=None, pipeline=None):
        """Fill short gaps with loop modeller."""
        pssm, templates = self.get_vals(data)

        try:
            tmpdir = tempfile.mkdtemp("-loop", "phyreengine-")
            self.logger.debug("Loop modelling using tmpdir: %s", tmpdir)

            out_dir = Path("loop")

            # Attempt to use existing models if the output directory exists.
            if not out_dir.exists():
                loop_pssm = Path(tmpdir, "loop.pssm")
                query_fasta = Path(tmpdir, "query.fasta")
                model_list = Path(tmpdir, "model.list")

                with model_list.open("w") as model_list_out:
                    for template in templates:
                        print(str(Path(template["model"]).resolve()),
                              file=model_list_out)
                with loop_pssm.open("w") as loop_out:
                    self.convert_ascii_pssm(pssm["ascii"], loop_out)
                with query_fasta.open("w") as query_out:
                    print(">{name}\n{sequence}\n".format(**data),
                          file=query_out)

                command_line = self.LOOP_MODELLER(
                    executable=(self.bin_dir, self.executable),
                    options={
                        "config": self.config,
                        "pssm": loop_pssm,
                        "query": query_fasta,
                        "model_list": model_list,
                        "out_dir": out_dir,
                    })
                self.logger.debug("Running %s", command_line)
                subprocess.run(command_line, check=True)

            # Replace "model" field of each template with the first loop model.
            for i, template in enumerate(templates):
                model_idx = i + 1
                model_path = (out_dir
                              / "model.{}".format(model_idx)
                              / "model.1.pdb")
                if not model_path.exists():
                    err_msg = "Loop-modelled file '{:s}' does not exist"
                    raise FileNotFoundError(err_msg.format(model_path))
                template["model"] = str(model_path)
        finally:
            shutil.rmtree(tmpdir)
        return data
