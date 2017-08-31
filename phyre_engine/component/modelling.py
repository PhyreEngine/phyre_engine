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
from phyre_engine.component.db.db import ChainPDBBuilder
from phyre_engine.tools.template import Template

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

            template_file = pdb.pdb_path(pdb_id, ".pdb", chain, self.chain_dir)
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

            model_file = self.model_name.format(**template)
            pdb_io.set_structure(model_structure)
            pdb_io.save(model_file)
            template["model"] = model_file
        return data
