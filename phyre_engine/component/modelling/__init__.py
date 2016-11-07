"""Module providing utilities for homology modelling."""

import Bio.SeqUtils
import Bio.PDB.Structure
import Bio.PDB.Model
import Bio.PDB.Chain
from Bio.PDB.Residue import Residue
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from phyre_engine.component import Component
import phyre_engine.tools.hhsuite.parser as parser
import os.path

class HomologyModeller(Component):
    """For each hit, generate a model.

    For the moment, we will do the parsing of hhsearch files from here. In
    actuality, we will have a separate component to parse those files so this
    can be more generic.
    """

    REQUIRED = ['sequence', 'hhsearch_report', 'hhsearch_atab']
    ADDS     = ['models']
    REMOVES  = []

    def __init__(self, mmcif_dir):
        """
        Create a new modeller.

        Args:
            ``mmcif_dir``: Directory containing MMCIF files.
        """
        self.mmcif_dir = mmcif_dir

    def run(self, data):
        """
        Run the modeller.
        """

        query_seq, report_file, atab_file = self.get_vals(data)

        #Parse report files to get hits, then combine the hits so that we have
        #the summary information and alignments in the same hit.
        report = parser.Report(report_file)
        atab   = parser.Tabular(atab_file)
        hits = []
        for r, a in zip(report.hits, atab.hits):
            hit = parser.Hit()
            hit.info = r.info
            hit.aln  = a.aln
            hits.append(hit)

        mmcif_parser = MMCIFParser()
        models = []
        for hit in hits:
            pdb_code, pdb_chain = self._parse_pdb_id(hit)
            mmcif_file = self._find_mmcif_file(pdb_code)
            structure = mmcif_parser.get_structure(hit.info["name"], mmcif_file)
            chain = structure[0][pdb_chain]

            # Sequence alignments index from 1, because they don't care what the
            # sequence identifiers the authors of a PDB file assigned to their
            # residues. Bio.PDB uses the author-assigned identifiers if they are
            # available, but the label_seq_id identifier would be much more
            # suitable. We build a map of label_seq_id IDs to auth_seq_ids so that
            # we can find a residue by its label_seq_id. If the structure does
            # not contain an auth_seq_id then Bio.PDB will use the label_seq_id,
            # which is mandatory, so we will just build an identity dict in that
            # case.
            mmcif_dict = MMCIF2Dict(mmcif_file)
            label_to_auth = {}
            if "_atom_site.auth_seq_id" in mmcif_dict:
                labels = mmcif_dict["_atom_site.label_seq_id"]
                auths  = mmcif_dict["_atom_site.auth_seq_id"]
                chains = mmcif_dict["_atom_site.auth_asym_id"]
                group  = mmcif_dict["_atom_site.group_PDB"]
                icodes = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]

                #Get model numbers of available or just use a default value of
                model_nums = None
                if "_atom_site.pdbx_PDB_model_num" in mmcif_dict:
                    model_nums = mmcif_dict["_atom_site.pdbx_PDB_model_num"]
                else:
                    model_nums = ["1"] * len(labels)

                #Filter to only contain ATOM records from the first model with
                #the desired chain. Store the icode as well as the auth_seq_id
                for id in zip(labels, auths, chains, model_nums, group, icodes):
                    if id[2] == pdb_chain and id[3] == "1" and id[4] == "ATOM":
                        icode = id[5] if id[5] != "?" else " "
                        label_to_auth[id[0]] = (' ', int(id[1]), icode)

            model_structure = Bio.PDB.Structure.Structure(
                    "model from {}_{}".format(pdb_code, pdb_chain))
            model_model = Bio.PDB.Model.Model(1)
            model_chain = Bio.PDB.Chain.Chain("A")
            model_model.add(model_chain)
            model_structure.add(model_model)

            for aln_pos in hit.aln:
                #indices i and j are the query and template indices
                #respectively. We want to create a model residue that is
                #identical to the template residue but has different identifiers
                #and residue type.
                i = aln_pos["i"]
                j = aln_pos["j"]

                #The residue may be missing from the template structure, in
                #which case we just skip it.
                auth_asym_id = label_to_auth.get(str(j), str(j))
                if not auth_asym_id in chain:
                    continue
                template_res = chain[auth_asym_id]

                #Remember that the sequence object is basically just a string,
                #and so is zero-indexed. Structure indices are based on
                #label_seq_id, so are one-indexed.
                query_res_type = Bio.SeqUtils.seq3(query_seq[i-1]).upper()
                query_res = Residue((" ", i, " "), query_res_type, " ")
                backbone_atoms = {"N", "C", "CA", "O"}
                for atom in template_res:
                    if atom.get_name() in backbone_atoms:
                        query_res.add(atom.copy())
                model_chain.add(query_res)
            models.append(model_structure)
        data["models"] = models
        return data

    def _parse_pdb_id(self, hit):
        """Return a tuple indicating the PDB ID and chain."""

        #The first word (i.e. before the first space char) of the name is the 
        #PDB code + chain. The PDB code and chain are separated by an
        #underscore.
        full_id, _ = hit.info["name"].split(maxsplit=1)
        return tuple(full_id.split("_"))

    def _find_mmcif_file(self, pdb_id):
        """Finds the MMCIF file corresponding to the ID given in pdb_id.

        We assume that the directory structure and names are identical to those
        on the RCSB's FTP server. That is, structure ``12as`` would be in
        ``<root>/2a/12as.cif``.

        Args:
            ``pdb_id``: PDB code.
        """

        middle = pdb_id[1:3].lower()
        file   = "{}.cif".format(pdb_id.lower())
        path   = os.path.join(self.mmcif_dir, middle, file)
        return path
