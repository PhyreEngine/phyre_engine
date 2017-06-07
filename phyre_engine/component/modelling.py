"""Module providing utilities for homology modelling."""

import Bio.SeqUtils
import Bio.PDB.Structure
import Bio.PDB.Model
import Bio.PDB.Chain
from Bio.PDB.Residue import Residue
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from phyre_engine.component import Component
import os.path

BACKBONE_ATOMS = {"N", "C", "CA", "O"} # Set of backbone atoms


class HomologyModeller(Component):
    """For each hit, generate a model.

    For the moment, we will do the parsing of hhsearch files from here. In
    actuality, we will have a separate component to parse those files so this
    can be more generic.
    """

    REQUIRED = ['sequence', 'hits']
    ADDS     = ['models']
    REMOVES  = []

    def __init__(self, mmcif_dir):
        """
        Create a new modeller.

        :param str mmcif_dir: Directory containing MMCIF files.
        """
        self.mmcif_dir = mmcif_dir

    def run(self, data):
        """
        Run the modeller.
        """

        query_seq, hits = self.get_vals(data)

        models = []
        for hit in hits:
            model_structure = self._map_aln_to_struc(query_seq, hit)
            models.append(model_structure)

        data["models"] = models
        return data

    def _parse_pdb_id(self, hit):
        """Return a tuple indicating the PDB ID and chain."""

        #The first word (i.e. before the first space char) of the name is the 
        #PDB code + chain. The PDB code and chain are separated by an
        #underscore.
        full_id = hit.info["name"].split(maxsplit=1)[0]
        return tuple(full_id.split("_"))

    def _find_mmcif_file(self, pdb_id):
        """Finds the MMCIF file corresponding to the ID given in pdb_id.

        We assume that the directory structure and names are identical to those
        on the RCSB's FTP server. That is, structure ``12as`` would be in
        ``<root>/2a/12as.cif``.

        :param str pdb_id: PDB code.
        """

        middle = pdb_id[1:3].lower()
        file   = "{}.cif".format(pdb_id.lower())
        path   = os.path.join(self.mmcif_dir, middle, file)
        return path


    def _map_aln_to_struc(self, query_seq, hit):
        """Map a sequence alignment to a structure.

        This method assumes that the sequence alignment is based on the
        *sequence* of a structure file; that is, the ``SEQRES`` records  or the
        ``label_seq_id``.
        """

        mmcif_parser = MMCIFParser()

        pdb_code, pdb_chain = self._parse_pdb_id(hit)
        mmcif_file = self._find_mmcif_file(pdb_code)
        structure = mmcif_parser.get_structure(hit.info["name"], mmcif_file)
        chain = structure[0][pdb_chain]

        # Init new structure
        model_structure = Bio.PDB.Structure.Structure(
                "model from {}_{}".format(pdb_code, pdb_chain))
        model_model = Bio.PDB.Model.Model(1)
        model_chain = Bio.PDB.Chain.Chain("A")
        model_model.add(model_chain)
        model_structure.add(model_model)

        # Get mapping of simple 1-based index to full atom IDs
        label_to_auth = self._auth_to_label_ids(mmcif_file, structure,
                pdb_chain)

        for aln_pos in hit.aln:
            #indices i and j are the query and template indices
            #respectively. We want to create a model residue that is
            #identical to the template residue but has different identifiers
            #and residue type.
            i = aln_pos["i"]
            j = aln_pos["j"]

            # The residue may be missing from the template structure, in which
            # case we just skip it. The default param for get() will be used if
            # auth_seq_id was not defined in the mmcif file, and is just the
            # label_seq_id.
            auth_asym_id = label_to_auth.get(j, j)
            if not auth_asym_id in chain:
                continue
            template_res = chain[auth_asym_id]

            #Remember that the sequence object is basically just a string,
            #and so is zero-indexed. Structure indices are based on
            #label_seq_id, so are one-indexed.
            query_res_type = Bio.SeqUtils.seq3(query_seq[i-1]).upper()
            query_res = Residue((" ", i, " "), query_res_type, " ")
            for atom in template_res:
                if atom.get_name() in BACKBONE_ATOMS:
                    query_res.add(atom.copy())
            model_chain.add(query_res)
        return model_structure


    def _auth_to_label_ids(self, mmcif_file, structure, pdb_chain):
        """
        Build a mapping of ``label_seq_id`` fields to ``auth_seq_id`` fields.

        Sequence alignments do not have any idea what identifier the author of
        a structure assigned to which residue, and are most sensibly indexed
        from 1. The author-assigned fields of structures may be ordered
        completely arbitrarily, and so we would prefer to use the
        ``label_seq_id`` field when working with sequence alignments. Biopython
        prefers author-assigned IDs, so this method allows us to build a
        mapping between a label ID and author ID.

        :param str mmcif_file: File from which to parse mappings.
        :param `Bio.PDB.Structure` structure: Structure object.
        :param str pdb_chain: The chain of interest.

        :return:
            Dictionary mapping ``label_seq_id`` to the full ID tuples of
            ``(group, auth_seq_id, pdbx_PDB_ins_code)`` used by BioPython to
            index residues.
        """

        # Sequence alignments index from 1, because they don't care what the
        # sequence identifiers the authors of a PDB file assigned to their
        # residues. Bio.PDB uses the author-assigned identifiers if they are
        # available, but the label_seq_id identifier would be much more
        # suitable. We build a map of label_seq_id IDs to auth_seq_ids so that
        # we can find a residue by its label_seq_id.
        #
        # If the structure does not contain an auth_seq_id then Bio.PDB will
        # use the label_seq_id, which is mandatory, so we will just build an
        # identity dict in that case.

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
            generator = zip(labels, auths, chains, model_nums, group, icodes)
            for res_id in generator:
                if (res_id[2] == pdb_chain
                    and res_id[3] == "1"
                    and res_id[4] == "ATOM"):

                    icode = res_id[5] if res_id[5] != "?" else " "
                    label_to_auth[int(res_id[0])] = (' ', int(res_id[1]), icode)
        return label_to_auth
