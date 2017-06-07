"""Module providing utilities for converting MMCIF files to PDB."""

import Bio.PDB
import os

class SimplifySelector(Bio.PDB.Select):
    """Class to select what should be output.

    The atoms to be output must satisfy the following constraints:

    * Part of an ``ATOM`` residue, not ``HETATM``.
    * One of the standard 20 amino acids.
    * If atoms have altlocs, the first location is chosen.

    Note that this class does *not* filter by model or chain. This is because
    we cannot be certain what the first model ID is in an MMCIF file. To ensure
    you only get a single model, pass only that model to PDBIO. To get only a
    single chain, pass that chain to PDBIO.
    """

    def accept_residue(self, res):
        """Keep ATOM residues only."""
        return res.get_id()[0] == ' '

    def accept_atom(self, atom):
        """Keep non-disordered atoms or the first encountered."""
        return (not atom.is_disordered()) or atom.get_altloc() == "A"

class MMCIFToPDBChains:
    """Used to build a valid PDB file for each chain of an MMCIF file.

    :param out_dir: Directory to which PDB files will be saved.
    """

    def __init__(self, out_dir):
        """Initialise a new converter."""
        self.out_dir = out_dir

    def convert(self, pdb_id, mmcif_file):
        """Convert an MMCIF file to many PDB files.

        Each chain in the MMCIF file will be given its own PDB file. The chain
        identifiers in MMCIF can contain more than one character, so the chain
        ID in each PDB file will be set to a single space character.

        Files with the PDB ID ``1xyz`` will be placed in the subdirectory
        ``xy``. This will always be converted to lowercase. Output files will
        be named like ``1xyz_AA.pdb``, where ``AA`` is the chain ID from which
        this file was generated.

        :param str pdb_id: ID of the structure for naming purposes.
        :param str mmcif_file: Path to an MMCIF file from which chains will be
            extracted.
        """

        pdbio = Bio.PDB.PDBIO()
        simple_selector = SimplifySelector()
        structure = Bio.PDB.MMCIFParser().get_structure(mmcif_file, mmcif_file)
        model = next(structure.get_models())

        base_dir = os.path.join(self.out_dir, pdb_id[1:3].lower())
        os.makedirs(base_dir, exist_ok=True)

        # There is a bug (or least behaviour that I find non-obvious) in
        # biopython. When calling pdbio.set_structure, the parent of the chain
        # is altered, which means that resetting the chain_id does not work and
        # biopython complains about duplicate chain IDs when setting the next
        # chain ID to " ". I'm going to hack around this bug by iterating over
        # the list of chains and deleting chains we've already seen from the
        # MMCIf model.

        chain_list = list(model)
        for chain in chain_list:
            file_name = "{}_{}.pdb".format(pdb_id.lower(), chain.get_id())
            out_path = os.path.join(base_dir, file_name)

            # Temporarily set chain ID to " "
            chain_id = chain.get_id()
            chain.id = " "
            pdbio.set_structure(chain)
            pdbio.save(out_path, simple_selector)
            # Delete chain from model as part of bugfix/hack mentioned above.
            del model[" "]
