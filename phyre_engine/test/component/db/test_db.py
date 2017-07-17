import json
import tempfile
import unittest
import phyre_engine.component.db.db as db
import phyre_engine.test
import phyre_engine.test.data
import phyre_engine.test.tools.test_pdb
from pathlib import Path
import shutil
import Bio.SeqIO

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACProtein

class TestStructureRetriever(unittest.TestCase):
    """Test STructureRetriever component."""
    _PIPE_STATE = {"templates": [
        {"PDB": "12as"},
        {"PDB": "4HHB"}
    ]}

    def test_retrieve(self):
        """Try and download some files."""
        types = ("pdb", "cif", db.StructureType.PDB, db.StructureType.MMCIF)

        for struc_type in types:
            with self.subTest("Getting structure", type=struc_type):
                with tempfile.TemporaryDirectory() as tmpdir:
                    retriever = db.StructureRetriever(struc_type, tmpdir)
                    retriever.run(self._PIPE_STATE)

                    if isinstance(struc_type, str):
                        suffix = struc_type
                    else:
                        suffix = struc_type.value

                    pdb_12as = Path(tmpdir, "2a/12as.{}.gz".format(suffix))
                    pdb_4hhb = Path(tmpdir, "hh/4hhb.{}.gz".format(suffix))
                    self.assertTrue(
                        pdb_12as.exists(),
                        "{!s} should exist".format(pdb_12as))
                    self.assertTrue(
                        pdb_4hhb.exists(),
                        "{!s} should exist".format(pdb_4hhb))


class TestChainPDBBuilder(unittest.TestCase):
    """Test ChainPDBBuilder class"""

    _12ASA_SEQ = (
        "AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
        "AVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDE"
        "DRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSE"
        "EFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGG"
        "KLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMG"
        "IRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTML"
        "LLQLPHIGQVQAGVWPAAVRESVPSLLN")

    def setUp(self):
        """Create a temporary directory for tests to use."""
        self.tmpdir = Path(tempfile.mkdtemp("-test", "chain-"))
        self.mmcif_dir = Path(phyre_engine.test.data.__file__).parent / "mmcif"

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(str(self.tmpdir))

    def test_build(self):
        """
        Try and extract a chain. Test that the extracted chain exists and that
        it contains REMARK 149 and ATOM lines matching the required sequence.
        """
        pdb_dir = self.tmpdir

        builder = db.ChainPDBBuilder(str(self.mmcif_dir), str(pdb_dir))
        results = builder.run({"templates": [{"PDB": "12as", "chain": "A"}]})

        pdb_12as_A_path = Path(results["templates"][0]["structure"])
        self.assertTrue(
            pdb_12as_A_path.exists(),
            "PDB file containing chain was created.")

        with pdb_12as_A_path.open("r") as pdb_in:
            self.assertEqual(
                str(Bio.SeqIO.read(pdb_in, "pdb-atom").seq), self._12ASA_SEQ,
                "Atom sequence matches read sequence")

        with pdb_12as_A_path.open("r") as pdb_in:
            json_str = ""
            for line in pdb_in:
                if not line.startswith("REMARK 149"):
                    continue
                json_str += line.replace("REMARK 149 ", "").rstrip()
            mapping = json.loads(json_str)
            map_dict = {m[0]: tuple(m[1]) for m in mapping}
            self.assertEqual(map_dict[1], (' ', 4, ' '))


    def test_build_with_selectors(self):
        """Write a structure with an active conformation selector."""
        class FilterHetatms:
            def select(self, chain):
                for res in list(chain.get_residues()):
                    if res.get_id()[0] != ' ':
                        del chain[res.get_id()]
                return chain

        pdb_dir = self.tmpdir

        builder = db.ChainPDBBuilder(
            str(self.mmcif_dir), str(pdb_dir),
            conf_sel=[FilterHetatms()])
        results = builder.run({"templates": [{"PDB": "12as", "chain": "A"}]})

        with Path(results["templates"][0]["structure"]).open("r") as pdb_in:
            # Results should match the same sequence without the terminating
            # residue, which is marked as a HETATM.
            self.assertEqual(
                str(Bio.SeqIO.read(pdb_in, "pdb-atom").seq),
                self._12ASA_SEQ[0:-1],
                "Atom sequence matches read sequence")

class TestPDBSequence(unittest.TestCase):
    """Tests for the PDBSequence component."""

    def test_sequence(self):
        """Sequence should match atom_seq function."""

        with tempfile.NamedTemporaryFile("w") as pdb_fh:
            pdb_fh.write(phyre_engine.test.tools.test_pdb.ATOM_SEQ_TEST_PDB)
            pdb_fh.flush()
            seq_parser = db.PDBSequence()
            results = seq_parser.run({
                "templates": [{"structure": pdb_fh.name}]
            })
            self.assertEqual(
                str(results["templates"][0]["sequence"].seq), "AG",
                "Sequence read correctly")

class TestReduce(unittest.TestCase):
    """Test the Reduce component."""

    @staticmethod
    def template_name(template):
        return (template["PDB"], template["chain"])

    _TEMPLATES = [
        {"sequence": SeqRecord(Seq("AGH", IUPACProtein)),
         "PDB": "1foo", "chain": "X"},
        {"sequence": SeqRecord(Seq("AGH", IUPACProtein)),
         "PDB": "1foo", "chain": "Y"},
        {"sequence": SeqRecord(Seq("AGH", IUPACProtein)),
         "PDB": "1bar", "chain": "A"},
        {"sequence": SeqRecord(Seq("HHH", IUPACProtein)),
         "PDB": "1baz", "chain": "B"},
        {"sequence": SeqRecord(Seq("HHH", IUPACProtein)),
         "PDB": "1baz", "chain": "C"},
        {"sequence": SeqRecord(Seq("HH", IUPACProtein)),
         "PDB": "1qux", "chain": "A"}
    ]

    def test_reduce(self):
        """Identical sequences should be combined."""
        reduce = db.Reduce()
        results = reduce.run({"templates": self._TEMPLATES})
        self.assertEqual(len(results["templates"]), 3, "Reduced to 3 seqs")

        self.assertSetEqual(
            set([str(t["sequence"].seq) for t in results["templates"]]),
            set(["AGH", "HHH", "HH"]),
            "New templates contains all unique seqs")

        self.assertEqual(len(results["reduction"]), 3, "Reduced to 3 seqs")
        self.assertSetEqual(
            set(results["reduction"].keys()),
            set(["AGH", "HHH", "HH"]),
            "reduction key is indexed by all unique seqs")

        self.assertSetEqual(
            set([self.template_name(t) for t in results["reduction"]["AGH"]]),
            set([self.template_name(t) for t in self._TEMPLATES[0:3]]),
            "Templates in results['reduction']['AGH'] correct")
        self.assertSetEqual(
            set([self.template_name(t) for t in results["reduction"]["HHH"]]),
            set([self.template_name(t) for t in self._TEMPLATES[3:5]]),
            "Templates in results['reduction']['HHH'] correct")
        self.assertSetEqual(
            set([self.template_name(t) for t in results["reduction"]["HH"]]),
            set([self.template_name(t) for t in self._TEMPLATES[5:]]),
            "Templates in results['reduction']['HH'] correct")

    def test_json(self):
        """Test JSON writer."""

        with tempfile.NamedTemporaryFile("w+") as json_file:
            reduce = db.Reduce(json_file.name)
            reduce.run({"templates": self._TEMPLATES})
            json_file.seek(0)
            json_content = json.load(json_file)

            self.assertIsInstance(json_content, dict, "Read map from JSON")
            self.assertEqual(len(json_content), 3, "Read 3 groups")

            self.assertSetEqual(
                set([tuple(t) for t in json_content["AGH"]]),
                set([self.template_name(t) for t in self._TEMPLATES[0:3]]),
                "AGH grouping correct")

            self.assertSetEqual(
                set([tuple(t) for t in json_content["HHH"]]),
                set([self.template_name(t) for t in self._TEMPLATES[3:5]]),
                "HHH grouping correct")

            self.assertSetEqual(
                set([tuple(t) for t in json_content["HH"]]),
                set([self.template_name(t) for t in self._TEMPLATES[5:]]),
                "H grouping correct")

class TestExpand(unittest.TestCase):
    """Expand templates list using identical templates."""

    _AGH_SEQ = SeqRecord(Seq("AGH"))
    _HHH_SEQ = SeqRecord(Seq("HHH"))
    _TEST_REDUCTION = {
        "AGH": [
            {"sequence": _AGH_SEQ,
             "PDB": "1foo", "chain": "A", "extra_var": "foo"},
            {"sequence": _AGH_SEQ,
             "PDB": "1foo", "chain": "B", "extra_var": "foo"},
        ]
    }
    _TEST_JSON_REDUCTION = {
        "AGH": [["1foo", "A"], ["1foo", "B"]]
    }
    _TEST_TEMPLATES = [
        {"sequence": _AGH_SEQ,
         "PDB": "1baz", "chain": "A", "another": "var"},
        {"sequence": _HHH_SEQ}
    ]
    _TEST_STATE = {
        "templates": _TEST_TEMPLATES,
        "reduction": _TEST_REDUCTION
    }

    @staticmethod
    def _convert_seqrecords(templates):
        # We can't compare SeqRecords directly, so use this to replace them
        # with their sequences.
        converted = []
        for template in templates:
            copy = template.copy()
            copy["sequence"] = str(copy["sequence"].seq)
            converted.append(copy)
        return converted

    def test_expand_reduction(self):
        """Expand using the pipeline 'reduction' key."""
        expand = db.Expand()
        state = {
            "templates": self._TEST_TEMPLATES.copy(),
            "reduction": self._TEST_REDUCTION.copy()
        }
        results = expand.run(state)
        templates = self._convert_seqrecords(results["templates"])

        self.assertEqual(len(templates), 4, "Added two sequences")
        self._test_original_remain(templates)

        self.assertIn(
            {"sequence": "AGH", "PDB": "1foo", "chain": "A",
             "extra_var": "foo", "another": "var"},
            templates,
            "Added 1foo_A including 'extra_var' and 'another'")

        self.assertIn(
            {"sequence": "AGH", "PDB": "1foo", "chain": "B",
             "extra_var": "foo", "another": "var"},
            templates,
            "Added 1foo_B including 'extra_var' and 'another'")


    def test_expand_json(self):
        """Test expand using the reduction_file parameter."""
        with tempfile.NamedTemporaryFile("w") as json_file:
            json.dump(self._TEST_JSON_REDUCTION, json_file)
            json_file.flush()
            json_file.seek(0)

            expand = db.Expand(json_file.name)
            state = {
                "templates": self._TEST_TEMPLATES.copy()
            }
            results = expand.run(state)
            templates = self._convert_seqrecords(results["templates"])

            self.assertEqual(len(templates), 4, "Added two sequences")
            self._test_original_remain(templates)

            self.assertIn(
                {"sequence": "AGH", "PDB": "1foo", "chain": "A",
                 "another": "var"},
                templates,
                "Added 1foo_A including 'another'")

            self.assertIn(
                {"sequence": "AGH", "PDB": "1foo", "chain": "B",
                 "another": "var"},
                templates,
                "Added 1foo_B including 'another'")


    def test_param_combination(self):
        """Expand requires ONE of reduction_file or reduction key."""
        with self.assertRaises(ValueError):
            db.Expand().run({"templates": self._TEST_TEMPLATES.copy()})

        with self.assertRaises(ValueError):
            db.Expand("foo").run({
                "templates": self._TEST_TEMPLATES.copy(),
                "reduction": self._TEST_REDUCTION.copy()
            })

    def _test_original_remain(self, results):
        self.assertIn(
            {"sequence": "AGH", "PDB": "1baz", "chain": "A", "another": "var"},
            results,
            "Retained 1baz_A")

        self.assertIn(
            {"sequence": "HHH"},
            results,
            "Retained template with no PDB or chain")

if __name__ == "__main__":
    unittest.main()
