"""Test components in :py:mod:`phyre_engine.component.foldlib`."""
import contextlib
import datetime
import pathlib
import sqlite3
import tempfile
import unittest
import unittest.mock

import phyre_engine.tools.template as ttools
import phyre_engine.component.foldlib as foldlib

SAMPLE_METADATA = {
    "deposition_date": datetime.date(1997, 12, 2),
    "release_date": datetime.date(1998, 12, 30),
    "last_update_date": datetime.date(2011, 7, 13),
    "method": "X-RAY DIFFRACTION",
    "resolution": 2.2,
    "organism_name": "Escherichia coli",
    "organism_id": 562,
    "title": (
        "ASPARAGINE SYNTHETASE MUTANT C51A, C315A COMPLEXED "
        "WITH L-ASPARAGINE AND AMP"),
    "descriptor": (
        "ASPARAGINE SYNTHETASE, L-ASPARAGINE, "
        "ADENOSINE MONOPHOSPHATE"),
}


@contextlib.contextmanager
def TemporaryDatabase():
    with tempfile.TemporaryDirectory("-template-db", "phyreengine-") as tmp:
        tmp = pathlib.Path(tmp)

        sql_db = tmp / "foldlib.db"
        chain_dir = tmp / "chains"
        ttools.TemplateDatabase.create(tmp / "foldlib.db")
        chain_dir.mkdir()
        yield ttools.TemplateDatabase(sql_db, chain_dir)

class TestOpenCopy(unittest.TestCase):
    """Test making copy of existing SQL database."""

    def test_copy(self):
        """Copy does not affect master."""
        with TemporaryDatabase() as tdb:
            copy_db = foldlib.OpenCopy(tdb.database, tdb.file_root)
            results = copy_db.run({})
            results["template_db"].add_pdb("12AS", SAMPLE_METADATA)
            results["template_db"].commit()

            expected_exception = ttools.TemplateDatabase.PdbNotFoundException
            with self.assertRaises(expected_exception):
                tdb.get_pdb("12AS")

    def test_dump(self):
        """SQL statements made on copy are dumped."""
        with TemporaryDatabase() as tdb:
            copy_db = foldlib.OpenCopy(tdb.database, tdb.file_root)
            results = copy_db.run({})
            results["template_db"].add_pdb("12AS", SAMPLE_METADATA)
            results["template_db"].commit()

            self.assertGreater(len(results["sql_dump"]), 1)
            self.assertTrue(any([sql.strip().startswith("INSERT INTO pdbs")
                            for sql in results["sql_dump"]]))


class TestRestoreSQLDump(unittest.TestCase):
    """Test RestoreSQLDump component."""

    STATEMENT = """\
        INSERT INTO pdbs (
           pdb_id, deposition_date, last_update_date, release_date,
           method, resolution,
           organism_name, organism_id,
           title, descriptor)
        VALUES (
           LOWER('12as'), '1997-12-02', '2011-07-13', '1998-12-30',
           'X-RAY DIFFRACTION', 2.2, 'Escherichia coli', 562,
           'ASPARAGINE SYNTHETASE MUTANT C51A, C315A COMPLEXED WITH L-ASPARAGINE AND AMP',
           'ASPARAGINE SYNTHETASE, L-ASPARAGINE, ADENOSINE MONOPHOSPHATE')
    """

    def test_restore(self):
        """Restore template DB from SQL statements."""
        with TemporaryDatabase() as db:
            restore = foldlib.RestoreSQLDump()
            results = restore.run({
                "sql_dump": [self.STATEMENT],
                "template_db": db,
            })
            self.assertEqual(
                results["template_db"].get_pdb("12AS")["organism_id"],
                562)
