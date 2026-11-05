#!/usr/bin/env python3
"""Synthetic transactional tests; no SGE, MAJIQ or real DB required.

Run: python3 -m unittest discover -s reads_processing/tests -p 'test_incremental_sqlite_writer.py'
"""
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from incremental_sqlite_writer import commit_microbatch

SCHEMA = """
PRAGMA foreign_keys=ON;
CREATE TABLE genome_info(genome_version TEXT PRIMARY KEY, species_name TEXT);
CREATE TABLE splicing_events(
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 junctionID TEXT UNIQUE, eventID TEXT, geneID TEXT, seqID TEXT, strand TEXT,
 event_type TEXT, description TEXT, start INTEGER, end INTEGER, coord TEXT,
 full_coord TEXT, upstream_exon TEXT, downstream_exon TEXT, species TEXT,
 genome_version TEXT,
 FOREIGN KEY(genome_version) REFERENCES genome_info(genome_version));
CREATE TABLE junction_sample_stats(
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 junctionID TEXT, psi_majiq REAL, psi_sgseq REAL, reads_majiq REAL,
 reads_sgseq_3p REAL, reads_sgseq_5p REAL, sra_id TEXT, species TEXT,
 UNIQUE(junctionID,sra_id),
 FOREIGN KEY(junctionID) REFERENCES splicing_events(junctionID));
CREATE TABLE junction_summary_stats(
 id INTEGER PRIMARY KEY AUTOINCREMENT, junctionID TEXT UNIQUE NOT NULL,
 mean_psi_majiq REAL DEFAULT 0, mean_psi_sgseq REAL DEFAULT 0,
 sum_reads_majiq REAL DEFAULT 0, sum_reads_sgseq_3p REAL DEFAULT 0,
 sum_reads_sgseq_5p REAL DEFAULT 0, n_of_samples_majiq INTEGER DEFAULT 0,
 n_of_samples_sgseq INTEGER DEFAULT 0,
 FOREIGN KEY(junctionID) REFERENCES splicing_events(junctionID));
CREATE TABLE sample_processing_status(
 sample_id TEXT NOT NULL, genome_version TEXT NOT NULL,
 status TEXT NOT NULL CHECK(status IN ('SUCCESS','FAILED')),
 failed_stage TEXT NOT NULL DEFAULT '', exit_code TEXT NOT NULL DEFAULT '',
 reason TEXT NOT NULL DEFAULT '', updated_at TEXT DEFAULT CURRENT_TIMESTAMP,
 PRIMARY KEY(sample_id,genome_version),
 FOREIGN KEY(genome_version) REFERENCES genome_info(genome_version));
CREATE TABLE sample_qc_metrics(
 sra_id TEXT NOT NULL, genome_version TEXT NOT NULL, tool TEXT NOT NULL,
 metric TEXT NOT NULL, numeric_value REAL, text_value TEXT,
 PRIMARY KEY(sra_id,genome_version,tool,metric),
 FOREIGN KEY(genome_version) REFERENCES genome_info(genome_version));
"""

EVENT_SQL = """
INSERT INTO splicing_events(
 junctionID,eventID,geneID,seqID,strand,event_type,description,start,end,
 coord,full_coord,upstream_exon,downstream_exon,species,genome_version)
 VALUES (?,?,'G','chr1','+','SE','I',10,20,'10-20','chr1:10-20',NULL,NULL,'Ath','ATH')
"""


def build(path, sample=None, psi=None, extra_sample=None):
    with sqlite3.connect(path) as conn:
        conn.executescript(SCHEMA)
        conn.execute("INSERT INTO genome_info VALUES ('ATH','Arabidopsis')")
        if sample:
            conn.execute(EVENT_SQL, ("ATH_J1", "E1"))
            conn.execute(EVENT_SQL, ("ATH_J2", "E2"))
            conn.execute(
                "INSERT INTO junction_sample_stats "
                "(junctionID,psi_majiq,psi_sgseq,reads_majiq,"
                "reads_sgseq_3p,reads_sgseq_5p,sra_id,species) "
                "VALUES ('ATH_J1',?, ?,10,2,3,?,'Ath')",
                (psi, psi, sample),
            )
            conn.execute(
                "INSERT INTO junction_sample_stats "
                "(junctionID,psi_majiq,psi_sgseq,reads_majiq,"
                "reads_sgseq_3p,reads_sgseq_5p,sra_id,species) "
                "VALUES ('ATH_J2',NULL,0.4,NULL,5,6,?,'Ath')",
                (sample,),
            )
            for tool in ("STAR", "BBDUK"):
                conn.execute(
                    "INSERT INTO sample_qc_metrics "
                    "(sra_id,genome_version,tool,metric,numeric_value) "
                    "VALUES (?,'ATH',?,'reads',100)",
                    (sample, tool),
                )
            if extra_sample:
                conn.execute(
                    "INSERT INTO sample_qc_metrics "
                    "(sra_id,genome_version,tool,metric,numeric_value) "
                    "VALUES (?,'ATH','STAR','reads',100)", (extra_sample,),
                )


def snapshot(conn):
    tables = ("genome_info", "splicing_events", "junction_sample_stats",
              "junction_summary_stats", "sample_qc_metrics", "sample_processing_status")
    return {
        name: conn.execute(f"SELECT * FROM {name} ORDER BY 1").fetchall()
        for name in tables
    }


class MicrobatchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dest = self.root / "dest.sqlite"
        self.stage_a = self.root / "A.sqlite"
        self.stage_b = self.root / "B.sqlite"
        build(self.dest)
        build(self.stage_a, "A", .2)
        build(self.stage_b, "B", .6)

    def test_two_microbatches_match_reference_sql_aggregation(self):
        self.assertEqual(commit_microbatch(self.dest, [("A", self.stage_a)], "ATH"), 1)
        self.assertEqual(commit_microbatch(self.dest, [("B", self.stage_b)], "ATH"), 1)
        with sqlite3.connect(self.dest) as c:
            self.assertEqual(c.execute("PRAGMA integrity_check").fetchone()[0], "ok")
            self.assertEqual(c.execute("PRAGMA foreign_key_check").fetchall(), [])
            self.assertEqual(c.execute("SELECT COUNT(*) FROM sample_processing_status").fetchone()[0], 2)
            actual = c.execute("SELECT mean_psi_majiq,mean_psi_sgseq,"
                               "sum_reads_majiq,sum_reads_sgseq_3p,sum_reads_sgseq_5p,"
                               "n_of_samples_majiq,n_of_samples_sgseq "
                               "FROM junction_summary_stats WHERE junctionID='ATH_J1'").fetchone()
            expected = c.execute("SELECT AVG(psi_majiq),AVG(psi_sgseq),"
                                 "SUM(reads_majiq),SUM(reads_sgseq_3p),SUM(reads_sgseq_5p),"
                                 "COUNT(psi_majiq),COUNT(psi_sgseq) "
                                 "FROM junction_sample_stats WHERE junctionID='ATH_J1'").fetchone()
            self.assertEqual(actual, expected)
            self.assertAlmostEqual(actual[0], .4)
            self.assertEqual(c.execute("SELECT COUNT(*) FROM junction_sample_stats").fetchone()[0], 4)

    def test_failure_rolls_back_both_samples_and_registry(self):
        with sqlite3.connect(self.dest) as c:
            before = snapshot(c)
        with self.assertRaisesRegex(RuntimeError, "TEST_INJECTED"):
            commit_microbatch(self.dest,
                              [("A", self.stage_a), ("B", self.stage_b)],
                              "ATH", fail_before_commit=True)
        with sqlite3.connect(self.dest) as c:
            self.assertEqual(snapshot(c), before)

    def test_duplicate_rejected_with_no_change(self):
        commit_microbatch(self.dest, [("A", self.stage_a)], "ATH")
        with sqlite3.connect(self.dest) as c:
            before = snapshot(c)
        with self.assertRaisesRegex(ValueError, "already registered"):
            commit_microbatch(self.dest, [("A", self.stage_a)], "ATH")
        with sqlite3.connect(self.dest) as c:
            self.assertEqual(snapshot(c), before)

    def test_contaminated_stage_rejected_with_no_change(self):
        bad = self.root / "bad.sqlite"
        build(bad, "A", .2, extra_sample="OTHER")
        with sqlite3.connect(self.dest) as c:
            before = snapshot(c)
        with self.assertRaisesRegex(ValueError, "Expected stage"):
            commit_microbatch(self.dest, [("A", bad)], "ATH")
        with sqlite3.connect(self.dest) as c:
            self.assertEqual(snapshot(c), before)

    def test_second_sample_failure_rolls_back_first(self):
        bad = self.root / "bad.sqlite"
        build(bad, "B", .6, extra_sample="OTHER")
        with sqlite3.connect(self.dest) as c:
            before = snapshot(c)
        with self.assertRaises(ValueError):
            commit_microbatch(self.dest, [("A", self.stage_a), ("B", bad)], "ATH")
        with sqlite3.connect(self.dest) as c:
            self.assertEqual(snapshot(c), before)


if __name__ == "__main__":
    unittest.main()
