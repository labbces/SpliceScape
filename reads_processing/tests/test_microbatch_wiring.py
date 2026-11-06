"""Static safety invariants for opt-in incremental original-parser workflow."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]

class MicrobatchWiringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workflow = (ROOT / "splicescape.nf").read_text()
        cls.modules = (ROOT / "modules.nf").read_text()
        cls.config = (ROOT / "splicescape-v2.config").read_text()

    def test_opt_in_and_existing_end_run_import_preserved(self):
        self.assertIn("if (params.db_batch_size != null && (params.db_batch_size as int) > 0)", self.workflow)
        self.assertIn("paired_batch = paired_results", self.workflow)
        self.assertIn("POPULATE_DB_BATCH (", self.workflow)
        self.assertIn("FINALIZE_SAMPLE(", self.workflow)

    def test_terminal_batch_includes_failure_and_remainder(self):
        self.assertIn("micro_success.mix(micro_failures)", self.workflow)
        self.assertIn(".collate(dbBatchSize)", self.workflow)
        self.assertIn("early_failures", self.workflow)
        self.assertIn("branch_terminal", self.workflow)

    def test_original_parser_and_finalizer_in_same_process(self):
        start = self.modules.index("process COMMIT_DB_MICROBATCH {")
        proc = self.modules[start:]
        self.assertIn("cache false", proc)
        self.assertIn("maxForks 1", proc)
        self.assertIn("bin/import_batch.py", proc)
        self.assertIn("bin/finalize_samples.py", proc)
        self.assertLess(proc.index("bin/import_batch.py"), proc.index("bin/finalize_samples.py"))
        self.assertIn("--existing-db", proc)

    def test_config_single_writer_is_explicit(self):
        self.assertEqual(self.config.count("withName: 'COMMIT_DB_MICROBATCH'"), 2)
        self.assertEqual(self.config.count("withName: 'COMMIT_DB_MICROBATCH' {\n        cpus = 8\n        memory = '24 GB'\n        maxForks = 1"), 2)

if __name__ == "__main__":
    unittest.main()
