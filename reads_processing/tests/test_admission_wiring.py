#!/usr/bin/env python3
"""Static guard against stale admission semaphore calls in the actual pipeline.

Run: python3 reads_processing/tests/test_admission_wiring.py
"""
from pathlib import Path
import re
import unittest

PIPELINE = Path(__file__).resolve().parents[1] / "splicescape.nf"


class AdmissionWiringTests(unittest.TestCase):
    def test_terminal_callbacks_use_controller(self):
        script = PIPELINE.read_text(encoding="utf-8")
        self.assertNotIn("admissionSlots", script,
                         "Legacy admission semaphore reference remains")
        self.assertNotIn("Semaphore(", script,
                         "Blocking semaphore was reintroduced")
        self.assertEqual(script.count("admissionController.terminal(sample)"), 3,
                         "Expected one early-failure and two analysis terminal callbacks")
        self.assertIn("admissionController.start()", script)

    def test_terminal_callbacks_are_attached_to_correct_outcomes(self):
        script = PIPELINE.read_text(encoding="utf-8")
        early = re.search(r"early_failures\.subscribe\s*\{.*?\n\s*\}", script, re.S)
        branch = re.search(r"branch_terminal\.subscribe\s*\{.*?\n\s*\}\s*\n\s*\}", script, re.S)
        self.assertIsNotNone(early)
        self.assertIsNotNone(branch)
        self.assertEqual(early.group().count("admissionController.terminal(sample)"), 1)
        self.assertEqual(branch.group().count("admissionController.terminal(sample)"), 2)


if __name__ == "__main__":
    unittest.main()
