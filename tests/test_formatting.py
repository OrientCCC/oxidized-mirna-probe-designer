"""Unit tests for Probe R display formatting helpers."""

from __future__ import annotations

import unittest

from probe_designer.formatting import build_probe_r_display, format_probe_r_last_two_as_rna_display
from probe_designer.models import ProbeDesignError


class ProbeDesignerFormattingTests(unittest.TestCase):
    def test_probe_r_display_formats_terminal_gc_as_rgrc(self) -> None:
        self.assertEqual(format_probe_r_last_two_as_rna_display("AAAGC"), "AAArGrC")

    def test_probe_r_display_formats_terminal_ac_as_rarc(self) -> None:
        self.assertEqual(format_probe_r_last_two_as_rna_display("TTAC"), "TTrArC")

    def test_probe_r_display_can_be_disabled(self) -> None:
        self.assertEqual(build_probe_r_display("TTAC", enabled=False), "TTAC")

    def test_probe_r_display_requires_at_least_two_bases(self) -> None:
        with self.assertRaisesRegex(ProbeDesignError, "at least 2 bases"):
            format_probe_r_last_two_as_rna_display("A")


if __name__ == "__main__":
    unittest.main()
