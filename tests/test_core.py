"""Unit tests for core probe-design behavior."""

from __future__ import annotations

import unittest

from probe_designer.constants import PROBE_L_ADAPTOR, PROBE_R_ADAPTOR
from probe_designer.core import design_probes
from probe_designer.models import ProbeDesignError
from probe_designer.utils import normalize_mirna_sequence, rna_to_dna_reverse_complement


class ProbeDesignerCoreTests(unittest.TestCase):
    def test_normalize_mirna_sequence_strips_whitespace_and_converts_t_to_u(self) -> None:
        cleaned, warnings = normalize_mirna_sequence(" au tg \nca ")
        self.assertEqual(cleaned, "AUUGCA")
        self.assertEqual(
            warnings,
            ["Input contained T bases; converted T -> U during RNA normalization."],
        )

    def test_rna_to_dna_reverse_complement_matches_orientation_rule(self) -> None:
        self.assertEqual(rna_to_dna_reverse_complement("AUGC"), "GCAT")

    def test_full_span_split_mode_assigns_oxidized_base_to_probe_r_by_default(self) -> None:
        result = design_probes(
            mirna_seq="AUGCUAGC",
            oxidized_pos=4,
        )
        self.assertEqual(result.split_mode, "oxidized_base_belongs_to_R")
        self.assertEqual((result.probe_l_target_start, result.probe_l_target_end), (1, 3))
        self.assertEqual((result.probe_r_target_start, result.probe_r_target_end), (4, 8))
        self.assertEqual(result.probe_l_target_rna, "AUG")
        self.assertEqual(result.probe_r_target_rna, "CUAGC")
        self.assertEqual(result.probe_l_binding_dna, "CAT")
        self.assertEqual(result.probe_r_binding_dna, "GCTAG")

    def test_full_span_split_mode_can_assign_oxidized_base_to_probe_l(self) -> None:
        result = design_probes(
            mirna_seq="AUGCUAGC",
            oxidized_pos=4,
            split_mode="oxidized_base_belongs_to_L",
            format_probe_r_last_two_as_rna=False,
        )
        self.assertEqual((result.probe_l_target_start, result.probe_l_target_end), (1, 4))
        self.assertEqual((result.probe_r_target_start, result.probe_r_target_end), (5, 8))
        self.assertEqual(result.probe_l_target_rna, "AUGC")
        self.assertEqual(result.probe_r_target_rna, "UAGC")
        self.assertEqual(result.probe_l_binding_dna, "GCAT")
        self.assertEqual(result.probe_r_binding_dna, "GCTA")

    def test_full_span_splitting_has_no_gap_and_no_overlap(self) -> None:
        result = design_probes(
            mirna_seq="AUGCUAGC",
            oxidized_pos=4,
            split_mode="oxidized_base_belongs_to_R",
            format_probe_r_last_two_as_rna=False,
        )
        probe_l_positions = set(range(result.probe_l_target_start, result.probe_l_target_end + 1))
        probe_r_positions = set(range(result.probe_r_target_start, result.probe_r_target_end + 1))
        self.assertEqual(probe_l_positions | probe_r_positions, set(range(1, 9)))
        self.assertEqual(probe_l_positions & probe_r_positions, set())

    def test_orientation_consistency_and_adaptor_assembly(self) -> None:
        result = design_probes(
            mirna_seq="AUGCAU",
            oxidized_pos=3,
            add_probe_l_phosphate_tag=False,
            format_probe_r_last_two_as_rna=False,
        )
        self.assertEqual(result.probe_l_target_rna, "AU")
        self.assertEqual(result.probe_r_target_rna, "GCAU")
        self.assertEqual(result.probe_l_binding_dna, "AT")
        self.assertEqual(result.probe_r_binding_dna, "ATGC")
        self.assertTrue(result.probe_r_full_raw.startswith(PROBE_R_ADAPTOR))
        self.assertEqual(result.probe_r_full_raw, f"{PROBE_R_ADAPTOR}ATGC")
        self.assertTrue(result.probe_l_full.endswith(PROBE_L_ADAPTOR))
        self.assertEqual(result.probe_l_full, f"AT{PROBE_L_ADAPTOR}")

    def test_invalid_inputs_raise_clear_errors(self) -> None:
        invalid_cases = [
            (
                {
                    "mirna_seq": "",
                    "oxidized_pos": 1,
                },
                "empty",
            ),
            (
                {
                    "mirna_seq": "AUGX",
                    "oxidized_pos": 2,
                },
                "Invalid characters",
            ),
            (
                {
                    "mirna_seq": "AUGC",
                    "oxidized_pos": 5,
                },
                "out of range",
            ),
            (
                {
                    "mirna_seq": "AUGC",
                    "oxidized_pos": 1,
                    "split_mode": "oxidized_base_belongs_to_R",
                },
                "Probe L would be empty",
            ),
            (
                {
                    "mirna_seq": "AUGC",
                    "oxidized_pos": 4,
                    "split_mode": "oxidized_base_belongs_to_L",
                },
                "Probe R would be empty",
            ),
        ]

        for kwargs, message_fragment in invalid_cases:
            with self.subTest(kwargs=kwargs):
                with self.assertRaisesRegex(ProbeDesignError, message_fragment):
                    design_probes(**kwargs)


if __name__ == "__main__":
    unittest.main()
