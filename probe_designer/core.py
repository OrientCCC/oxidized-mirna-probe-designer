"""Core rule-based miRNA probe design logic."""

from __future__ import annotations

from probe_designer.constants import (
    DEFAULT_ADD_PROBE_L_PHOSPHATE_TAG,
    DEFAULT_FORMAT_PROBE_R_LAST_TWO_AS_RNA,
    DEFAULT_SPLIT_MODE,
    LIGATED_PRODUCT_ORIENTATION,
    PROBE_L_ADAPTOR,
    PROBE_R_ADAPTOR,
    SPLIT_MODE_DESCRIPTIONS,
    SPLIT_MODE_OXIDIZED_TO_R,
    SUPPORTED_SPLIT_MODES,
    UNIVERSAL_PRIMER_F,
    UNIVERSAL_PRIMER_R,
)
from probe_designer.formatting import build_probe_r_display
from probe_designer.models import DesignResult, ProbeDesignError
from probe_designer.utils import format_range, normalize_mirna_sequence, rna_to_dna_reverse_complement
from probe_designer.utils import slice_one_based


def _validate_split_mode(split_mode: str) -> None:
    if split_mode not in SUPPORTED_SPLIT_MODES:
        supported = ", ".join(SUPPORTED_SPLIT_MODES)
        raise ProbeDesignError(
            f"Unsupported split_mode '{split_mode}'. Supported values: {supported}."
        )


def _calculate_windows(
    seq_len: int,
    oxidized_pos: int,
    split_mode: str,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Resolve inclusive 1-based full-span miRNA windows for Probe L and Probe R.

    The current assay rules require both probes to retain at least one
    miRNA-complementary base. Splits that would leave either side empty are
    rejected explicitly.
    """

    if not 1 <= oxidized_pos <= seq_len:
        raise ProbeDesignError(
            "oxidized_pos is out of range for the cleaned miRNA sequence. "
            f"Received {oxidized_pos}; valid range is 1-{seq_len}."
        )

    if split_mode == DEFAULT_SPLIT_MODE:
        probe_l_start = 1
        probe_l_end = oxidized_pos - 1
        probe_r_start = oxidized_pos
        probe_r_end = seq_len
    else:
        probe_l_start = 1
        probe_l_end = oxidized_pos
        probe_r_start = oxidized_pos + 1
        probe_r_end = seq_len

    if probe_l_start > probe_l_end:
        raise ProbeDesignError(
            "The selected split leaves Probe L with zero miRNA-complementary bases. "
            f"With split_mode='{split_mode}' and oxidized_pos={oxidized_pos}, "
            "Probe L would be empty."
        )

    if probe_r_start > probe_r_end:
        raise ProbeDesignError(
            "The selected split leaves Probe R with zero miRNA-complementary bases. "
            f"With split_mode='{split_mode}' and oxidized_pos={oxidized_pos}, "
            "Probe R would be empty."
        )

    return (probe_l_start, probe_l_end), (probe_r_start, probe_r_end)


def _build_schematic(
    cleaned_mirna_seq_rna: str,
    split_mode: str,
    probe_l_start: int,
    probe_l_end: int,
    probe_r_start: int,
    probe_r_end: int,
) -> str:
    return "\n".join(
        [
            f"miRNA (5'->3'): {cleaned_mirna_seq_rna}",
            f"Split rule: {SPLIT_MODE_DESCRIPTIONS[split_mode]}",
            f"Probe L covers: positions {probe_l_start}..{probe_l_end}",
            f"Probe R covers: positions {probe_r_start}..{probe_r_end}",
            f"Ligated product: {LIGATED_PRODUCT_ORIENTATION}",
        ]
    )


def design_probes(
    mirna_seq: str,
    oxidized_pos: int,
    split_mode: str = DEFAULT_SPLIT_MODE,
    add_probe_l_phosphate_tag: bool = DEFAULT_ADD_PROBE_L_PHOSPHATE_TAG,
    format_probe_r_last_two_as_rna: bool = DEFAULT_FORMAT_PROBE_R_LAST_TWO_AS_RNA,
) -> DesignResult:
    """Design rule-based Probe L and Probe R sequences from an RNA miRNA target.

    The input miRNA is treated as RNA written 5'->3'. Probe target-binding regions
    are returned as DNA oligos written 5'->3', meaning they are the DNA
    reverse-complements of the assigned miRNA windows.
    """

    _validate_split_mode(split_mode)

    cleaned_mirna_seq_rna, normalization_warnings = normalize_mirna_sequence(mirna_seq)
    (probe_l_start, probe_l_end), (probe_r_start, probe_r_end) = _calculate_windows(
        seq_len=len(cleaned_mirna_seq_rna),
        oxidized_pos=oxidized_pos,
        split_mode=split_mode,
    )

    probe_l_target_rna = slice_one_based(cleaned_mirna_seq_rna, probe_l_start, probe_l_end)
    probe_r_target_rna = slice_one_based(cleaned_mirna_seq_rna, probe_r_start, probe_r_end)

    probe_l_binding_dna = rna_to_dna_reverse_complement(probe_l_target_rna)
    probe_r_binding_dna = rna_to_dna_reverse_complement(probe_r_target_rna)

    probe_l_prefix = "/PO4/" if add_probe_l_phosphate_tag else ""
    probe_l_full = f"{probe_l_prefix}{probe_l_binding_dna}{PROBE_L_ADAPTOR}"
    probe_r_full_raw = f"{PROBE_R_ADAPTOR}{probe_r_binding_dna}"
    probe_r_full_display = build_probe_r_display(
        probe_r_full_raw,
        enabled=format_probe_r_last_two_as_rna,
    )

    notes = [
        "Input miRNA is interpreted as RNA written 5'->3'.",
        "Probe L automatically covers the full left miRNA segment defined by the split.",
        "Probe R automatically covers the full right miRNA segment defined by the split.",
        "Probe target-binding regions are DNA reverse complements of the assigned RNA windows.",
        "The two probes together span the full miRNA with no gap; split mode only decides which probe contains the oxidized base.",
        f"Ligated probe product orientation is {LIGATED_PRODUCT_ORIENTATION}.",
        (
            "Version 1 default split mode is "
            f"'{SPLIT_MODE_OXIDIZED_TO_R}', which assigns the oxidized base to Probe R."
        ),
        f"Selected split mode: {SPLIT_MODE_DESCRIPTIONS[split_mode]}.",
        "Designs are rejected if either probe would have zero miRNA-complementary bases.",
        (
            "Probe R display formatting "
            + (
                "annotates the final two displayed bases as rN."
                if format_probe_r_last_two_as_rna
                else "is disabled, so displayed Probe R matches the raw sequence."
            )
        ),
        (
            "Probe L ordering output "
            + (
                "includes the literal /PO4/ tag."
                if add_probe_l_phosphate_tag
                else "does not include the /PO4/ tag."
            )
        ),
    ]

    return DesignResult(
        cleaned_mirna_seq_rna=cleaned_mirna_seq_rna,
        oxidized_pos=oxidized_pos,
        split_mode=split_mode,
        probe_l_target_start=probe_l_start,
        probe_l_target_end=probe_l_end,
        probe_r_target_start=probe_r_start,
        probe_r_target_end=probe_r_end,
        probe_l_target_rna=probe_l_target_rna,
        probe_r_target_rna=probe_r_target_rna,
        probe_l_binding_dna=probe_l_binding_dna,
        probe_r_binding_dna=probe_r_binding_dna,
        probe_l_full=probe_l_full,
        probe_r_full_raw=probe_r_full_raw,
        probe_r_full_display=probe_r_full_display,
        universal_primer_f=UNIVERSAL_PRIMER_F,
        universal_primer_r=UNIVERSAL_PRIMER_R,
        notes=notes,
        warnings=normalization_warnings,
        schematic=_build_schematic(
            cleaned_mirna_seq_rna=cleaned_mirna_seq_rna,
            split_mode=split_mode,
            probe_l_start=probe_l_start,
            probe_l_end=probe_l_end,
            probe_r_start=probe_r_start,
            probe_r_end=probe_r_end,
        ),
    )


def build_design_report(result: DesignResult) -> str:
    """Render a copy-friendly text report for UI display or download."""

    lines = [
        "Rule-Based miRNA Probe Design Report",
        "",
        result.schematic,
        "",
        "Inputs",
        f"- Cleaned miRNA sequence (RNA, 5'->3'): {result.cleaned_mirna_seq_rna}",
        f"- Oxidized position (1-based): {result.oxidized_pos}",
        f"- Split mode: {result.split_mode}",
        "",
        "Probe windows",
        (
            "- Probe L assigned positions: "
            f"{format_range(result.probe_l_target_start, result.probe_l_target_end)} "
            f"| RNA subsequence: {result.probe_l_target_rna} "
            f"| DNA binding: {result.probe_l_binding_dna}"
        ),
        (
            "- Probe R assigned positions: "
            f"{format_range(result.probe_r_target_start, result.probe_r_target_end)} "
            f"| RNA subsequence: {result.probe_r_target_rna} "
            f"| DNA binding: {result.probe_r_binding_dna}"
        ),
        "",
        "Full oligos",
        f"- Probe R full raw: {result.probe_r_full_raw}",
        f"- Probe R full display: {result.probe_r_full_display}",
        f"- Probe L full: {result.probe_l_full}",
        "",
        "Universal primers",
        f"- Universal primer F: {result.universal_primer_f}",
        f"- Universal primer R: {result.universal_primer_r}",
    ]

    if result.warnings:
        lines.extend(["", "Warnings"])
        lines.extend(f"- {warning}" for warning in result.warnings)

    if result.notes:
        lines.extend(["", "Notes"])
        lines.extend(f"- {note}" for note in result.notes)

    return "\n".join(lines)
