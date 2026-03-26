"""Typed models used by the rule-based probe designer."""

from __future__ import annotations

from dataclasses import dataclass, field


class ProbeDesignError(ValueError):
    """Raised when the requested probe design is invalid."""


@dataclass(slots=True)
class DesignResult:
    """Structured output returned by :func:`probe_designer.core.design_probes`."""

    cleaned_mirna_seq_rna: str
    oxidized_pos: int
    split_mode: str
    probe_l_target_start: int
    probe_l_target_end: int
    probe_r_target_start: int
    probe_r_target_end: int
    probe_l_target_rna: str
    probe_r_target_rna: str
    probe_l_binding_dna: str
    probe_r_binding_dna: str
    probe_l_full: str
    probe_r_full_raw: str
    probe_r_full_display: str
    universal_primer_f: str
    universal_primer_r: str
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    schematic: str = ""

    def as_flat_dict(self) -> dict[str, str | int]:
        """Return a flat mapping suitable for CSV export."""

        return {
            "cleaned_mirna_seq_rna": self.cleaned_mirna_seq_rna,
            "oxidized_pos": self.oxidized_pos,
            "split_mode": self.split_mode,
            "probe_l_target_start": self.probe_l_target_start,
            "probe_l_target_end": self.probe_l_target_end,
            "probe_r_target_start": self.probe_r_target_start,
            "probe_r_target_end": self.probe_r_target_end,
            "probe_l_target_rna": self.probe_l_target_rna,
            "probe_r_target_rna": self.probe_r_target_rna,
            "probe_l_binding_dna": self.probe_l_binding_dna,
            "probe_r_binding_dna": self.probe_r_binding_dna,
            "probe_l_full": self.probe_l_full,
            "probe_r_full_raw": self.probe_r_full_raw,
            "probe_r_full_display": self.probe_r_full_display,
            "universal_primer_f": self.universal_primer_f,
            "universal_primer_r": self.universal_primer_r,
            "notes": " | ".join(self.notes),
            "warnings": " | ".join(self.warnings),
            "schematic": self.schematic,
        }
