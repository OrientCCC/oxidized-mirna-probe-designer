"""Display-only formatting helpers for designed probes."""

from __future__ import annotations

from probe_designer.models import ProbeDesignError

_DNA_TO_RNA_DISPLAY = {
    "A": "rA",
    "C": "rC",
    "G": "rG",
    "T": "rT",
}


def format_probe_r_last_two_as_rna_display(seq: str) -> str:
    """Annotate the final two DNA bases of Probe R in ``rN`` display format."""

    if len(seq) < 2:
        raise ProbeDesignError(
            "Probe R full sequence must be at least 2 bases long to annotate the "
            "final two nucleotides."
        )

    tail = seq[-2:]
    invalid_chars = sorted({char for char in tail if char not in _DNA_TO_RNA_DISPLAY})
    if invalid_chars:
        invalid_display = ", ".join(invalid_chars)
        raise ProbeDesignError(
            "Probe R display formatting expects canonical DNA bases in the final "
            f"two positions. Invalid characters: {invalid_display}."
        )

    formatted_tail = "".join(_DNA_TO_RNA_DISPLAY[base] for base in tail)
    return f"{seq[:-2]}{formatted_tail}"


def build_probe_r_display(seq: str, enabled: bool = True) -> str:
    """Return either the raw Probe R sequence or the display-formatted version."""

    if not enabled:
        return seq
    return format_probe_r_last_two_as_rna_display(seq)
