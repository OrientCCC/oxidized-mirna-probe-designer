"""Normalization and sequence utility helpers."""

from __future__ import annotations

from probe_designer.models import ProbeDesignError

_RNA_BASES = frozenset("AUCG")
_RNA_TO_DNA_COMPLEMENT = {
    "A": "T",
    "U": "A",
    "C": "G",
    "G": "C",
}


def normalize_mirna_sequence(seq: str) -> tuple[str, list[str]]:
    """Normalize user input into an uppercase RNA sequence.

    Whitespace is stripped. Any ``T`` characters are converted to ``U`` and
    reported back as a warning so the normalization is transparent.
    """

    collapsed = "".join(seq.split()).upper()
    if not collapsed:
        raise ProbeDesignError("miRNA sequence is empty after removing whitespace.")

    warnings: list[str] = []
    if "T" in collapsed:
        collapsed = collapsed.replace("T", "U")
        warnings.append(
            "Input contained T bases; converted T -> U during RNA normalization."
        )

    invalid_chars = sorted({char for char in collapsed if char not in _RNA_BASES})
    if invalid_chars:
        invalid_display = ", ".join(invalid_chars)
        raise ProbeDesignError(
            "miRNA sequence must contain only RNA bases A/U/C/G after normalization. "
            f"Invalid characters: {invalid_display}."
        )

    return collapsed, warnings


def rna_to_dna_reverse_complement(rna_seq: str) -> str:
    """Convert an RNA target window into a DNA reverse-complement oligo."""

    if not rna_seq:
        raise ProbeDesignError("Target RNA window cannot be empty.")

    invalid_chars = sorted({char for char in rna_seq if char not in _RNA_BASES})
    if invalid_chars:
        invalid_display = ", ".join(invalid_chars)
        raise ProbeDesignError(
            "RNA target window must contain only A/U/C/G. "
            f"Invalid characters: {invalid_display}."
        )

    return "".join(_RNA_TO_DNA_COMPLEMENT[base] for base in reversed(rna_seq))


def slice_one_based(seq: str, start: int, end: int) -> str:
    """Slice a sequence using inclusive 1-based coordinates."""

    return seq[start - 1 : end]


def format_range(start: int, end: int) -> str:
    """Render inclusive coordinates in a compact human-readable form."""

    if start == end:
        return str(start)
    return f"{start}-{end}"
