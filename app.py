"""Streamlit web app for the rule-based miRNA probe designer."""

from __future__ import annotations

import csv
import io
from pathlib import Path

import streamlit as st

from probe_designer import (
    DEFAULT_ADD_PROBE_L_PHOSPHATE_TAG,
    DEFAULT_FORMAT_PROBE_R_LAST_TWO_AS_RNA,
    DEFAULT_SPLIT_MODE,
    LIGATED_PRODUCT_ORIENTATION,
    PROBE_L_ADAPTOR,
    SPLIT_MODE_DESCRIPTIONS,
    SUPPORTED_SPLIT_MODES,
    ProbeDesignError,
    build_design_report,
    design_probes,
)
from probe_designer.utils import format_range

SCHEMATIC_IMAGE_PATH = Path(__file__).resolve().parent / "assets" / "Probe_Design_Scheme.png"


def _inject_styles() -> None:
    st.markdown(
        """
        <style>
        .block-container {
            max-width: 1120px;
            padding-top: 2.0rem;
            padding-bottom: 2.0rem;
        }
        .hero {
            padding: 1.25rem 1.5rem;
            border-radius: 18px;
            background:
                radial-gradient(circle at top right, rgba(0, 102, 128, 0.18), transparent 34%),
                linear-gradient(135deg, #eef7f6 0%, #f6f3eb 100%);
            border: 1px solid rgba(0, 102, 128, 0.14);
            margin-bottom: 1rem;
        }
        .hero h1 {
            margin: 0 0 0.4rem 0;
            color: #12343b;
            letter-spacing: -0.02em;
        }
        .hero p {
            margin: 0;
            color: #35525a;
            line-height: 1.45;
        }
        .section-note {
            padding: 0.9rem 1rem;
            border-left: 4px solid #006680;
            background: #f5fafb;
            border-radius: 0.4rem;
            margin: 0.5rem 0 1rem 0;
        }
        .primary-note {
            min-height: 5.6rem;
            color: #8a929d;
            font-size: 0.95rem;
            line-height: 1.65;
            margin: 0.35rem 0 0.85rem 0;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def _render_header() -> None:
    st.markdown(
        """
        <div class="hero">
            <h1>Site-Specific Oxidized miRNA Probe Designer</h1>
            <p>
                Input a miRNA in RNA 5'->3' orientation, choose the oxidized guanine
                position, and generate Probe L / Probe R binding regions plus assembled
                oligos with transparent, configurable split rules.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.markdown(
        """
        <div class="section-note">
            Probe L automatically covers the full left miRNA segment, Probe R
            automatically covers the full right segment, binding regions are DNA
            reverse complements, and version 1 defaults to assigning the oxidized
            base to Probe R.
        </div>
        """,
        unsafe_allow_html=True,
    )


def _format_split_mode(mode: str) -> str:
    return SPLIT_MODE_DESCRIPTIONS[mode]


def _window_summary(label: str, start: int, end: int, rna_seq: str, binding_dna: str) -> None:
    st.markdown(f"**{label} assigned miRNA positions**")
    st.code(format_range(start, end), language="text")
    st.markdown(f"**{label} assigned RNA subsequence (5'->3')**")
    st.code(rna_seq, language="text")
    st.markdown(f"**{label} DNA binding region (5'->3')**")
    st.code(binding_dna, language="text")


def _build_probe_l_variants(binding_dna: str) -> tuple[str, str]:
    probe_l_without_phosphate = f"{binding_dna}{PROBE_L_ADAPTOR}"
    probe_l_with_phosphate = f"/PO4/{probe_l_without_phosphate}"
    return probe_l_without_phosphate, probe_l_with_phosphate


def _render_primary_note(lines: list[str]) -> None:
    joined = "<br>".join(lines)
    st.markdown(f'<div class="primary-note">{joined}</div>', unsafe_allow_html=True)


def _result_to_csv_bytes(result) -> bytes:
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["field", "value"])
    for key, value in result.as_flat_dict().items():
        writer.writerow([key, value])
    return buffer.getvalue().encode("utf-8")


def _result_to_txt_bytes(result) -> bytes:
    return build_design_report(result).encode("utf-8")


def main() -> None:
    st.set_page_config(
        page_title="Oxidized miRNA Probe Designer",
        layout="wide",
    )
    _inject_styles()
    _render_header()

    with st.sidebar:
        st.header("Design Controls")
        split_mode = st.radio(
            "Split mode",
            options=list(SUPPORTED_SPLIT_MODES),
            format_func=_format_split_mode,
            index=list(SUPPORTED_SPLIT_MODES).index(DEFAULT_SPLIT_MODE),
            help="Version 1 defaults to assigning the oxidized base to Probe R.",
        )
        add_probe_l_phosphate_tag = st.checkbox(
            "Add /PO4/ prefix to Probe L output",
            value=DEFAULT_ADD_PROBE_L_PHOSPHATE_TAG,
        )
        format_probe_r_last_two_as_rna = st.checkbox(
            "Annotate final two Probe R bases as rN",
            value=DEFAULT_FORMAT_PROBE_R_LAST_TWO_AS_RNA,
            help="Display-only annotation of the final two Probe R bases, for example ...GC -> ...rGrC.",
        )

    with st.form("probe_designer_form"):
        mirna_seq = st.text_area(
            "miRNA sequence (RNA, 5'->3')",
            placeholder="Example: UGAGGUAGUAGGUUGUAUAGUU",
            height=120,
        )
        oxidized_pos = st.number_input(
            "Oxidized guanine position (1-based index)",
            min_value=1,
            value=1,
            step=1,
        )
        submitted = st.form_submit_button("Generate Probes", use_container_width=True)

    if not submitted:
        with st.expander("Design assumptions", expanded=True):
            st.markdown(
                "\n".join(
                    [
                        "- Input miRNA is RNA written 5'->3'.",
                        "- Probe L automatically covers the full left miRNA segment defined by the split.",
                        "- Probe R automatically covers the full right miRNA segment defined by the split.",
                        "- Probe L and Probe R binding regions are DNA reverse complements of their assigned RNA subsequences.",
                        f"- Ligated product orientation is {LIGATED_PRODUCT_ORIENTATION}.",
                        "- Version 1 defaults to split mode `oxidized_base_belongs_to_R`.",
                        "- The app rejects split positions that would leave Probe L or Probe R with zero miRNA-complementary bases.",
                        "- Probe L can optionally include a literal `/PO4/` ordering tag.",
                        "- Probe R can optionally annotate its final two displayed bases in `rN` notation without changing the underlying raw sequence.",
                        "- Version 1 does not model Tm optimization, mismatch scoring, ligation efficiency, or secondary structure.",
                    ]
                )
            )
        return

    try:
        result = design_probes(
            mirna_seq=mirna_seq,
            oxidized_pos=int(oxidized_pos),
            split_mode=split_mode,
            add_probe_l_phosphate_tag=add_probe_l_phosphate_tag,
            format_probe_r_last_two_as_rna=format_probe_r_last_two_as_rna,
        )
    except ProbeDesignError as exc:
        st.error(str(exc))
        return

    st.success("Probe design generated.")

    top_left, top_mid, top_right = st.columns([2.2, 1, 1])
    with top_left:
        st.markdown("**Cleaned miRNA sequence (RNA, 5'->3')**")
        st.code(result.cleaned_mirna_seq_rna, language="text")
    with top_mid:
        st.metric("Oxidized position", result.oxidized_pos)
    with top_right:
        st.markdown("**Split mode**")
        st.code(SPLIT_MODE_DESCRIPTIONS[result.split_mode], language="text")

    st.subheader("Schematic")
    st.code(result.schematic, language="text")
    if SCHEMATIC_IMAGE_PATH.exists():
        with st.expander("Assay schematic image", expanded=False):
            st.info(
                "This figure is a conceptual schematic for understanding probe "
                "orientation, adaptors, /PO4/ notation, and Probe R rN display. "
                "It is not a literal or fully realistic experimental diagram."
            )
            st.image(
                str(SCHEMATIC_IMAGE_PATH),
                caption=(
                    "Conceptual schematic only. Used to help interpret Probe L / "
                    "Probe R orientation, adaptors, /PO4/ notation, and Probe R "
                    "rN display; not intended as a literal real-world drawing."
                ),
                use_container_width=True,
            )

    probe_l_col, probe_r_col = st.columns(2)
    with probe_l_col:
        st.subheader("Probe L")
        probe_l_without_phosphate, probe_l_with_phosphate = _build_probe_l_variants(
            result.probe_l_binding_dna
        )
        selected_probe_l_variant = (
            "Current ordering option: with /PO4/"
            if add_probe_l_phosphate_tag
            else "Current ordering option: without /PO4/"
        )
        st.markdown("**Recommended full sequence (5'->3')**")
        st.code(probe_l_with_phosphate, language="text")
        _render_primary_note(
            [
                "Recommended ordering display with /PO4/.",
                selected_probe_l_variant,
            ]
        )
        st.markdown("---")
        _window_summary(
            label="Probe L",
            start=result.probe_l_target_start,
            end=result.probe_l_target_end,
            rna_seq=result.probe_l_target_rna,
            binding_dna=result.probe_l_binding_dna,
        )
        st.markdown("**Probe L full raw sequence (5'->3')**")
        st.code(probe_l_without_phosphate, language="text")
        if not add_probe_l_phosphate_tag:
            st.caption(selected_probe_l_variant)

    with probe_r_col:
        st.subheader("Probe R")
        st.markdown("**Recommended full sequence (5'->3')**")
        st.code(result.probe_r_full_display, language="text")
        _render_primary_note(
            [
                "Recommended display output. Final two displayed bases can be annotated as rN without changing the raw sequence.",
                (
                    "Current display option: last two bases shown as rN."
                    if format_probe_r_last_two_as_rna
                    else "Current display option: raw canonical bases shown."
                ),
            ]
        )
        st.markdown("---")
        _window_summary(
            label="Probe R",
            start=result.probe_r_target_start,
            end=result.probe_r_target_end,
            rna_seq=result.probe_r_target_rna,
            binding_dna=result.probe_r_binding_dna,
        )
        st.markdown("**Probe R full raw sequence (5'->3')**")
        st.code(result.probe_r_full_raw, language="text")

    primer_col, download_col = st.columns([1.5, 1])
    with primer_col:
        st.subheader("Universal Primers")
        st.markdown("**Universal primer F (5'->3')**")
        st.code(result.universal_primer_f, language="text")
        st.markdown("**Universal primer R (5'->3')**")
        st.code(result.universal_primer_r, language="text")

    with download_col:
        st.subheader("Downloads")
        st.download_button(
            label="Download TXT report",
            data=_result_to_txt_bytes(result),
            file_name="mirna_probe_design_report.txt",
            mime="text/plain",
            use_container_width=True,
        )
        st.download_button(
            label="Download CSV summary",
            data=_result_to_csv_bytes(result),
            file_name="mirna_probe_design_summary.csv",
            mime="text/csv",
            use_container_width=True,
        )

    if result.warnings:
        for warning in result.warnings:
            st.warning(warning)

    with st.expander("Design assumptions", expanded=False):
        for note in result.notes:
            st.markdown(f"- {note}")

    with st.expander("Copy-friendly text report", expanded=False):
        st.code(build_design_report(result), language="text")


if __name__ == "__main__":
    main()
