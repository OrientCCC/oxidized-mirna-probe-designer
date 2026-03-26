# Rule-Based miRNA Probe Designer

This project implements a small, transparent Python tool for designing miRNA Probe L / Probe R oligos from:

- a miRNA sequence in RNA 5'->3' orientation
- an oxidized guanine position using a 1-based index
- a configurable split rule

It includes a lightweight Streamlit app and a small test suite.

## What the tool outputs

For each design request, the tool returns:

- the Probe R target-binding region
- the Probe L target-binding region
- the full Probe R sequence after adding its fixed universal adaptor
- the full Probe L sequence after adding its fixed universal adaptor
- the fixed universal primer F and universal primer R
- a text summary of sequence orientation, coordinates, assumptions, notes, and warnings

## Fixed sequences

- Probe L universal adaptor: `CTCTATGGGCAGTCGGTGATATCGACCTGTACCTCT`
- Probe R universal adaptor: `CCATCTCATCCCTGCGTGTCCATCATCCGTACGTGT`
- Universal primer F: `CCATCTCATCCCTGCGTGTC`
- Universal primer R: `AGAGGTACAGGTCGATATCA`

## Orientation conventions

- Input miRNA is treated as RNA written 5'->3'.
- The ligated single-stranded probe product is defined as `5' - Probe R - Probe L - 3'`.
- Probe target-binding regions are output as DNA oligos written 5'->3', so each binding region is the DNA reverse complement of its assigned miRNA target window.

## Split rule

Version 1 makes the split rule explicit and configurable instead of hardcoding it.

Supported modes:

- `oxidized_base_belongs_to_R`
- `oxidized_base_belongs_to_L`

Version 1 defaults to `oxidized_base_belongs_to_R` in both the core API and the Streamlit UI.

In the current version, the probes together cover the full miRNA with no gap. The oxidized base is included in exactly one probe depending on the selected split mode:

- If the oxidized base belongs to Probe R, Probe L covers positions `1..p-1` and Probe R covers positions `p..N`.
- If the oxidized base belongs to Probe L, Probe L covers positions `1..p` and Probe R covers positions `p+1..N`.

The only split decision is which probe contains the oxidized base.

Designs are rejected if either probe would end up with zero miRNA-complementary bases. For example:

- `oxidized_pos = 1` with `oxidized_base_belongs_to_R` is rejected because Probe L would be empty.
- `oxidized_pos = N` with `oxidized_base_belongs_to_L` is rejected because Probe R would be empty.

## Probe assembly

Version 1 assembles probes as:

- Probe R full sequence = `[Probe R adaptor] + [Probe R target-binding region]`
- Probe L full sequence = `[Probe L target-binding region] + [Probe L adaptor]`

Probe L can optionally be rendered with a literal ordering prefix:

- enabled: `/PO4/[Probe L target-binding region][Probe L adaptor]`
- disabled: `[Probe L target-binding region][Probe L adaptor]`

## Probe R `rN` display formatting

Probe R supports an optional display-only formatting mode for the final two nucleotides of the full probe sequence.

Examples:

- `...AC -> ...rArC`
- `...GC -> ...rGrC`

This does not change the underlying designed sequence. The result object includes both:

- `probe_r_full_raw`
- `probe_r_full_display`

Version 1 defaults to enabling this display formatting in both the core API and the Streamlit UI.

## Configurable assumptions in version 1

- split mode
- Probe L `/PO4/` ordering tag
- Probe R final-two-base `rN` display annotation

## Biological rules not modeled in version 1

This first version is intentionally simple and rule-based. It does not model:

- Tm optimization
- mismatch scoring
- ligation efficiency
- secondary structure
- off-target analysis
- empirical chemistry constraints beyond the explicit rules above

## Installation

Use Python 3.11+.

Recommended with `mamba`:

```bash
mamba env create -f environment.yml
mamba activate mirna-probe-designer
```

If you prefer to install into an existing environment, you can still use:

```bash
mamba install -c conda-forge streamlit
```

## Run the Streamlit app

```bash
mamba activate mirna-probe-designer
streamlit run app.py
```

## Run the tests

```bash
mamba activate mirna-probe-designer
python -m unittest discover -s tests -v
```

## Core API

The main function is:

```python
from probe_designer import design_probes

result = design_probes(
    mirna_seq="UGAGGUAGUAGGUUGUAUAGUU",
    oxidized_pos=8,
    split_mode="oxidized_base_belongs_to_R",
    add_probe_l_phosphate_tag=True,
    format_probe_r_last_two_as_rna=True,
)
```

The returned `DesignResult` includes cleaned inputs, target-window coordinates, RNA windows, DNA binding sequences, full probes, universal primers, notes, warnings, and a text schematic.
