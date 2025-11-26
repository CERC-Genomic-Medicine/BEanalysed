# BEscreened : CRISPR Base Editing Screen Analysis Pipeline

A comprehensive toolkit for analyzing CRISPR Base Editing screens derived from crispr-BEasy.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Input Files](#input-files)
- [Pipeline Steps](#pipeline-steps)
  1. [Consolidating Files](#1-consolidating_filespy)
  2. [Repeat on Repeat](#2-repeatonrepeatpy)
  3. [RAUC Analysis](#3-raucpy)
  4. [Lollipop Plots](#4-bescreen_lollipop_plotpy)
  5. [Scatter Plots](#5-scatter_plot_bescreenpy-optional)
  6. [ChimeraX Visualization](#6-chimerax_scoringpy-and-chimeraxsh)

---

## Overview

This pipeline provides tools to analyze CRISPR Base Editing screens, from consolidating raw MaGeCK outputs with variant annotations to generating publication-ready visualizations and 3D structure mappings.

---


---

## Recommendation 

If you are working with files from an earlier version and want to annotate them with the current version (which is often preferable), run it without filtering (your original filtering will still be applied). As there have been changes between updates that can impact your results. 
___

## Requirements

### Software

- Python 3.10+
- UCSF ChimeraX

### Python Packages

See `requirements.txt` for complete list.

---

## Input Files

### Key File

A tab-delimited file with the following required columns:

| Column | Description |
|--------|-------------|
| `ID_guide` | Unique identifier for each guide |
| `Sequence` | Protospacer sequence |

### MaGeCK Gene Summary Files

Standard MaGeCK gene summary output files containing columns:

- `id`
- `pos|rank`, `neg|rank`
- `pos|p-value`, `neg|p-value`
- `pos|fdr`, `neg|fdr`
- `pos|lfc`

### MaGeCK sgRNA Summary Files

Tab-separated files with:

- An `sgrna` column containing guide identifiers
- Replicate identifiers encoded as `_r1`, `_r2`, etc. in the sgRNA names
- The analysis column (default: `LFC`)

**Example:**

| sgrna | LFC | score |
|-------|-----|-------|
| GENE1_guide1_r1 | -0.5 | 2.3 |
| GENE1_guide1_r2 | -0.4 | 2.1 |
| GENE1_guide2_r1 | 1.2 | 0.8 |
| GENE1_guide2_r2 | 1.1 | 0.9 |

### Annotation Excel

Excel file with VEP annotations including:

- `#Uploaded_variation`
- `Location`
- `Consequence`
- `Protein_position`
- `IMPACT`
- `Amino_acids`
- `Feature` (for isoform selection)
- `PICK` (if using `--Pick`)

And main annotation including:

- `ID`
- `protospacer`

### Domain BED File

Tab-delimited file with 4 columns:

| Column | Description |
|--------|-------------|
| `name` | Domain name |
| `start` | Start position |
| `end` | End position |
| `color` | Color (optional) |

### 3D Model

`.cif` file, usually obtained from [PDB](https://www.rcsb.org/)

### LiftOver File for 3D Model

Tab-delimited file with 4 columns:

| Column | Description |
|--------|-------------|
| `Protein` | Protein identifier |
| `Chain` | PDB chain identifier |
| `Start` | Start position |
| `End` | End position |

---

## Pipeline Steps

---

### 1. Consolidating_files.py

Consolidates MaGeCK results with VEP annotations into a single Excel summary.

#### Usage

```bash
python3 Consolidating_files.py -k FILE -I FILE [FILE ...] -a FILE [FILE ...] \
    -s str -X VARIABLE_NAMES [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-k`, `--key` | Tab-delimited key file containing columns `ID_guide` and `Sequence` (protospacer sequence) |
| `-I`, `--input` | One or more MaGeCK output files (per-gene results) |
| `-a`, `--annotation_excel` | Excel file(s) containing predicted variant consequences |
| `-s`, `--sheet` | Name of the sheet to use from the annotation Excel file (requires format `editor - XXXX`) |
| `-X`, `--xvar` | Comma-separated experimental conditions that should be reflected in the MaGeCK filenames (e.g., `UNT/TREAT,KO/WT`) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-l`, `--lib_sheet` | `Library` | Name of the library sheet in the annotation Excel file |
| `--Isoform` | — | Specific transcript isoform(s) to select for annotation |
| `--Pick` | `False` | Use VEP's PICK flag (column `PICK`) to select canonical annotations |
| `--e` | `False` | Mark sgRNAs with no predicted mutation as negative controls |
| `-n`, `--negative_Control` | — | List of protein/region names to designate as negative controls |
| `-p`, `--positive_Control` | — | List of protein/region names to designate as positive controls |
| `--out` | `summary` | Output filename (without extension) |

#### Experimental Condition Naming (`-X`)

The `-X` argument defines how experimental conditions are encoded in filenames. Use `/` to separate paired conditions and `,` to separate independent variables.

**Example:**

```bash
-X "UNT/TREAT,KO/WT"
```

This expects filenames containing exactly one of `UNT` or `TREAT` **and** exactly one of `KO` or `WT`.

**Valid filenames:**
- `results_UNT_KO.gene_summary.txt`
- `results_TREAT_WT.gene_summary.txt`

**Invalid filenames:**
- `results_UNT_TREAT.gene_summary.txt` (contains both conditions from a pair)

#### Output

The script generates an Excel file (`<output>.xlsx`) with:
- One sheet per experimental condition
- Columns including guide IDs, sequences, LFC values, p-values, FDR, consequence annotations, and control designations

#### Examples

**Basic Usage:**

```bash
python create_base_summary.py \
  -k guides_key.txt \
  -I mageck_UNT.gene_summary.txt mageck_TREAT.gene_summary.txt \
  -a vep_annotations.xlsx \
  -s "VEP Results" \
  -X "UNT/TREAT" \
  --out my_summary
```

**With Controls and Isoform Selection:**

```bash
python create_base_summary.py \
  -k guides_key.txt \
  -I results/*.gene_summary.txt \
  -a annotations.xlsx \
  -s "VEP Sheet" \
  -X "UNT/TREAT,KO/WT" \
  --Pick \
  -n AAVS1 ROSA26 \
  -p TP53 \
  --e \
  --out full_summary
```

**Using Specific Isoforms:**

```bash
python create_base_summary.py \
  -k guides_key.txt \
  -I results/*.gene_summary.txt \
  -a annotations.xlsx \
  -s "VEP Sheet" \
  -X "UNT/TREAT" \
  --Isoform ENST00000123456 ENST00000789012 \
  --out isoform_summary
```

#### Notes

- The `--Pick` and `--Isoform` options are mutually exclusive
- When using `--e`, guides targeting empty windows are automatically designated as negative controls
- Biological significance analysis requires negative controls to be defined

---

### 2. RepeatOnRepeat.py

Generates scatter plot grids comparing replicate measurements to assess reproducibility.

#### Description

This script generates scatter plot grids comparing replicate measurements (e.g., `_r1` vs `_r2`) from CRISPR screen data. It calculates Pearson correlations and produces publication-ready figures to assess replicate reproducibility.

#### Usage

```bash
python repeat_on_repeat.py -F <input_files> --xvar <conditions> [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--xvar` | Comma-separated experimental conditions encoded in filenames (e.g., `UNT/TREAT,KO/WT`) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-F`, `--files` | — | Input TSV files containing replicate data |
| `--var` | `LFC` | Column name to analyze (e.g., `LFC`, `score`) |
| `-o`, `--out` | `repeat_on_repeat.pdf` | Output image path |
| `--pivot` | `False` | Pivot/transpose the grid layout |

#### Output

The script generates a PDF file containing:
- A grid of scatter plots comparing all replicate pairs
- Pearson correlation coefficient (r) displayed on each plot
- A diagonal reference line (1:1) for visual comparison

#### Notes

- Replicate IDs are extracted from sgRNA names using the pattern `_r<number>` (e.g., `_r1`, `_r2`)
- Guides are matched by their base name (sgRNA name without replicate suffix)
- Only guides present in both replicates are included in the correlation

---

### 3. RAUC.py

Generates ROC-AUC curves for screen quality assessment.

#### Usage

```bash
python RAUC_excel.py -I EXCEL_FILE -V {pos,neg} [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-I`, `--input` | Excel file containing MaGeCK results (output from `create_base_summary.py`) |
| `-V`, `--value` | Direction to analyze: `pos` (positive selection) or `neg` (negative selection) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-R`, `--Remove_sheet` | — | Sheet name(s) to exclude from analysis |
| `--Gene_positive` | — | Gene name(s) to designate as positive controls |
| `-P`, `--Positive_consequence` | `["non-sense", "splice"]` | Consequence annotations to treat as positive controls |
| `--Gene_negative` | — | Gene name(s) to designate as negative controls |
| `-N`, `--Negative_consequence` | `["No predicted Mutation"]` | Consequence annotations to treat as negative controls |
| `-O`, `--out` | `RAUC.png` | Output image path |

#### Control Classification

**Positive Controls (Expected Hits):**

Guides classified as positive controls if any of the following are true:
- `Consequence` matches values in `-P` (default: `non-sense`, `splice`)
- `proteins` matches values in `--Gene_positive`
- `Controls` column equals `positive_control`

**Negative Controls (Expected Non-Hits):**

Guides classified as negative controls if any of the following are true:
- `Consequence` matches values in `-N` (default: `No predicted Mutation`)
- `proteins` matches values in `--Gene_negative`
- `Controls` column equals `negative_control`

#### Output

The script generates:

1. **ROC-AUC Plot** (`<out>.pdf`): Combined ROC curves for all sheets with AUC values
2. **Cumulative Rank Plots** (`<sheet>_<out>_rank_cumulative.pdf`): Per-sheet cumulative distribution of ranks by consequence category

---

### 4. BEscreen_lollipop_plot.py

Creates lollipop plots with protein domain annotations.

#### Description

This script creates publication-ready lollipop plots that visualize guide-level effects (log fold change) along a protein sequence. It overlays protein domain/feature annotations from a BED-like file and supports various visual encodings for statistical and biological significance.

#### Usage

```bash
python BEscreen_lollipop_plot_simplified.py -b BED_FILE -i EXCEL_FILE \
    --Bio_threshold BIOLOGICAL_THRESHOLD [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-b`, `--bed` | BED-like file containing protein features (start, end, name, protein) |
| `-i`, `--input` | Excel file from MaGeCK/VEP output (from `create_base_summary.py`) |
| `--Bio_threshold` | Biological significance threshold (two-sided quantile, e.g., `0.05`) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--stat_method` | `quantile` | Method for biological significance: `quantile`, `binom_sign`, or `sign_test` |
| `--scheme_location` | `top` | Position of domain scheme: `top`, `middle`, or `bottom` |
| `--histogram` | `False` | Add histogram of guide positions above the plot |
| `--violin` | `False` | Add violin/box plots showing LFC distribution by consequence |
| `--violin_detail` | `low` | Violin plot grouping: `low` (by consequence) or `high` (by detailed consequence) |
| `-p`, `--Prob_Threshold` | `1` | P-value threshold for marker size encoding |
| `--no_stem` | `False` | Remove stem lines from lollipop plot |
| `-F` | `False` | Use FDR instead of p-value for significance |
| `--highlight_region` | — | Highlight specific regions using `protein-feature` format |

#### Output

The script generates PDF files for each protein in each experimental condition:

```
BEscreen_<sheet_name>_<protein>.pdf
```

Each plot includes:
- Lollipop plot with LFC on y-axis and amino acid position on x-axis
- Protein domain/feature annotations
- Legend for consequences, p-values, and biological significance
- Optional histogram and/or violin plots

---

### 5. scatter_plot_BEscreen.py (Optional)

Generates scatter plot comparisons between paired experimental conditions.

#### Description

This script produces scatter plot grids comparing measurements (e.g., log fold change) between paired experimental conditions. It visualizes correlations with options for color-coding by consequence type, transparency by significance, and filtering by protein or consequence category.

#### Usage

```bash
python scatter_plot_BEscreen.py -I FILE -X string [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-I`, `--input` | Excel file containing MaGeCK/VEP results (output from `create_base_summary.py`) |
| `-X`, `--comparison` | Comparison to plot, encoded as paired conditions (e.g., `UNT/TREAT`) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-R`, `--Remove_sheet` | — | Sheet name(s) to exclude from analysis |
| `-C`, `--color` | `None` | Color points by attribute: `None` (all black) or `Consequence` |
| `-A`, `--alpha` | `None` | Transparency mode: `None` (opaque), `Significance`, or `all_little_bit` |
| `-B`, `--biological_line` | — | Plot biological significance threshold lines (quantile, e.g., `0.05`) |
| `-S`, `--Significance_threshold` | — | P-value threshold for significance-based transparency |
| `-e`, `--elements_to_plot` | `all` | Filter elements: `all` or `coding_only` |
| `-F`, `--Square_Format` | `False` | Force square 1:1 aspect ratio with matched axis ranges |
| `-P`, `--protein` | — | Protein name(s) to filter and plot |
| `-V`, `--Variable` | `lfc` | Column to plot (e.g., `lfc`, `p_value`) |
| `-O`, `--Output` | — | Output filename (without extension) |

#### Output

The script generates a PDF file containing:
- Scatter plot grid comparing paired conditions
- Pearson correlation coefficient (r) on each plot
- Optional biological significance threshold lines
- Legend for consequence types and significance

**Default naming:**
- With `-O`: `<output>.pdf`
- Without `-O`: `<comparison>_<proteins>_<variable>.pdf`

---

### 6. ChimeraX_scoring.py and chimeraX.sh

Maps screen results to 3D protein structures for visualization in ChimeraX.

#### Description

This script produces `.defattr` attribute files compatible with UCSF ChimeraX, enabling visualization of CRISPR screen results (e.g., log fold change) directly on protein 3D structures. It maps guide-level scores to residue positions and handles coordinate system conversion between screen data and PDB structures.

#### Usage

```bash
python ChimeraX_scoring.py -i FILE -b FILE [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Excel file containing MaGeCK/VEP results (output from `create_base_summary.py`) |
| `-b`, `--lift_over` | Tab-delimited file mapping protein positions to PDB chain coordinates |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--Bio_threshold` | `1` | Biological significance p-value threshold (two-sided) |
| `--duplicate_strategy` | `median` | Strategy for multiple guides at same position: `max`, `mean`, or `median` |
| `-p`, `--Prob_Threshold` | `1` | P-value threshold for inclusion |
| `-R`, `--Remove_sheet` | — | Sheet name(s) to exclude from processing |

#### Output Files

For each sheet in the Excel file, the script generates two `.defattr` files:

**Significant Residues (`<sheet>_Sig.defattr`):**

Contains LFC values for residues meeting significance thresholds:
- Biological p-value ≤ `--Bio_threshold`
- Statistical p-value ≤ `-p`

**Non-Significant Residues (`<sheet>_nonSig.defattr`):**

Contains a flag (value = 1) for residues with data but not meeting significance criteria.

#### Duplicate Handling

When multiple guides target the same residue position, the `--duplicate_strategy` determines the reported value:

| Strategy | Description |
|----------|-------------|
| `median` | Median LFC across guides (default, robust to outliers) |
| `mean` | Average LFC across guides |
| `max` | Maximum absolute LFC (most extreme effect) |

#### ChimeraX Visualization

```bash
bash chimeraX.sh <your_3D_model> <your_attribute_file>
```

---

## Author

Vincent Chapdelaine (vincent.chapdelaine@mcgill.ca)

## Version

1.1 (2024)
