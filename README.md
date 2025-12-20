# BEanalysed: CRISPR Base Editing Screen Analysis Pipeline

A comprehensive toolkit for analyzing CRISPR Base Editing screens derived from crispr-BEasy.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Input Files](#input-files)
- [Pipeline Steps](#pipeline-steps)
  1. [Consolidating Files](#1-consolidating_filespy)
  2. [Repeat on Repeat](#2-repeatonrepeatpy)
  3. [RAUC Analysis](#3-BE_RAUCedpy)
  4. [RAUC Analysis per Gene](#4-BE_rauc_PerGenepy)
  5. [Lollipop Plots](#5-BEanalyzed_lollipop_plotpy)
  6. [Scatter Plots](#6-scatter_plot_BEanalysedpy-optional)
  7. [ChimeraX Visualization](#7-chimerax_scoringpy-and-chimeraxsh)

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

A tab/comma/space-delimited file with the following required columns:

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

### Annotation Excel files (s)

Excel file(s) (`.xlsx`, `.xls`, `.xlsm`) with VEP annotations including:

- `#Uploaded_variation`
- `Location`
- `Consequence`
- `Protein_position`
- `IMPACT`
- `Amino_acids`
- `Feature` (for isoform selection)
- `PICK` (if using `--Pick`)
- `Replicate` or column named after the editor (derived from the sheet name `"editor - YourEditorName"`)

And main annotation including:

- `ID`
- `protospacer`

**Note:** Files must be valid Excel format. If you have a `.zip` file, extract it first before running the pipeline.

### Domain BED-like File

Tab-delimited file with 4 columns (5 with optional color) :


| Column | Description |
|--------|-------------|
| `proteins` | Proteins asociated with the domains |
| `name` | Domain name |
| `start` | Start position |
| `end` | End position |
| `color` | Color (optional) |

**Note:** this file may contain multiple proteins (all of which will be represented on separate figures).
**Note:** The File must also contain a segment with the end of the protein (Domain name can be None), otherwise the protein's limit are infered from the maximum of the domains and the represented mutations.

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


## Caveat about current behaviour

```diff
- Important!
Current version requieres Isoform filtering of ALL annotations (including controls negative and positive).
 This may change in the future. Intergenic variants generally holds a non discript feature class '-' that may be used as Isoform.
 Also, as the guides with no predicted mutation are not annotated (since they do not induce a mutation) they pass this filter no matter where they are.
 **Please verify that your controls are indeed treated as you wish them to be treated**. (see Isoform selection section below] )
```

#### Usage

```bash
python3 Consolidating_files.py -k FILE -I FILE [FILE ...] -a FILE [FILE ...] \
    -s str -X VARIABLE_NAMES [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-k`, `--key` | Tab/comma/space-delimited key file containing columns `ID_guide` and `Sequence` (protospacer sequence) |
| `-I`, `--input` | One or more MaGeCK output files (per-gene results)  **space delimited**  |
| `-a`, `--annotation_excel` | Excel file(s) containing predicted variant consequences  **space delimited**  |
| `-s`, `--sheet` | Name of the sheet to use from the annotation Excel file (naming consideration see [Empty Sheet Detection](#empty-sheet-detection) )  |
| `-X`, `--xvar` | Comma-separated experimental conditions that should be reflected in the MaGeCK filenames (e.g., `UNT/TREAT,KO/WT`) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-e`, `--empty_sheet` | Auto-detect | Sheet name for guides without predicted mutation (see [Empty Sheet Detection](#empty-sheet-detection)) |
| `-l`, `--lib_sheet` | `Library` | Name of the library sheet in the annotation Excel file |
| `--Isoform` | — | Specific transcript isoform(s) to select for annotation including Positive and Negative controls (if used) |
| `--Pick` | `False` | Use VEP's PICK flag (column `PICK`) to select canonical annotations |
| `--Empty_controls` | `False` | Mark sgRNAs with no predicted mutation as negative controls (even in target genes) |
| `-F` | `False` | Bypass some validation checks (force mode) |
| `--Negative_Control_Genes` | — |  **space delimited**  List of protein/region names to designate as negative controls |
| `--Positive_Control_Genes` | — |  **space delimited**  List of protein/region names to designate as positive controls |
| `--Negative_Control_Consequences` | — |  **space delimited**  List of consequence types to designate as negative controls |
| `--Positive_Control_Consequences` | — |  **space delimited**  List of consequence types to designate as positive controls |
| `--out` | `summary` | Output filename (without extension) |

#### Empty Sheet Detection

The `-e/--empty_sheet` option specifies which sheet contains guides without predicted mutations (empty editing windows).

**Auto-detection logic** (when `-e` not provided):

1. Attempts to replace `"editor"` with `"no_mutation"` in the annotation sheet name
   - Example: `"editor - YourEditorName"` → `"no_mutation - YourEditorName"`
2. If the auto-detected sheet exists, it is used
3. If the auto-detected sheet does not exist, or if `"editor"` was not in the sheet name, guides missing from VEP annotations are inferred as having no predicted mutation

** Reminder '' allows a string of character including spaces to be used as a single argument.

**Explicit specification:**

```bash
-e "VEP_no_mutation"
```

When `-e` is explicitly specified and the sheet is not found, an error is raised (no fallback).

#### Replicate Window

The script will also output information considering replicate windows (i.e. which guides offer predicted replicate window).
**Auto-detection logic**
1. column named 'Replicate'
2. Column named after you editor, derived from sheet name with `"editor - YourEditorName"` (including spaces)

#### Isoform selection

Isoform selection is requiered either upstream by user or through the --Isoform or --Pick options. Isoform correspond to the Feature column in the VEP file (or CRISPR-BEasy annotation), Pick correspond to the column of the same name. This is requiered to categorise which consequence the mutation will hold in the consolidated file (and downstream analysis). In case of intergenic variant '-' may be used as an Isoform.

#### Validation

The script performs strict validation to ensure data completeness:

1. **Guide coverage:** All guides in the key file must be present in the library annotation sheet (unless `-F` is used)
2. **Annotation completeness:** All guides must be present in either the VEP annotation sheet or the empty guides sheet. Missing guides will raise an error:
   ```
   Error: The following 2 guide(s) are missing from both VEP annotations 
   ("VEP editor Replicate") and the empty guides sheet ("VEP no_mutation Replicate"):
   ['BRCA1_001', 'TP53_002']
   ```
3. **File format:** Only valid Excel files are accepted as annotation file.

#### Experimental Condition Naming (`-X`)

The `-X` argument defines how experimental conditions are encoded in filenames. Use `/` to separate paired conditions and `,` to separate independent variables.

**Example:**

```bash
-X "UNT/TREAT,KO/WT"
```

This expects filenames containing exactly one of `UNT` or `TREAT` **and** exactly one of `KO` or `WT' per 'UNT' and 'TREAT'

**Valid filenames:**
- `results_UNT_KO.gene_summary.txt`
- `results_TREAT_WT.gene_summary.txt`
- `results_UNT_WT.gene_summary.txt`
- `results_TREAT_KO.gene_summary.txt`

**Additionnal Consideration**
Condition names should not contain one another. such as TREATED/UNTREATED as TREATED is found in UNTREATED.

**Invalid filenames:**
- `results_UNT_TREAT.gene_summary.txt` (contains both conditions from a pair)

#### Output

The script generates an Excel file (`<out>.xlsx`) with:
- One sheet per experimental condition
- Columns including guide IDs, sequences, LFC values, p-values, FDR, consequence annotations, and control designations

#### Examples

**Basic Usage:**

```bash
python Consolidating_files.py \
  -k guides_key.txt \
  -I mageck_UNT.gene_summary.txt mageck_TREAT.gene_summary.txt \
  -a vep_annotations.xlsx \
  -s "VEP editor Replicate" \
  -X "UNT/TREAT" \
  --out my_summary
```

**With Explicit Empty Sheet:**

```bash
python Consolidating_files.py \
  -k guides_key.txt \
  -I mageck_UNT.gene_summary.txt mageck_TREAT.gene_summary.txt \
  -a vep_annotations.xlsx \
  -s "VEP editor Replicate" \
  -e "VEP no_mutation Replicate" \
  -X "UNT/TREAT" \
  --out my_summary
```

**With Controls and Isoform Selection:**

```bash
python Consolidating_files.py \
  -k guides_key.txt \
  -I results/*.gene_summary.txt \
  -a annotations.xlsx \
  -s "VEP editor Replicate" \
  -X "UNT/TREAT,KO/WT" \
  --Pick \
  --Negative_Control_Genes AAVS1 ROSA26 \
  --Positive_Control_Genes TP53 \
  --Empty_controls \
  --out full_summary
```

**Using Specific Isoforms:**

```bash
python Consolidating_files.py \
  -k guides_key.txt \
  -I results/*.gene_summary.txt \
  -a annotations.xlsx \
  -s "VEP editor Replicate" \
  -X "UNT/TREAT" \
  --Isoform ENST00000123456 ENST00000789012 \
  --out isoform_summary
```

#### Python Module Usage

The script can also be imported and used as a Python module:

```python
from Consolidating_files import (
    BaseSummaryGenerator,
    ControlConfig
)

# Configure controls
controls = ControlConfig(
    negative_genes={"AAVS1", "NonTargeting"},
    positive_genes={"Essential1"},
    empty_as_negative=True
)

# Create generator
generator = BaseSummaryGenerator(
    key_file="guides.csv",
    annotation_files=["annotations.xlsx"],
    annotation_sheet="VEP editor Replicate",
    library_sheet="Library",
    empty_sheet="VEP no_mutation Replicate",  # Optional
    control_config=controls
)

# Generate summary
generator.generate_summary(
    mageck_files=["mageck_TREAT_KO.txt", "mageck_UNT_WT.txt"],
    treatment_spec="TREAT/UNT,KO/WT",
    output_path="summary.xlsx"
)
```

#### Notes

- The `--Pick` and `--Isoform` options are mutually exclusive
- When using `--Empty_controls`, guides targeting empty windows are automatically designated as negative controls
- Biological significance analysis requires negative controls to be defined
- All guides must be accounted for in either the VEP annotation sheet or the empty guides sheet

---

### 2. RepeatOnRepeat.py

Generates scatter plot grids comparing replicate measurements to assess reproducibility.

#### Description

This script creates scatter plot grids comparing log fold change (or other metrics) between biological/technical replicates. It calculates Pearson correlation coefficients and visualizes the reproducibility of the screen.

#### Usage

```bash
python RepeatOnRepeat.py -F FILE [FILE ...] [options]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-F`, `--files` | — | Input TSV files containing replicate data ( **space delimited** ) |
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

### 3. BE_ROC.py

Generates ROC-AUC curves for screen quality assessment and a representation of the cumulative fraction per rank for each `Consequence_Detail`.

#### Usage

```bash
python BE_ROC.py -I EXCEL_FILE -V {pos,neg} [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-I`, `--input` | Excel file containing MaGeCK results (output from `Consolidating_files.py`) |
| `-V`, `--value` | Direction to analyze: `pos` (positive selection) or `neg` (negative selection) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-R`, `--Remove_sheet` | — | Sheet name(s) to exclude from analysis |
| `-O`, `--out` | `ROC` | Output image path **Do not include extensions** |

#### Control Classification

Controls are as defined in the consolidation script.

#### Output

The script generates:

1. **ROC-AUC Plot** (`<out>.pdf`): Combined ROC curves for all sheets with AUC values
2. **Cumulative Rank Plots** (`<sheet>_<out>_rank_cumulative.pdf`): Per-sheet cumulative distribution of ranks by consequence category

---


### 4. BE_ROC_PerGene.py


Generates ROC-AUC curves for screen quality assessment based on the predicted consequences within each gene independently. With a second figure detailling the `Consequence Detail` ROC-AUC breakdown per gene x analysis combinations.

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-I`, `--input` | Excel file containing MaGeCK results (output from `Consolidating_files.py`) |
| `-V`, `--value` | Direction to analyze: `pos` (positive selection) or `neg` (negative selection) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-R`, `--Remove_sheet` | — | Sheet name(s) to exclude from analysis ( **space delimited** ) |
| `-O`, `--out` | `ROC_per_Gene` | Output images path ( **Do not include extension.** )  |
| `-G`, `--Genes` | `` (all genes) | **space delimited** list of Gene names |
|`--Negative_Control_Consequences`|`No predicted Mutation`| **space delimited** list of negative control consequences (accepted 'splice','non-sense','missense','synonymous','non-coding','regulatory', 'No predicted Mutation','N/A')|
|`--Positive_Control_Consequences`|`'splice' 'non-sense'`| **space delimited** list of negative control consequences (accepted 'splice','non-sense','missense','synonymous','non-coding','regulatory', 'No predicted Mutation','N/A')|


### 5. BEanalyzed_lollipop_plot.py

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
| `-i`, `--input` | Excel file from MaGeCK/VEP output (from `Consolidating_files.py`) |
| `--Bio_threshold` | Biological significance threshold (two-sided quantile, e.g., `0.05`) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--stat_method` | `quantile` | Method for biological significance: `quantile`, `binom_sign`, or `sign_test` |
| `--scheme_location` | `top` | Position of domain scheme: `top` or `bottom` |
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
BEscreened_<sheet_name>_<protein>.pdf
```

Each plot includes:
- Lollipop plot with LFC on y-axis and amino acid position on x-axis
- Protein domain/feature annotations
- Legend for consequences, p-values, and biological significance
- Optional histogram and/or violin plots

---

### 6. scatter_plot_BEanalysed.py (Optional)

Generates scatter plot comparisons between paired experimental conditions. Considered optional as the exact comparison made and choice is highly dependant on the biological question ( Example : technical assessement spCas9 vs spCas9-NG or biological question WT vs KO)

#### Description

This script produces scatter plot grids comparing measurements (e.g., log fold change) between paired experimental conditions. It visualizes correlations with options for color-coding by consequence type, transparency by significance, and filtering by protein or consequence category.

#### Usage

```bash
python scatter_plot_BEscreen.py -I FILE -X string [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `-I`, `--input` | Excel file containing MaGeCK/VEP results (output from `Consolidating_files.py`) |
| `-X`, `--comparison` | Comparison to plot, encoded as paired conditions (e.g., `UNT/TREAT`) |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-R`, `--Remove_sheet` | — | Sheet name(s) to exclude from analysis ( **space delimited** ) |
| `-C`, `--color` | `None` | Color points by attribute: `None` (all black) or `Consequence` |
| `-A`, `--alpha` | `None` | Transparency mode: `None` (opaque), `Significance`, or `all_little_bit` |
| `-B`, `--biological_line` | — | Plot biological significance threshold lines (quantile, e.g., `0.05`) |
| `-S`, `--Significance_threshold` | — | P-value threshold for significance-based transparency |
| `-e`, `--elements_to_plot` | `all` | Filter elements: `all` or `coding_only` |
| `-F`, `--Square_Format` | `False` | Force square 1:1 aspect ratio with matched axis ranges |
| `-P`, `--protein` | — | Protein name(s) to filter and plot ( **space delimited** ) |
| `-V`, `--Variable` | `lfc` | Column to plot (e.g., `lfc`, `p_value`) |
| `-O`, `--Output` | — | Output filename (without extension) |

#### Output

The script generates a PDF file containing:
- Scatter plot grid comparing paired conditions
- Pearson correlation coefficient (r) on each plot
- Optional biological significance threshold lines
- Legend for consequence types and significance

**Default naming:**
- With `-O`: `<o>.pdf`
- Without `-O`: `<comparison>_<proteins>_<variable>.pdf`

---

### 7. ChimeraX_scoring.py and chimeraX.sh

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
| `-i`, `--input` | Excel file containing MaGeCK/VEP results (output from `Consolidating_files.py`) |
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

## Error Messages

The pipeline provides informative error messages for common issues:

| Error | Cause | Solution |
|-------|-------|----------|
| `Unsupported file format ".zip"` | Annotation file is a zip archive | Extract the `.xlsx` file from the archive |
| `Annotation file not found` | File path is incorrect | Check the file path |
| `Sheet "X" not found` | Sheet name doesn't exist | Check available sheets listed in the error |
| `guide(s) are missing from both VEP annotations and empty sheet` | Guides not annotated | Ensure all guides have VEP annotations or are listed in the empty sheet |
| `Specified empty sheet "X" not found` | Explicit `-e` sheet doesn't exist | Check sheet name or remove `-e` to use auto-detection |
| Extra arguments | space within an argument | specify the arguments with '' |

---

## Author

Vincent Chapdelaine (vincent.chapdelaine@mcgill.ca)

## Version

2.0 (2024)
