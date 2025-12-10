#!/usr/bin/env python3
"""
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 2.0
YEAR: 2024

Create Base Summary for CRISPR Base Editor Screen Analysis

This module provides functionality to create concatenated summaries for
analyses performed with MaGeCK, combining guide information with VEP
annotations and control designations.

Can be used as a command-line tool or imported as a module.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import warnings
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path
from typing import Iterator, Optional, Sequence, Union

import numpy as np
import pandas as pd


# =============================================================================
# CONSTANTS & CONFIGURATION
# =============================================================================

VALID_CONSEQUENCE_CATEGORIES = frozenset([
    'splice', 'non-sense', 'missense', 'synonymous',
    'non-coding', 'regulatory', 'No predicted Mutation', 'N/A'
])

VARIANT_CONSEQUENCES_MAPPING = {
    'missense_variant': 'missense',
    'intron_variant': 'non_coding',
    'downstream_gene_variant': 'non_coding',
    'NMD_transcript_variant': 'non-sense',
    'upstream_gene_variant': 'non_coding',
    '3_prime_UTR_variant': 'non_coding',
    'synonymous_variant': 'synonymous',
    'non_coding_transcript_exon_variant': 'non_coding',
    'splice_region_variant': 'non_coding',
    'splice_polypyrimidine_tract_variant': 'non_coding',
    'stop_gained': 'non-sense',
    'coding_sequence_variant': 'coding',
    '5_prime_UTR_variant': 'non_coding',
    'regulatory_region_variant': 'regulatory',
    'splice_donor_variant': 'splice',
    'splice_acceptor_variant': 'splice',
    'non_coding_transcript_variant': 'non_coding',
    'splice_donor_region_variant': 'non_coding',
    'splice_donor_5th_base_variant': 'non_coding',
    'TF_binding_site_variant': 'regulatory',
    'start_lost': 'non-sense',
    'stop_lost' : 'non-sense' ,
    'stop_retained_variant' : 'stop_retained_variant', # Will be explored below (can be missense or synonymous)
    'incomplete_terminal_codon_variant': 'non-sense',
    'NA/NotFound': 'none',
    'None': 'none',
    'No predicted Mutation': 'No predicted Mutation'
}

# Consequence priority for classification (higher index = higher priority)
CONSEQUENCE_PRIORITY = ['non_coding', 'regulatory', 'synonymous', 'missense', 'non-sense', 'splice']


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ControlConfig:
    """Configuration for positive and negative controls."""
    negative_genes: Optional[set[str]] = None
    positive_genes: Optional[set[str]] = None
    negative_consequences: Optional[set[str]] = None
    positive_consequences: Optional[set[str]] = None
    empty_as_negative: bool = False

    def validate(self, force: bool = False) -> None:
        """Validate control configuration for conflicts."""
        if self.negative_genes and self.positive_genes:
            overlap = self.negative_genes & self.positive_genes
            if overlap:
                raise ValueError(f'Control Genes cannot overlap: {overlap}')

        if self.negative_consequences and self.positive_consequences and not force:
            overlap = self.negative_consequences & self.positive_consequences
            if overlap:
                raise ValueError(
                    f'Control consequences should not overlap: {overlap}\n'
                    'Use force=True to override'
                )


@dataclass
class AnnotationResult:
    """Result of parsing a single VEP annotation."""
    guide_id: str
    consequence: Optional[str]
    consequence_detail: str
    impact: Optional[str]
    amino_acid_position: Optional[int]
    mutations: str
    replicates: str


@dataclass
class GuideAnnotation:
    """Complete annotation for a guide RNA."""
    guide_id: str
    sequence: str
    protein: str
    consequence: Optional[str] = None
    consequence_detail: Optional[str] = None
    impact: Optional[str] = None
    position: Optional[int] = None
    mutations: Optional[str] = None
    replicates: Optional[str] = None
    control_type: Optional[str] = None


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def safe_int(value: str) -> Optional[int]:
    """Safely convert a string to integer, returning None on failure."""
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def parse_treatment_levels(xvar: str) -> list[list[str]]:
    """
    Parse treatment variable specification into sorted level groups.
    
    Args:
        xvar: Comma-separated groups of slash-separated levels
              (e.g., "UNT/TREAT,KO/WT")
    
    Returns:
        List of level groups, each sorted by length (descending)
    """
    return [
        sorted(group.split('/'), key=len, reverse=True)
        for group in xvar.split(',')
    ]


def extract_labels_from_filename(filename: str, levels: list[list[str]]) -> list[str]:
    """
    Extract treatment labels from filename based on level specifications.
    
    Args:
        filename: Path to file
        levels: Parsed treatment levels from parse_treatment_levels()
    
    Returns:
        List of matched labels from the filename
    """
    fname = os.path.basename(filename)
    labels = []
    
    for group in product(*levels):
        hits = [lvl for lvl in group if re.search(rf'{re.escape(lvl)}', fname)]
        if len(hits) == len(group):
            labels.extend(hits)
    
    return labels


def read_delimited_file(
    filepath: Union[str, Path],
    usecols: Optional[list[str]] = None
) -> pd.DataFrame:
    """
    Read a delimited file with automatic delimiter detection.
    
    Args:
        filepath: Path to the file
        usecols: Columns to read (optional)
    
    Returns:
        DataFrame with file contents
    
    Raises:
        ValueError: If file cannot be parsed
    """
    try:
        return pd.read_csv(filepath, sep=None, engine='python', usecols=usecols)
    except Exception as e:
        raise ValueError(
            f'Could not parse {filepath}. Expected tab, comma, or space '
            f'delimited file with columns: {usecols}'
        ) from e


def read_excel_sheets(
    excel_files: Sequence[Union[str, Path]],
    sheet_name: str,
    dtype: type = str
) -> pd.DataFrame:
    """
    Read and concatenate a sheet from multiple Excel files.
    
    Args:
        excel_files: Sequence of Excel file paths
        sheet_name: Name of sheet to read
        dtype: Data type for columns
    
    Returns:
        Concatenated DataFrame from all files
    """
    dfs = []
    for excel_file in excel_files:
        df = pd.read_excel(excel_file, sheet_name=sheet_name, dtype=dtype)
        dfs.append(df.fillna(''))
    return pd.concat(dfs, axis=0, ignore_index=True)


# =============================================================================
# CORE CLASSES
# =============================================================================

class KeyFileHandler:
    """Handles loading and mapping of guide key files."""
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize with a key file.
        
        Args:
            filepath: Path to key file with ID_guide and Sequence columns
        """
        self.filepath = Path(filepath)
        self._df = read_delimited_file(filepath, usecols=['ID_guide', 'Sequence'])
        
        self.sequence_to_id = dict(zip(self._df['Sequence'], self._df['ID_guide']))
        self.id_to_sequence = dict(zip(self._df['ID_guide'], self._df['Sequence']))
    
    @property
    def sequences(self) -> set[str]:
        """Set of all guide sequences."""
        return set(self._df['Sequence'])
    
    @property
    def guide_ids(self) -> set[str]:
        """Set of all guide IDs."""
        return set(self._df['ID_guide'])
    
    def get_id(self, sequence: str) -> Optional[str]:
        """Get guide ID for a sequence."""
        return self.sequence_to_id.get(sequence)
    
    def get_sequence(self, guide_id: str) -> Optional[str]:
        """Get sequence for a guide ID."""
        return self.id_to_sequence.get(guide_id)


class MaGeCKParser:
    """Parser for MaGeCK output files."""
    
    @staticmethod
    def parse(filepath: Union[str, Path], value_column: Optional[str] = None) -> pd.DataFrame:
        """
        Parse a MaGeCK output file.
        
        Args:
            filepath: Path to MaGeCK file
            value_column: If specified, only extract this column with id
        
        Returns:
            DataFrame with parsed MaGeCK data
        """
        mage = pd.read_csv(filepath, sep='\t', header=0)
        
        if value_column:
            return pd.DataFrame({
                'id': mage['id'],
                value_column: mage[value_column]
            })
        
        result = pd.DataFrame({
            'id': mage['id'],
            'lfc': mage['pos|lfc'],
            'pos|rank': mage['pos|rank'],
            'neg|rank': mage['neg|rank'],
            'pos|p-value': mage['pos|p-value'],
            'neg|p-value': mage['neg|p-value'],
            'pos|fdr': mage['pos|fdr'],
            'neg|fdr': mage['neg|fdr'],
        })
        
        # Select p-value and FDR based on LFC direction
        result['p_value'] = np.where(
            mage['pos|lfc'] < 0,
            mage['neg|p-value'],
            mage['pos|p-value']
        )
        result['fdr'] = np.where(
            mage['pos|lfc'] < 0,
            mage['neg|fdr'],
            mage['pos|fdr']
        )
        
        return result


class IsoformFilter:
    """Handles isoform-based filtering of annotations."""
    
    def __init__(self, annotations_df: pd.DataFrame):
        """
        Initialize with annotations DataFrame.
        
        Args:
            annotations_df: DataFrame with Feature, Gene, and ID columns
        """
        self.df = annotations_df
        self._build_mappings()
    
    def _build_mappings(self) -> None:
        """Build feature-gene and gene-feature mappings."""
        self.feature_to_genes = self.df.groupby('Feature')['Gene'].apply(set).to_dict()
        self.gene_to_features = self.df.groupby('Gene')['Feature'].apply(
            lambda x: list(set(x))
        ).to_dict()
    
    def filter_by_isoforms(self, isoforms: Sequence[str]) -> pd.DataFrame:
        """
        Filter annotations to specified isoforms.
        
        Args:
            isoforms: List of isoform/feature IDs to keep
        
        Returns:
            Filtered DataFrame
        
        Raises:
            ValueError: If isoforms not found or create conflicts
        """
        isoform_set = set(isoforms)
        
        # Check for missing isoforms
        missing = isoform_set - set(self.df['Feature'])
        if missing:
            raise ValueError(f"Isoforms not found in annotations: {missing}")
        
        # Check for gene conflicts
        genes_selected = pd.Series([self.feature_to_genes[i] for i in isoforms])
        conflicts = genes_selected[genes_selected.duplicated(keep='first')]
        
        if len(conflicts) > 0:
            error_lines = []
            for gene in conflicts:
                features = ' '.join([
                    f for f in self.gene_to_features.get(gene, [])
                    if f in isoform_set
                ])
                error_lines.append(f"  Gene '{gene}' is shared by features: {features}")
            raise ValueError(
                "Multiple features map to the same gene:\n" + "\n".join(error_lines)
            )
        
        # Apply filter
        filtered = self.df[self.df['Feature'].isin(isoform_set)].copy()
        
        # Report dropped annotations
        not_taken = self.df[~self.df['ID'].isin(set(filtered['ID']))]
        if not not_taken.empty:
            print(f'Warning: {len(set(not_taken["ID"]))} guides eliminated due to isoform selection')
            print(f'Eliminated guides: {set(not_taken["ID"])}')
            
            orphan_genes = set(not_taken['Gene']) - set(self.gene_to_features.keys())
            if orphan_genes:
                print(f'Guides for these genes were not selected: {orphan_genes}')
        
        # Verify no ambiguity remains
        if len(filtered['ID']) != len(set(filtered['ID'])):
            ambiguous = filtered[filtered.duplicated(subset=['ID'], keep=False)]
            raise ValueError(
                f"Annotation remains ambiguous after isoform selection: {ambiguous}"
            )
        
        return filtered
    
    def filter_by_pick(self) -> pd.DataFrame:
        """Filter annotations using VEP's PICK column."""
        return self.df[self.df['PICK'] == '1'].copy()


class VEPAnnotationParser:
    """Parser for VEP annotation Excel files."""
    
    def __init__(
        self,
        consequence_mapping: dict[str, str] = None
    ):
        """
        Initialize parser.
        
        Args:
            consequence_mapping: Mapping from VEP consequences to categories
        """
        self.consequence_mapping = consequence_mapping or VARIANT_CONSEQUENCES_MAPPING
    
    def parse(
        self,
        excel_files: Sequence[Union[str, Path]],
        sheet_name: str,
        id_mapping: dict[str, str],
        isoform_selection: Optional[Union[str, list[str]]] = None,
        empty_sheet: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Parse VEP annotations from Excel files.
        
        Args:
            excel_files: Paths to Excel files
            sheet_name: Sheet name containing annotations
            id_mapping: Mapping from annotation IDs to guide IDs
            isoform_selection: 'PickH' or list of isoforms, or None
            empty_sheet: Sheet name for guides without predicted mutation.
                         If None, tries 'editor' -> 'no_mutation' replacement,
                         then falls back to inferring from missing guides.
        
        Returns:
            DataFrame with parsed annotations
        """
        # Load main annotations
        vep = read_excel_sheets(excel_files, sheet_name)
        
        # Try to extract replicate column name from sheet name pattern
        # Expected format: "VEP editor Replicate" -> "Replicate"
        # Falls back to searching for a column containing 'Replicate' or 'replicate'
        
        if [c for c in vep.columns if 'replicate' in c.lower()] :
            replicate_cols = [c for c in vep.columns if 'replicate' in c.lower()]
            replicate_col = replicate_cols[0] if replicate_cols else 'Replicate'
        else :
            sheet_parts = sheet_name.split(' ')
            if len(sheet_parts) != 3 :
                raise ValueError(
                f"The annotations must have a Replicate column"
            )
            replicate_col = sheet_parts[2]

        
        # Filter to relevant guides
        vep = vep[vep['#Uploaded_variation'].isin(set(id_mapping.keys()))].copy()
        vep['ID'] = vep['#Uploaded_variation'].map(id_mapping)
        
        # Load empty/no-mutation guides
        empties = self._load_empty_guides(excel_files, sheet_name, id_mapping, vep, empty_sheet)
        
        # Create no-mutation DataFrame
        no_mutation_df = pd.DataFrame({
            'ID': empties,
            'Location': None,
            'Consequence': 'No predicted Mutation',
            'Protein_position': None,
            'IMPACT': None,
            'Amino_acids': None,
            replicate_col: None
        })
        
        # Apply isoform selection
        vep = self._apply_isoform_selection(vep, isoform_selection)
        
        # Combine with no-mutation guides
        vep = pd.concat([vep, no_mutation_df], ignore_index=True)
        
        # Parse annotations
        results = list(self._parse_annotations(vep, replicate_col, id_mapping))
        
        return pd.DataFrame(results, columns=[
            'ID', 'consequence', 'consequence_detail', 'impact',
            'POS_AA', 'mutations_AA', 'replicate'
        ])
    
    def _load_empty_guides(
        self,
        excel_files: Sequence[Union[str, Path]],
        sheet_name: str,
        id_mapping: dict[str, str],
        vep_df: pd.DataFrame,
        empty_sheet: Optional[str] = None
    ) -> list[str]:
        """
        Load guides with no predicted mutations.
        
        Args:
            excel_files: Paths to Excel files
            sheet_name: Main VEP annotation sheet name
            id_mapping: Mapping from annotation IDs to guide IDs
            vep_df: DataFrame with VEP annotations (already filtered)
            empty_sheet: Explicit sheet name for empty guides.
                         If None, tries 'editor' -> 'no_mutation' replacement.
        
        Returns:
            List of guide IDs with no predicted mutation
        """
        # Determine which sheet to try for empty guides
        if empty_sheet is not None:
            # Explicit sheet name provided
            target_sheet = empty_sheet
            fallback_reason = f'Specified empty sheet "{empty_sheet}" not found.'
        else:
            # Try replacing 'editor' with 'no_mutation'
            target_sheet = sheet_name.replace('editor', 'no_mutation')
            
            # Check if replacement actually changed anything
            if target_sheet == sheet_name:
                # 'editor' wasn't in the sheet name, fall back immediately
                print(
                    f'Warning: Sheet name "{sheet_name}" does not contain "editor". '
                    'Guides not in VEP annotations will be marked as no predicted mutation.'
                )
                return list(set(id_mapping.values()) - set(vep_df['ID']))
            
            fallback_reason = f'Empty guides sheet "{target_sheet}" not found.'
        
        # Try to load the empty guides sheet
        try:
            no_vep = read_excel_sheets(excel_files, target_sheet)
            no_vep = no_vep[no_vep['Guides'].isin(set(id_mapping.keys()))].copy()
            no_vep['ID'] = no_vep['Guides'].map(id_mapping)
            
            # Validate completeness (compare mapped IDs)
            all_mapped_ids = set(id_mapping.values())
            missing = all_mapped_ids - set(vep_df['ID']) - set(no_vep['ID'])
            if missing:
                raise ValueError(
                    f'Guides without any annotations: {missing}'
                )
            
            return no_vep['ID'].tolist()
            
        except Exception as e:
            error_str = str(e).lower()
            if 'missing' in error_str or 'not found' in error_str or 'no sheet' in error_str:
                print(
                    f'Warning: {fallback_reason} '
                    'Guides not in VEP annotations will be marked as no predicted mutation.'
                )
                return list(set(id_mapping.values()) - set(vep_df['ID']))
            raise
    
    def _apply_isoform_selection(
        self,
        vep: pd.DataFrame,
        selection: Optional[Union[str, list[str]]]
    ) -> pd.DataFrame:
        """Apply isoform selection to VEP DataFrame."""
        already_filtered = len(set(vep['ID'])) == len(vep['ID'])
        
        if already_filtered and selection:
            print('Warning: Annotations appear to already be filtered')
        
        if selection == 'PickH':
            return vep[vep['PICK'] == '1'].copy()
        elif isinstance(selection, list):
            isoform_filter = IsoformFilter(vep)
            return isoform_filter.filter_by_isoforms(selection)
        
        return vep
    
    def _parse_annotations(
        self,
        vep: pd.DataFrame,
        replicate_col: str,
        id_mapping: dict[str, str]
    ) -> Iterator[tuple]:
        """Parse individual annotations from VEP DataFrame."""
        for _, row in vep.iterrows():
            guide_id = row['ID']
            prot_pos = row.get('Protein_position', '')
            aa = row.get('Amino_acids', '')
            cons = row.get('Consequence', '')
            impact = row.get('IMPACT')
            replicate = row.get(replicate_col)
            
            # Handle NaN values
            prot_pos = '' if pd.isna(prot_pos) else str(prot_pos)
            aa = '' if pd.isna(aa) else str(aa)
            cons = '' if pd.isna(cons) else str(cons)
            
            # Parse replicates
            replicate_str = self._parse_replicates(replicate, id_mapping)
            
            # Parse amino acid position and mutations
            pos_aa, mutations_aa = self._parse_amino_acid_changes(prot_pos, aa)
            
            # Determine consequence category
            consequence, consequence_detail = self._classify_consequence(cons,mutations_aa)
            
            yield (guide_id, consequence, consequence_detail, impact,
                   pos_aa, mutations_aa, replicate_str)
    
    def _parse_replicates(
        self,
        replicate: Optional[str],
        id_mapping: dict[str, str]
    ) -> str:
        """Parse replicate string, mapping IDs."""
        if pd.isna(replicate) or replicate is None:
            return ""
        
        return ",".join(
            id_mapping[i.strip()]
            for i in str(replicate).split(',')
            if i.strip() in id_mapping
        )
    
    def _parse_amino_acid_changes(
        self,
        prot_pos: str,
        amino_acids: str
    ) -> tuple[Optional[int], str]:
        """Parse protein position and amino acid changes."""
        if '?' in prot_pos or '/' not in amino_acids:
            part = prot_pos.split('-')[0]
            part = part if part != '?' else prot_pos.split('-')[-1]
            return safe_int(part), ""
        
        if '/' not in amino_acids:
            return None, ""
        
        aa_parts = amino_acids.split('/')
        if len(aa_parts) != 2:
            return None, ""
        
        # Parse position range
        if '-' in prot_pos:
            pos_parts = prot_pos.split('-')
            pos_initial = int(pos_parts[0])
            pos_final = int(pos_parts[1])
        else:
            pos_initial = pos_final = int(prot_pos)
        
        # Find differing positions
        arr_a = np.frombuffer(aa_parts[0].encode(), dtype='S1')
        arr_b = np.frombuffer(aa_parts[1].encode(), dtype='S1')
        
        # Handle length differences
        min_len = min(len(arr_a), len(arr_b))
        diff_idx = np.where(arr_a[:min_len] != arr_b[:min_len])[0].tolist()
        
        if not diff_idx:
            return int(np.round((pos_initial + pos_final) / 2)), ""
        
        mutations = [
            f"{aa_parts[0][i]}_{i + pos_initial}"
            for i in diff_idx
        ]
        
        start, end = min(diff_idx), max(diff_idx)
        pos_aa = int(np.round(pos_initial + (start + end) / 2))
        
        return pos_aa, ",".join(mutations)
    
    def _classify_consequence(self, consequence_str: str, AA: str) -> tuple[Optional[str], str, str]:
        """Classify VEP consequence string into category."""
        if not consequence_str:
            return None, 'N/A'
        
        mapped = [
    ('missense' if AA else 'synonymous') if self.consequence_mapping.get(c) == 'stop_retained_variant' else self.consequence_mapping.get(c)
    for c in consequence_str.split(',')
]



        if 'splice' in mapped:
            return 'splice', 'splice'
        elif 'non-sense' in mapped:
            return 'non-sense', 'non-sense'
        elif 'missense' in mapped:
            return 'missense', 'missense'
        elif 'synonymous' in mapped:
            return 'synonymous', 'synonymous'
        elif 'regulatory' in mapped:
            return None, 'regulatory'
        elif 'non_coding' in mapped:
            return None, 'non-coding'
        elif 'No predicted Mutation' in mapped:
            return 'No predicted Mutation', 'No predicted Mutation'
        else:
            return None, 'N/A'


class ControlAssigner:
    """Assigns control labels to guides."""
    
    def __init__(self, config: ControlConfig):
        """
        Initialize with control configuration.
        
        Args:
            config: Control configuration
        """
        self.config = config
        self._universal_matcher = type('', (), {'__contains__': lambda s, i: True})()
    
    def assign_controls(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Assign control labels to guides in DataFrame.
        
        Args:
            df: DataFrame with 'proteins', 'Consequence' columns
        
        Returns:
            DataFrame with 'Controls' column added
        """
        df = df.copy()
        df['Controls'] = None
        
        # Empty controls
        if self.config.empty_as_negative:
            mask = df['Consequence'] == 'No predicted Mutation'
            df.loc[mask, 'Controls'] = 'negative_control'
        
        # Negative controls
        df = self._assign_control_type(
            df,
            genes=self.config.negative_genes,
            consequences=self.config.negative_consequences,
            control_label='negative_control',
            consequence_label='Negative Control Gene'
        )
        
        # Positive controls
        df = self._assign_control_type(
            df,
            genes=self.config.positive_genes,
            consequences=self.config.positive_consequences,
            control_label='positive_control',
            consequence_label='Positive Control Gene'
        )
        
        return df
    
    def _assign_control_type(
        self,
        df: pd.DataFrame,
        genes: Optional[set[str]],
        consequences: Optional[set[str]],
        control_label: str,
        consequence_label: str
    ) -> pd.DataFrame:
        """Assign a specific control type to matching guides."""
        if not genes and not consequences:
            return df
        
        # Use universal matcher if one criterion not specified
        gene_matcher = genes if genes else self._universal_matcher
        cons_matcher = consequences if consequences else self._universal_matcher
        # Find matching guides
        mask = df.apply(
            lambda row: (
                row['Consequence'] in cons_matcher and
                row['proteins'] in gene_matcher
            ),
            axis=1
        )
        
        if not mask.any():
            raise ValueError(f'{control_label} options did not yield any controls')
        if sum(mask) < 10 :
            print(f'Warning : fewer than 10 {consequence_label}s were detected this may cause an issue.')
        if set(df.loc[mask,'proteins']) != set(genes) :
            raise ValueError(f'{set(genes) - set(df.loc[mask,'proteins'])} gene(s) / region(s) did not yield any controls')

        

        
        # Remove non-matching guides within control genes (if both specified)
        if genes and consequences:
            exclude_mask = df.apply(
                lambda row: (
                    not (row['Consequence'] in cons_matcher) and
                    row['proteins'] in gene_matcher
                ),
                axis=1
            )
            if any(exclude_mask):
                print(f'Warning: {sum(exclude_mask)} presumed controls were removed since they are within control regions/genes but did not have the desired consequence')
            df = df[~exclude_mask]
        
        df.loc[mask, 'Controls'] = control_label
        df.loc[mask, 'Consequence'] = consequence_label
        
        return df


class BaseSummaryGenerator:
    """Main class for generating base summaries."""
    
    def __init__(
        self,
        key_file: Union[str, Path],
        annotation_files: Sequence[Union[str, Path]],
        annotation_sheet: str,
        library_sheet: str = 'Library',
        empty_sheet: Optional[str] = None,
        isoform_selection: Optional[Union[str, list[str]]] = None,
        control_config: Optional[ControlConfig] = None,
        force: bool = False
    ):
        """
        Initialize the summary generator.
        
        Args:
            key_file: Path to guide key file
            annotation_files: Paths to VEP annotation Excel files
            annotation_sheet: Sheet name for VEP annotations
            library_sheet: Sheet name for library information
            empty_sheet: Sheet name for guides without predicted mutation.
                         If None, tries 'editor' -> 'no_mutation' replacement,
                         then falls back to inferring from missing guides.
            isoform_selection: 'PickH', list of isoforms, or None
            control_config: Control configuration
            force: Bypass validation checks
        """
        self.key_handler = KeyFileHandler(key_file)
        self.annotation_files = annotation_files
        self.annotation_sheet = annotation_sheet
        self.library_sheet = library_sheet
        self.empty_sheet = empty_sheet
        self.isoform_selection = isoform_selection
        self.control_config = control_config or ControlConfig()
        self.force = force
        
        # Validate control config
        self.control_config.validate(force=force)
        
        # Load and prepare data
        self._prepare_data()
    
    def _prepare_data(self) -> None:
        """Load and prepare annotation data."""
        # Load library sheet
        library_df = read_excel_sheets(
            self.annotation_files,
            self.library_sheet
        )
        
        # Validate key file coverage
        missing = self.key_handler.sequences - set(library_df['protospacer'])
        if missing:
            if not self.force:
                missing_ids = [self.key_handler.get_id(s) for s in missing]
                raise ValueError(
                    f'Unannotated guides found: {missing_ids}'
                )
            print(f"Warning: Removing {len(missing)} unannotated guides (force mode)")
        
        # Filter to relevant guides
        library_df = library_df[
            library_df['protospacer'].isin(self.key_handler.sequences)
        ].copy()
        
        # Create ID mappings
        library_df['original_ID'] = library_df['ID']
        library_df['ID'] = library_df['protospacer'].map(self.key_handler.sequence_to_id)
        
        self.id_mapping = dict(zip(library_df['original_ID'], library_df['ID']))
        self.id_to_sequence = dict(zip(library_df['ID'], library_df['protospacer']))
        self.protein_mapping = dict(zip(library_df['ID'], library_df['Protein']))
        
        # Parse VEP annotations
        vep_parser = VEPAnnotationParser()
        self.annotations = vep_parser.parse(
            self.annotation_files,
            self.annotation_sheet,
            self.id_mapping,
            isoform_selection=self.isoform_selection,
            empty_sheet=self.empty_sheet
        )
        
        # Build annotation lookup dictionaries
        self._build_annotation_lookups()
    
    def _build_annotation_lookups(self) -> None:
        """Build dictionaries for quick annotation lookup."""
        self.consequence_map = dict(zip(
            self.annotations['ID'],
            self.annotations['consequence']
        ))
        self.consequence_detail_map = dict(zip(
            self.annotations['ID'],
            self.annotations['consequence_detail']
        ))
        self.impact_map = dict(zip(
            self.annotations['ID'],
            self.annotations['impact']
        ))
        self.position_map = dict(zip(
            self.annotations['ID'],
            self.annotations['POS_AA']
        ))
        self.mutations_map = dict(zip(
            self.annotations['ID'],
            self.annotations['mutations_AA']
        ))
        self.replicate_map = dict(zip(
            self.annotations['ID'],
            self.annotations['replicate']
        ))
    
    def process_mageck_file(
        self,
        mageck_file: Union[str, Path],
        treatment_levels: list[list[str]]
    ) -> tuple[str, pd.DataFrame]:
        """
        Process a single MaGeCK file.
        
        Args:
            mageck_file: Path to MaGeCK output file
            treatment_levels: Parsed treatment levels for label extraction
        
        Returns:
            Tuple of (sheet_label, processed_dataframe)
        """
        # Parse MaGeCK file
        data = MaGeCKParser.parse(mageck_file)
        
        # Validate guide coverage
        missing = set(self.annotations['ID']) - set(data['id'])
        if missing:
            print(f'Warning: {len(missing)} Guides were not found in the MageCK file : {" ".join(missing)} \n')
            data = data[~data['id'].isin(missing)]
        
        # Filter to annotated guides
        data = data[data['id'].isin(set(self.annotations['ID']))].copy()
        
        # Add sequence
        data.insert(1, 'sgRNA_seq', data['id'].map(self.id_to_sequence))
        
        # Extract protein from guide ID
        data['proteins'] = data['id'].map(self.protein_mapping)
        
        # Add annotations
        data['Consequence'] = data['id'].map(self.consequence_map)
        data['Consequence_Detail'] = data['id'].map(self.consequence_detail_map)
        data['Position'] = data['id'].map(self.position_map)
        data['Mutations'] = data['id'].map(self.mutations_map)
        data['replicates_window'] = data['id'].map(self.replicate_map)
        
        # Assign controls
        control_assigner = ControlAssigner(self.control_config)
        data = control_assigner.assign_controls(data)
        
        # Generate sheet label
        label = "_".join(extract_labels_from_filename(
            str(mageck_file),
            treatment_levels
        ))
        
        return label, data
    
    def generate_summary(
        self,
        mageck_files: Sequence[Union[str, Path]],
        treatment_spec: str,
        output_path: Union[str, Path]
    ) -> None:
        """
        Generate summary Excel file from MaGeCK files.
        
        Args:
            mageck_files: Paths to MaGeCK output files
            treatment_spec: Treatment specification string
            output_path: Path for output Excel file
        """
        treatment_levels = parse_treatment_levels(treatment_spec)
        output_path = Path(output_path)
        
        if not output_path.suffix:
            output_path = output_path.with_suffix('.xlsx')
        
        print('Generating summary...')
        
        with pd.ExcelWriter(output_path) as writer:
            for mageck_file in mageck_files:
                print(f'Processing: {mageck_file}')
                
                try:
                    label, data = self.process_mageck_file(
                        mageck_file,
                        treatment_levels
                    )
                    data.to_excel(writer, sheet_name=label, index=False)
                except Exception as e:
                    raise ValueError(
                        f'Error processing {mageck_file}: {e}'
                    ) from e
        
        print(f'Summary written to: {output_path}')


# =============================================================================
# CLI
# =============================================================================

def create_argument_parser() -> argparse.ArgumentParser:
    """Create and configure argument parser."""
    parser = argparse.ArgumentParser(
        description='Create base concatenated summary for MaGeCK analyses',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        '-k', '--key',
        metavar='FILE',
        required=True,
        dest='key_file',
        help='Tab/comma/space delimited key file with ID_guide and Sequence columns'
    )
    parser.add_argument(
        '-I', '--input',
        metavar='FILE',
        required=True,
        dest='mageck_files',
        nargs='+',
        help='MaGeCK output files'
    )
    parser.add_argument(
        '-a', '--annotation_excel',
        metavar='FILE',
        dest='annotation_files',
        nargs='+',
        required=True,
        help='Excel file(s) with VEP annotations'
    )
    parser.add_argument(
        '-s', '--sheet',
        metavar='STR',
        dest='annotation_sheet',
        required=True,
        help='Sheet name for VEP annotations'
    )
    parser.add_argument(
        '-X', '--xvar',
        dest='treatment_spec',
        required=True,
        help='Comma-separated treatment conditions (e.g., UNT/TREAT,KO/WT)'
    )
    
    # Optional arguments
    parser.add_argument(
        '-e', '--empty_sheet',
        metavar='STR',
        dest='empty_sheet',
        required=False,
        default=None,
        help='Sheet name for guides without predicted mutation. If not provided, '
             'tries replacing "editor" with "no_mutation" in annotation sheet name, '
             'otherwise infers from guides missing in VEP annotations.'
    )
    parser.add_argument(
        '-l', '--lib_sheet',
        metavar='STR',
        dest='library_sheet',
        default='Library',
        help='Sheet name for library information (default: Library)'
    )
    parser.add_argument(
        '--Isoform',
        metavar='STR',
        dest='isoforms',
        nargs='+',
        help='Isoforms to select'
    )
    parser.add_argument(
        '--Pick',
        dest='use_pick',
        action='store_true',
        help="Use VEP's PICK column for annotation selection"
    )
    parser.add_argument(
        '--Empty_controls',
        dest='empty_as_negative',
        action='store_true',
        help='Mark guides with no predicted mutation as negative controls'
    )
    parser.add_argument(
        '-F',
        dest='force',
        action='store_true',
        help='Bypass validation checks'
    )
    parser.add_argument(
        '--Negative_Control_Genes',
        metavar='STR',
        dest='negative_genes',
        nargs='+',
        help='Negative control genes/regions'
    )
    parser.add_argument(
        '--Positive_Control_Genes',
        metavar='STR',
        dest='positive_genes',
        nargs='+',
        help='Positive control genes/regions'
    )
    parser.add_argument(
        '--Negative_Control_Consequences',
        metavar='STR',
        dest='negative_consequences',
        nargs='+',
        choices=list(VALID_CONSEQUENCE_CATEGORIES),
        help='Negative control consequence types'
    )
    parser.add_argument(
        '--Positive_Control_Consequences',
        metavar='STR',
        dest='positive_consequences',
        nargs='+',
        choices=list(VALID_CONSEQUENCE_CATEGORIES),
        help='Positive control consequence types'
    )
    parser.add_argument(
        '--out',
        dest='output',
        default='summary',
        help='Output file name (default: summary)'
    )
    
    return parser


def main(args: Optional[list[str]] = None) -> int:
    """
    Main entry point.
    
    Args:
        args: Command line arguments (uses sys.argv if None)
    
    Returns:
        Exit code (0 for success)
    """
    parser = create_argument_parser()
    parsed = parser.parse_args(args)
    
    # Validate isoform selection
    if parsed.use_pick and parsed.isoforms:
        parser.error('Cannot use both --Pick and --Isoform')
    
    # Determine isoform selection mode
    if parsed.use_pick:
        isoform_selection = 'PickH'
    elif parsed.isoforms:
        isoform_selection = parsed.isoforms
    else:
        isoform_selection = None
    
    # Build control configuration
    control_config = ControlConfig(
        negative_genes=set(parsed.negative_genes) if parsed.negative_genes else None,
        positive_genes=set(parsed.positive_genes) if parsed.positive_genes else None,
        negative_consequences=set(parsed.negative_consequences) if parsed.negative_consequences else None,
        positive_consequences=set(parsed.positive_consequences) if parsed.positive_consequences else None,
        empty_as_negative=parsed.empty_as_negative
    )
    
    try:
        # Create generator
        generator = BaseSummaryGenerator(
            key_file=parsed.key_file,
            annotation_files=parsed.annotation_files,
            annotation_sheet=parsed.annotation_sheet,
            library_sheet=parsed.library_sheet,
            empty_sheet=parsed.empty_sheet,
            isoform_selection=isoform_selection,
            control_config=control_config,
            force=parsed.force
        )
        
        # Generate summary
        generator.generate_summary(
            mageck_files=parsed.mageck_files,
            treatment_spec=parsed.treatment_spec,
            output_path=parsed.output
        )
        
        return 0
        
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        raise


if __name__ == "__main__":
    sys.exit(main())