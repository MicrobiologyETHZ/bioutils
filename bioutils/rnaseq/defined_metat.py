"""
Metatranscriptomic Data Processing Package
For analyzing defined bacterial communities with known composition
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pyranges as pr
import plotly.express as px
from sklearn.decomposition import PCA
from typing import Dict, List, Optional, Union
import warnings
import json
import yaml


class StrainMapper:
    """
    Handles mapping between chromosome/contig names and strain identifiers
    """

    def __init__(self, strain_mapping: Optional[Dict[str, str]] = None,
                 mapping_file: Optional[Union[str, Path]] = None,
                 config_file: Optional[Union[str, Path]] = None):
        """
        Initialize strain mapper

        Args:
            strain_mapping: Dictionary mapping contig/chromosome -> strain name
            mapping_file: Path to CSV file with columns ['contig', 'strain']
            config_file: Path to config file (JSON/YAML) with strain mapping
        """
        if strain_mapping:
            self.strain_map = strain_mapping
        elif mapping_file:
            self.strain_map = self._load_from_csv(mapping_file)
        elif config_file:
            self.strain_map = self._load_from_config(config_file)
        else:
            # Try to load from default config file
            self.strain_map = self._load_default_config()

    def _load_from_csv(self, mapping_file: Union[str, Path]) -> Dict[str, str]:
        """Load strain mapping from CSV file"""
        df = pd.read_csv(mapping_file)
        if 'contig' not in df.columns or 'strain' not in df.columns:
            raise ValueError(
                "CSV file must contain 'contig' and 'strain' columns")
        return dict(zip(df['contig'], df['strain']))

    def _load_from_config(self, config_file: Union[str, Path]) -> Dict[str, str]:
        """Load strain mapping from config file (JSON or YAML)"""
        config_path = Path(config_file)

        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        try:
            with open(config_path, 'r') as f:
                if config_path.suffix.lower() in ['.yaml', '.yml']:
                    config = yaml.safe_load(f)
                elif config_path.suffix.lower() == '.json':
                    config = json.load(f)
                else:
                    raise ValueError("Config file must be JSON or YAML format")

            # Extract strain mapping from config
            if 'strain_mapping' in config:
                return config['strain_mapping']
            elif 'strains' in config:
                return config['strains']
            else:
                # Assume entire config is the strain mapping
                return config

        except Exception as e:
            raise ValueError(f"Error loading config file {config_path}: {e}")

    def _load_default_config(self) -> Dict[str, str]:
        """Load strain mapping from default config file"""
        # Look for config files in common locations
        module_dir = Path(__file__).parent
        default_locations = [
            module_dir / 'strain_config.yaml',
            module_dir / 'strain_config.yml',
            module_dir / 'strain_config.json',
            module_dir / 'config' / 'strain_mapping.yaml',
            module_dir / 'config' / 'strain_mapping.yml',
            module_dir / 'config' / 'strain_mapping.json',
            Path.cwd() / 'strain_config.yaml',
            Path.cwd() / 'strain_config.yml',
            Path.cwd() / 'strain_config.json',
            Path.cwd() / 'config' / 'strain_mapping.yaml',
            Path.cwd() / 'config' / 'strain_mapping.yml',
            Path.cwd() / 'config' / 'strain_mapping.json',
            Path.home() / '.config' / 'metatranscriptomics' / 'strain_mapping.yaml'
        ]

        for config_path in default_locations:
            if config_path.exists():
                print(f"Loading strain mapping from: {config_path}")
                return self._load_from_config(config_path)

        # If no config file found, create a template and use fallback
        self._create_config_template()
        return self._get_fallback_mapping()

    def _create_config_template(self):
        """Create a template config file for user reference"""
        template_path = Path.cwd() / 'strain_config_template.yaml'

        if not template_path.exists():
            template_config = {
                'strain_mapping': self._get_fallback_mapping(),
                '_info': {
                    'description': 'Strain mapping configuration',
                    'format': 'contig_id: strain_name',
                    'usage': 'Rename this file to strain_config.yaml to use as default'
                }
            }

            with open(template_path, 'w') as f:
                yaml.dump(template_config, f,
                          default_flow_style=False, sort_keys=False)

            print(f"Created config template at: {template_path}")
            print("Rename to 'strain_config.yaml' and modify as needed")

    def _get_fallback_mapping(self) -> Dict[str, str]:
        """Fallback strain mapping if no config file is found"""
        print("Warning: No config file found, using fallback strain mapping")
        return {
            'CP015399.2': 'YL32', 'CP015400.2': 'KB18', 'CP015401.2': 'I48',
            'CP015402.2': 'YL27', 'CP015403.2': 'YL45', 'CP015404.2': 'I46',
            'CP015405.2': 'YL58', 'CP015406.2': 'YL31', 'CP015407.2': 'YL2',
            'CP015408.2': 'I49', 'CP015409.2': 'YL44', 'CP015410.2': 'KB1',
            'GCF_000364265': 'ASF519', 'FQ312003.1': 'SL1344',
            'HE654725.1': 'SL1344', 'HE654726.1': 'SL1344', 'HE654724.1': 'SL1344',
            'CP097573.1': 'ASF500', 'NZ_CP097810.1': 'ASF356',
            'NZ_CP097561.1': 'ASF361', 'NZ_CP097562.1': 'ASF457',
        }

    def save_config(self, config_file: Union[str, Path], format: str = 'yaml'):
        """
        Save current strain mapping to config file

        Args:
            config_file: Path to save config file
            format: File format ('yaml' or 'json')
        """
        config_path = Path(config_file)

        config_data = {
            'strain_mapping': self.strain_map,
            '_metadata': {
                'description': 'Strain mapping configuration',
                'total_strains': len(set(self.strain_map.values())),
                'total_contigs': len(self.strain_map)
            }
        }

        config_path.parent.mkdir(parents=True, exist_ok=True)

        with open(config_path, 'w') as f:
            if format.lower() == 'yaml':
                yaml.dump(config_data, f, default_flow_style=False,
                          sort_keys=False)
            elif format.lower() == 'json':
                json.dump(config_data, f, indent=2)
            else:
                raise ValueError("Format must be 'yaml' or 'json'")

        print(f"Strain mapping saved to: {config_path}")

    def get_strain(self, contig: str) -> str:
        """Get strain name for a given contig"""
        return self.strain_map.get(contig, f'Unknown_{contig}')

    def add_mapping(self, contig: str, strain: str):
        """Add new contig -> strain mapping"""
        self.strain_map[contig] = strain

    def get_strains(self, include_unknown: bool = True) -> List[str]:
        strains = sorted(list(set(self.strain_map.values())))
        if not include_unknown:
            strains = [s for s in strains if not s.startswith('Unknown_')]
        return strains

    def get_contigs_for_strain(self, strain: str) -> List[str]:
        """Get list of contigs for a specific strain"""
        return [contig for contig, s in self.strain_map.items() if s == strain]


class GenomeAnnotation:
    """
    Handles GFF file processing and gene annotation
    """

    def __init__(self, gff_file: Union[str, Path], feature_type: str = "CDS"):
        """
        Initialize genome annotation

        Args:
            gff_file: Path to GFF3 file
            feature_type: Type of features to extract (gene, CDS, etc.)
        """
        self.gff_file = Path(gff_file)
        self.feature_type = feature_type
        self.annotation = self._process_gff()

    def _process_gff(self) -> pd.DataFrame:
        """Process GFF file and extract relevant features"""
        try:
            # Read GFF3 file
            gff = pr.read_gff3(self.gff_file).as_df()

            # Filter for desired feature type
            features = gff[gff['Feature'] == self.feature_type].copy()

            # Extract essential columns
            annotation_cols = ['Chromosome',
                               'Feature', 'Start', 'End', 'Strand']

            # Add ID column (try different attribute names)
            if 'ID' in features.columns:
                annotation_cols.append('ID')
            elif 'gene_id' in features.columns:
                features['ID'] = features['gene_id']
                annotation_cols.append('ID')
            elif 'locus_tag' in features.columns:
                features['ID'] = features['locus_tag']
                annotation_cols.append('ID')

            # Add other useful columns if they exist
            optional_cols = ['Name', 'locus_tag', 'gene_biotype', 'product']
            for col in optional_cols:
                if col in features.columns:
                    annotation_cols.append(col)

            return features[annotation_cols].reset_index(drop=True)

        except Exception as e:
            raise ValueError(f"Error processing GFF file {self.gff_file}: {e}")

    def get_gene_info(self, gene_id: str) -> pd.Series:
        """Get annotation info for a specific gene"""
        gene_info = self.annotation[self.annotation['ID'] == gene_id]
        if gene_info.empty:
            return pd.Series()
        return gene_info.iloc[0]


class MetatranscriptomicProcessor:
    """
    Main class for processing metatranscriptomic featureCounts data
    """

    def __init__(self, strain_mapper: StrainMapper,
                 genome_annotation: Optional[GenomeAnnotation] = None):
        """
        Initialize processor

        Args:
            strain_mapper: StrainMapper instance
            genome_annotation: GenomeAnnotation instance (optional)
        """
        self.strain_mapper = strain_mapper
        self.genome_annotation = genome_annotation
        self.count_data = pd.DataFrame()
        self.summary_stats = pd.DataFrame()

    def merge_featurecounts_files(self, count_dir: Union[str, Path],
                                  pattern: str = "*.count.txt") -> pd.DataFrame:
        """
        Merge multiple featureCounts output files into one table

        Args:
            count_dir: Directory containing featureCounts output files
            pattern: File pattern to match

        Returns:
            Merged count table
        """
        count_dir = Path(count_dir)
        files = list(count_dir.rglob(pattern))

        if not files:
            raise ValueError(
                f"No files found matching pattern {pattern} in {count_dir}")

        print(f"Found {len(files)} featureCounts files")

        # Process first file to establish the index structure
        first_file = files[0]
        sample_name = first_file.stem.split(".count")[0]
        print(f"Processing {sample_name}")

        df = pd.read_csv(first_file, sep='\t', comment='#')
        expected_cols = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
        if len(df.columns) == 7:
            df.columns = expected_cols + [sample_name]
        else:
            raise ValueError(f"Unexpected number of columns in {first_file}")

        # Set multi-column index for efficient joining
        merged_df = df[['Geneid', 'Chr', 'Length', sample_name]
                       ].set_index(['Geneid', 'Chr', 'Length'])

        # Process remaining files
        for file_path in files[1:]:
            sample_name = file_path.stem.split(".count")[0]
            print(f"Processing {sample_name}")

            df = pd.read_csv(file_path, sep='\t', comment='#')
            if len(df.columns) == 7:
                df.columns = expected_cols + [sample_name]
            else:
                raise ValueError(
                    f"Unexpected number of columns in {file_path}")

            # Set same index and join (much faster than merge)
            df_indexed = df[['Geneid', 'Chr', 'Length', sample_name]].set_index(
                ['Geneid', 'Chr', 'Length'])
            merged_df = merged_df.join(df_indexed, how='outer')

        # Reset index to get columns back
        merged_df = merged_df.reset_index()

        # Fill NaN values with 0
        sample_cols = [col for col in merged_df.columns if col not in [
            'Geneid', 'Chr', 'Length']]
        merged_df[sample_cols] = merged_df[sample_cols].fillna(0)

        # Add strain information
        merged_df['strain'] = merged_df['Chr'].apply(
            self.strain_mapper.get_strain)

        # Add gene annotation if available
        if self.genome_annotation:
            annotation_df = self.genome_annotation.annotation[[
                'ID', 'Chromosome']].copy()
            annotation_df = annotation_df.rename(
                columns={'ID': 'Geneid', 'Chromosome': 'Chr'})
            merged_df = merged_df.merge(
                annotation_df, on=['Geneid', 'Chr'], how='left')

        self.count_data = merged_df
        print(f"Merged data shape: {merged_df.shape}")
        return merged_df

    def parse_featurecounts_summary(self, count_dir: Union[str, Path],
                                    pattern="*.summary") -> pd.DataFrame:
        """
        Parse featureCounts summary files to get mapping statistics

        Args:
            count_dir: Directory containing .summary files

        Returns:
            Summary statistics DataFrame
        """
        count_dir = Path(count_dir)
        summary_files = list(count_dir.rglob(pattern))

        if not summary_files:
            warnings.warn("No .summary files found")
            return pd.DataFrame()

        summary_dfs = []
        for file_path in summary_files:
            # Extract sample name
            sample_name = file_path.stem.split(".count")[0]

            # Read summary file
            df = pd.read_csv(file_path, sep='\t')
            df['sample_id'] = sample_name

            # Rename columns for clarity
            df.columns = ['status', 'read_counts', 'sample_id']
            summary_dfs.append(df)

        # Combine all summaries
        combined_summary = pd.concat(summary_dfs, ignore_index=True)

        # Pivot to get stats per sample
        summary_pivot = combined_summary.pivot(
            index='sample_id', columns='status', values='read_counts')
        summary_pivot = summary_pivot.reset_index().fillna(0)

        # Calculate percentages
        summary_pivot['total_reads'] = summary_pivot.select_dtypes(
            include=[np.number]).sum(axis=1)

        for col in summary_pivot.columns:
            if col not in ['sample_id', 'total_reads']:
                summary_pivot[f'pct_{col.lower()}'] = (
                    summary_pivot[col] / summary_pivot['total_reads'] * 100).round(2)

        self.summary_stats = summary_pivot
        return summary_pivot

    def get_strain_summary(self) -> pd.DataFrame:
        """
        Generate summary statistics per strain

        Returns:
            DataFrame with strain-level statistics
        """
        if self.count_data.empty:
            raise ValueError(
                "No count data available. Run merge_featurecounts_files first.")

        # Get sample columns (exclude metadata columns)
        sample_cols = [col for col in self.count_data.columns
                       if col not in ['Geneid', 'Chr', 'Length', 'strain']]

        strain_stats = []

        for strain in self.count_data['strain'].unique():

            strain_data = self.count_data[self.count_data['strain'] == strain]

            strain_info = {
                'strain': strain,
                'num_genes': len(strain_data),
                'total_length': strain_data['Length'].sum()
            }

            # Calculate total reads per sample for this strain
            for sample in sample_cols:
                strain_info[f'{sample}_reads'] = strain_data[sample].sum()
                strain_info[f'{sample}_genes_with_reads'] = (
                    strain_data[sample] > 0).sum()

            strain_stats.append(strain_info)

        return pd.DataFrame(strain_stats)

    def calculate_tpm(self, pseudocount: float = 0.5) -> pd.DataFrame:
        """
        Calculate TPM (Transcripts Per Million) values

        Args:
            pseudocount: Value to add before log transformation

        Returns:
            DataFrame with TPM values
        """
        if self.count_data.empty:
            raise ValueError("No count data available")

        sample_cols = [col for col in self.count_data.columns
                       if col not in ['Geneid', 'Chr', 'Length', 'strain']]

        tpm_data = self.count_data.copy()

        for sample in sample_cols:
            # Calculate TPM: (reads / length_kb) / sum(reads / length_kb) * 1e6
            rpk = tpm_data[sample] / (tpm_data['Length'] / 1000)
            tpm_data[f'{sample}_tpm'] = rpk / rpk.sum() * 1e6
            tpm_data[f'{sample}_log2tpm'] = np.log2(
                tpm_data[f'{sample}_tpm'] + pseudocount)

        return tpm_data

    def get_composition_data(self, normalize: bool = True) -> pd.DataFrame:
        """
        Get strain composition data for plotting

        Args:
            normalize: Whether to normalize to proportions

        Returns:
            DataFrame suitable for composition plotting
        """
        if self.count_data.empty:
            raise ValueError("No count data available")

        sample_cols = [col for col in self.count_data.columns
                       if col not in ['Geneid', 'Chr', 'Length', 'strain']]

        # Sum reads per strain per sample
        composition = self.count_data.groupby('strain')[sample_cols].sum().T
        composition.index.name = 'sample_id'
        composition = composition.reset_index()

        # Melt for plotting
        comp_melted = composition.melt(
            id_vars='sample_id', var_name='strain', value_name='reads')

        if normalize:
            # Calculate proportions
            total_reads = comp_melted.groupby(
                'sample_id')['reads'].transform('sum')
            comp_melted['proportion'] = comp_melted['reads'] / total_reads

        return comp_melted
