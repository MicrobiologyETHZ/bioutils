"""
Metatranscriptomic Data Analysis and Visualization Module

Functions for common analyses on normalized count data including:
- PCA analysis
- Hierarchical clustering  
- Correlation analysis
- Visualization utilities
"""

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
from typing import Dict, List, Optional, Union, Tuple
import warnings


class MetatranscriptomicAnalyzer:
    """
    Analysis and visualization tools for metatranscriptomic data
    """
    
    def __init__(self, processor=None):
        """
        Initialize analyzer
        
        Args:
            processor: MetatranscriptomicProcessor instance (optional)
        """
        self.processor = processor
        
    def prepare_data_for_analysis(self, data_df: pd.DataFrame, 
                                 value_type: str = 'clr',
                                 metadata: Optional[pd.DataFrame] = None,
                                 sample_id_col: str = 'sample_id') -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Prepare normalized data for analysis by extracting sample matrix and metadata
        
        Args:
            data_df: DataFrame with normalized data (from calculate_clr/calculate_tpm)
            value_type: Type of values to extract ('clr', 'tpm', 'log2tpm', 'raw')
            metadata: Sample metadata DataFrame (optional)
            sample_id_col: Column name in metadata containing sample IDs
            
        Returns:
            Tuple of (sample_matrix, sample_metadata)
            - sample_matrix: genes x samples matrix for analysis
            - sample_metadata: metadata for samples
        """
        # Identify sample columns based on value_type
        if value_type == 'clr':
            sample_cols = [col for col in data_df.columns if col.endswith('_clr')]
            sample_names = [col[:-4] for col in sample_cols]  # Remove '_clr' suffix
        elif value_type == 'tpm':
            sample_cols = [col for col in data_df.columns if col.endswith('_tpm') and not col.endswith('_log2tpm')]
            sample_names = [col[:-4] for col in sample_cols]  # Remove '_tpm' suffix
        elif value_type == 'log2tpm':
            sample_cols = [col for col in data_df.columns if col.endswith('_log2tpm')]
            sample_names = [col[:-8] for col in sample_cols]  # Remove '_log2tpm' suffix
        elif value_type == 'raw':
            # Identify raw sample columns (not metadata, not transformed)
            metadata_cols = ['Geneid', 'Chr', 'Length', 'strain']
            transformed_suffixes = ['_clr', '_tpm', '_log2tpm']
            sample_cols = [col for col in data_df.columns 
                          if col not in metadata_cols and 
                          not any(col.endswith(suffix) for suffix in transformed_suffixes)]
            sample_names = sample_cols
        else:
            raise ValueError(f"Unknown value_type: {value_type}. Use 'clr', 'tpm', 'log2tpm', or 'raw'")
            
        if not sample_cols:
            raise ValueError(f"No {value_type} columns found in data")
            
        # Create sample matrix (genes x samples)
        sample_matrix = data_df[sample_cols].copy()
        sample_matrix.columns = sample_names  # Clean column names
        sample_matrix.index = data_df.get('Geneid', range(len(data_df)))
        
        # Prepare sample metadata
        if metadata is not None:
            # Filter metadata to samples present in data
            sample_metadata = metadata[metadata[sample_id_col].isin(sample_names)].copy()
            # Reorder to match sample matrix column order
            sample_metadata = sample_metadata.set_index(sample_id_col).reindex(sample_names).reset_index()
        elif self.processor and hasattr(self.processor, 'metadata') and self.processor.metadata is not None:
            # Use metadata from processor
            metadata_df = self.processor.metadata
            sample_metadata = metadata_df[metadata_df[sample_id_col].isin(sample_names)].copy()
            sample_metadata = sample_metadata.set_index(sample_id_col).reindex(sample_names).reset_index()
        else:
            # Create minimal metadata
            sample_metadata = pd.DataFrame({sample_id_col: sample_names})
            
        print(f"Prepared {value_type} data: {sample_matrix.shape[0]} genes x {sample_matrix.shape[1]} samples")
        
        return sample_matrix, sample_metadata
        
    def calculate_sample_correlations(self, sample_matrix: pd.DataFrame, 
                                    method: str = 'spearman') -> pd.DataFrame:
        """
        Calculate pairwise correlations between samples
        
        Args:
            sample_matrix: genes x samples matrix
            method: Correlation method ('spearman', 'pearson')
            
        Returns:
            samples x samples correlation matrix
        """
        if method == 'spearman':
            # Calculate Spearman correlation
            corr_matrix = sample_matrix.corr(method='spearman')
        elif method == 'pearson':
            corr_matrix = sample_matrix.corr(method='pearson')
        else:
            raise ValueError(f"Unknown method: {method}. Use 'spearman' or 'pearson'")
            
        return corr_matrix
        
    def perform_pca(self, sample_matrix: pd.DataFrame, 
                   sample_metadata: pd.DataFrame,
                   n_components: int = 3,
                   standardize: bool = True) -> Dict:
        """
        Perform PCA analysis on sample data
        
        Args:
            sample_matrix: genes x samples matrix
            sample_metadata: sample metadata DataFrame
            n_components: Number of principal components to calculate
            standardize: Whether to standardize features before PCA
            
        Returns:
            Dictionary with PCA results
        """
        # Transpose for PCA (samples x genes)
        X = sample_matrix.T
        
        # Standardize if requested
        if standardize:
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
        else:
            X_scaled = X.values
            
        # Perform PCA
        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X_scaled)
        
        # Create results DataFrame
        pc_columns = [f'PC{i+1}' for i in range(n_components)]
        pca_df = pd.DataFrame(X_pca, columns=pc_columns, index=sample_matrix.columns)
        
        # Merge with metadata
        pca_results = pca_df.reset_index().rename(columns={'index': 'sample_id'})
        if 'sample_id' in sample_metadata.columns:
            pca_results = pca_results.merge(sample_metadata, on='sample_id', how='left')
        
        # Calculate explained variance
        explained_variance = pca.explained_variance_ratio_
        
        results = {
            'pca_data': pca_results,
            'explained_variance': explained_variance,
            'loadings': pd.DataFrame(
                pca.components_.T, 
                columns=pc_columns,
                index=sample_matrix.index
            ),
            'pca_object': pca
        }
        
        print(f"PCA completed. Explained variance: {explained_variance[:3]}")
        
        return results
        
    def hierarchical_clustering(self, sample_matrix: pd.DataFrame,
                              method: str = 'average',
                              metric: str = 'correlation') -> Dict:
        """
        Perform hierarchical clustering on samples
        
        Args:
            sample_matrix: genes x samples matrix  
            method: Linkage method ('average', 'single', 'complete', 'ward')
            metric: Distance metric ('correlation', 'euclidean', 'manhattan')
            
        Returns:
            Dictionary with clustering results
        """
        # Transpose for clustering (samples x genes)
        X = sample_matrix.T
        
        # Calculate distance matrix
        if metric == 'correlation':
            # Use 1 - correlation as distance
            corr_matrix = np.corrcoef(X)
            distance_matrix = 1 - corr_matrix
            # Convert to condensed form for linkage
            distances = squareform(distance_matrix, checks=False)
        else:
            distances = pdist(X, metric=metric)
            
        # Perform hierarchical clustering
        linkage_matrix = linkage(distances, method=method)
        
        results = {
            'linkage_matrix': linkage_matrix,
            'sample_names': sample_matrix.columns.tolist(),
            'distance_matrix': squareform(distances) if metric != 'correlation' else distance_matrix
        }
        
        return results
        
    def plot_pca(self, pca_results: Dict, 
                 color_by: Optional[str] = None,
                 size_by: Optional[str] = None,
                 title: str = "PCA Analysis") -> go.Figure:
        """
        Create interactive PCA plot
        
        Args:
            pca_results: Results from perform_pca()
            color_by: Column name to color points by
            size_by: Column name to size points by  
            title: Plot title
            
        Returns:
            Plotly figure
        """
        pca_data = pca_results['pca_data']
        explained_var = pca_results['explained_variance']
        
        # Create 2D plot
        fig = px.scatter(
            pca_data,
            x='PC1', y='PC2',
            color=color_by,
            size=size_by,
            hover_data=pca_data.columns.tolist(),
            title=title,
            labels={
                'PC1': f'PC1 ({explained_var[0]:.1%} variance)',
                'PC2': f'PC2 ({explained_var[1]:.1%} variance)'
            }
        )
        
        fig.update_layout(
            width=800, height=600,
            showlegend=True
        )
        
        return fig
        
    def plot_correlation_heatmap(self, sample_matrix: pd.DataFrame,
                               sample_metadata: Optional[pd.DataFrame] = None,
                               method: str = 'spearman',
                               cluster: bool = True,
                               figsize: Tuple[int, int] = (10, 8)) -> plt.Figure:
        """
        Plot correlation heatmap with optional clustering
        
        Args:
            sample_matrix: genes x samples matrix
            sample_metadata: sample metadata for annotation
            method: Correlation method ('spearman', 'pearson')
            cluster: Whether to cluster samples
            figsize: Figure size
            
        Returns:
            Matplotlib figure
        """
        # Calculate correlations
        corr_matrix = self.calculate_sample_correlations(sample_matrix, method)
        
        # Set up the plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create heatmap
        if cluster:
            sns.clustermap(
                corr_matrix, 
                method='average',
                metric='correlation',
                cmap='RdBu_r',
                center=0,
                figsize=figsize,
                cbar_kws={'label': f'{method.title()} Correlation'}
            )
            plt.tight_layout()
            return plt.gcf()
        else:
            sns.heatmap(
                corr_matrix,
                annot=False,
                cmap='RdBu_r',
                center=0,
                square=True,
                ax=ax,
                cbar_kws={'label': f'{method.title()} Correlation'}
            )
            
        ax.set_title(f'Sample Correlation Heatmap ({method.title()})')
        plt.tight_layout()
        
        return fig
        
    def plot_dendrogram(self, clustering_results: Dict,
                       sample_metadata: Optional[pd.DataFrame] = None,
                       color_threshold: Optional[float] = None,
                       figsize: Tuple[int, int] = (12, 6)) -> plt.Figure:
        """
        Plot dendrogram from hierarchical clustering
        
        Args:
            clustering_results: Results from hierarchical_clustering()
            sample_metadata: Sample metadata for annotation
            color_threshold: Threshold for coloring clusters
            figsize: Figure size
            
        Returns:
            Matplotlib figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        dendrogram(
            clustering_results['linkage_matrix'],
            labels=clustering_results['sample_names'],
            ax=ax,
            color_threshold=color_threshold,
            leaf_rotation=45,
            leaf_font_size=10
        )
        
        ax.set_title('Hierarchical Clustering Dendrogram')
        ax.set_xlabel('Samples')
        ax.set_ylabel('Distance')
        
        plt.tight_layout()
        return fig