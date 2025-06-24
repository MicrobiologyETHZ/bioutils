import pandas as pd
from pathlib import Path
from typing import Union, List, Dict, Tuple, Optional
import numpy as np
from numpy.random import RandomState
import warnings
import plotly.express as px
import plotly.graph_objects as go

# Saturation curves


def rarefy(x, depth=1000, iterations=1, seed=42):
    """
    Rarefies a count or frequency vector 'x' by randomly subsampling elements.

    Parameters:
    - x (numpy.ndarray): Input count or frequency vector to be rarefied. Meant to represent gene or species counts.
    - depth (int, optional): The desired rarefaction depth, i.e., the number of elements to subsample.
                             Default is 1000.
    - iterations (int, optional): The number of iterations to perform rarefaction.
                                  Default is 1. If > 1, random list of seeds is generated. Overules seed param.
    - seed (int, optional): Seed for reproducibility of random sampling. Default is 42.

    Returns:
    numpy.ndarray: Rarefied vector with the same length as the input vector 'x'.
                  The result is the mean of rarefied counts over multiple iterations.
                  If the number of iterations exceeds 100000, a warning is printed, and the
                  number of iterations is set to 100000.
                  If the rarefaction depth exceeds the total count in 'x', an array of NaNs with the
                  same length as 'x' is returned.
                  If 'x' has zero counts or length zero, an array of NaNs with the same length as 'x' is returned.
    """

    # Convert pandas Series or lists to numpy array if needed
    if not isinstance(x, np.ndarray):
        x = np.array(x)

    noccur = np.sum(x)
    nvar = len(x)

    # Check for invalid count vectors
    if noccur == 0:
        warnings.warn("Input vector has zero counts")
        return np.array([np.nan] * nvar)

    if nvar == 0:
        warnings.warn("Input vector has zero length")
        return np.array([])

    # Check if the number of iterations is within a reasonable range
    if iterations <= 0:
        raise ValueError("Number of iterations must be positive")

    if iterations > 100000:
        warnings.warn(
            'Max number of iterations allowed is 100000, setting to 100000')
        iterations = 100000

    # Check if the rarefaction depth exceeds the total count in 'x'
    if depth <= 0:
        raise ValueError("Rarefaction depth must be positive")

    if depth > noccur:
        warnings.warn(
            f"Rarefaction depth ({depth}) exceeds total count in vector ({noccur})")
        return np.array([np.nan]*nvar)

     # Calculate probability vector
    p = x/noccur
    seeds = np.random.choice(
        100000, size=iterations) if iterations > 1 else [seed]

    # Initialize results array
    results = np.zeros((iterations, nvar))

    # Perform rarefaction for each iteration
    for i, s in enumerate(seeds):
        prng = RandomState(s)
        choice = prng.choice(nvar, size=depth, p=p)
        results[i, :] = np.bincount(choice, minlength=nvar)
    # Return the mean of rarefied counts over multiple iterations
    return np.nanmean(results, axis=0)  # or np.mean?


class SaturationCurve:

    def __init__(self, df: pd.DataFrame, depths: List[int], cutoffs: List[int] = [5, 10], iterations: int = 1, seed: int = 42) -> None:

        self.feature_df = df
        # Set column name for sample IDs
        self.feature_df.columns.name = 'sample_id'
        self.depths = depths
        self.cutoffs = cutoffs
        self.iterations = iterations
        self.seed = seed
        self.rare_df = pd.DataFrame()
        self.sat_curve_df = pd.DataFrame()

    def rarefy_dataframe(self):
        """
        Perform rarefaction on counts or frequencies in a DataFrame across multiple rarefaction depths

        Returns:
        pd.DataFrame: A melted DataFrame containing the rarefied counts or frequencies for each sample
                    at different rarefaction depths. The resulting DataFrame has the following columns:
                    - The original column name representing the sample ID (before melting).
                    - 'depth': The rarefaction depth applied to the counts or frequencies.
                    - 'sample_id': The melted column name representing the sample ID.
                    - 'counts': The rarefied counts or frequencies corresponding to the 'sample_id' at the given 'depth'.
        """
        df_list = []
        # Iterate over specified rarefaction depths
        for depth in self.depths:
            rare_df = self.feature_df.apply(
                rarefy, depth=depth, iterations=self.iterations, seed=self.seed).assign(depth=depth)
            df_list.append(rare_df)
        # Concatenate the rarefied DataFrames and reset the index
        rdf = pd.concat(df_list).reset_index()
        self.rare_df = rdf.melt(
            id_vars=[rdf.columns[0], 'depth'], var_name='sample_id', value_name='counts')

    def saturation_curves(self):

        # Rarefy dataframe
        self.rarefy_dataframe()

        # Calculate saturation curves based on cutoffs
        self.sat_curve_df = (self.rare_df.groupby(['sample_id', 'depth'])
                             .agg({'counts': [lambda x, c=c: (x > c).sum() for c in self.cutoffs]})
                             .reset_index())
        self.sat_curve_df.columns = [
            'sample_id', 'depth'] + [f'>{c} reads' for c in self.cutoffs]
        self.sat_curve_df = self.sat_curve_df.melt(id_vars=['sample_id', 'depth'], var_name=[
                                                   'cutoffs'], value_name=['num_features'])

    def plot_saturation_curve(self, sample_id: str):
        if self.sat_curve_df.empty:
            self.saturation_curves()
        fig = px.scatter(self.sat_curve_df[self.sat_curve_df.sample_id == sample_id].sort_values('depth'),
                         x='depth', y='num_features', color='cutoffs', height=600, width=800)
