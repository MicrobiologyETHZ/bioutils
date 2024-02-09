import pandas as pd
from pathlib import Path
from typing import Union
import pandas as pd
import numpy as np
from numpy.random import RandomState
import warnings
from typing import List
import plotly.express as px



def validate_sample_data(sample_data: pd.DataFrame):

    """
    Make sure sample data dataframe is compatible with downstream use
    
    """
    return 'sample_id' in sample_data.columns


def parse_fastq_screen_results(screen_dir: Union[str, Path], sample_data: pd.DataFrame):
    """
    Look at all the files with _screen.txt extension, merge with sample data
    Sample data format: has to contain sample_id column
    Bonus: look for relationships with metatdata vars
    
    """
    assert validate_sample_data(sample_data)
    files = list(screen_dir.rglob("*screen.txt"))
    if len(files) == 0:
        print("No Fastq_screen results found")
        return pd.DataFrame()
    dfs = []
    for file in files:
        df = pd.read_table(file, skiprows=1)
        not_mapped = df.iloc[-1,0]
        df = df[:-1]
        to_keep = ['Genome', "%One_hit_one_genome", "%Multiple_hits_one_genome", "%One_hit_multiple_genomes", "%Multiple_hits_multiple_genomes"]
        df = (df[to_keep].melt(id_vars=['Genome'],var_name='Hit', value_name='Percent')
              .assign(sample_id=file.stem.split("_screen")[0]))
        sample_id_expanded = df.sample_id.str.split("_", expand=True)
        df['orient'] = sample_id_expanded.iloc[:,-1]
        df['sample_id'] = sample_id_expanded[[0,1]].agg('_'.join, axis=1)
        dfs.append(df)
    return pd.concat(dfs).merge(sample_data, on='sample_id', how='outer')


def parse_star_log_out():
    pass

def parse_salmon_log_out():
    pass

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

    res = None
    noccur = np.sum(x)
    nvar = len(x)
    # Check for invalid count vector
    if noccur == 0 or nvar == 0:
        print("Invalid count vector x")
        return np.array([np.nan]*nvar)
    
    # Check if the number of iterations is within a reasonable range
    if iterations > 100000:
        warnings.warn(UserWarning('Max number of iterations allowed is 100000, setting to 100000'))
        iterations = 100000

    # Check if the rarefaction depth exceeds the total count in 'x'
    if depth > noccur:
        return np.array([np.nan]*nvar)
    p = x/noccur
    seeds = np.random.choice(100000, size=iterations) if iterations > 1 else [seed]
    
    # Perform rarefaction for each iteration
    for seed in seeds:
        prng = RandomState(seed)
        choice = prng.choice(nvar, size=depth, p=p)

        # Concatenate the results for each iteration
        if res is None:
            res = np.bincount(choice, minlength=nvar)[np.newaxis,:]
        else:
            res = np.concatenate((res, np.bincount(choice, minlength=nvar)[np.newaxis, :]))
    # Return the mean of rarefied counts over multiple iterations
    return np.nanmean(res, axis=0)


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
            rare_df = self.feature_df.apply(rarefy, depth=depth, iterations=self.iterations, seed=self.seed).assign(depth=depth)
            df_list.append(rare_df)
        # Concatenate the rarefied DataFrames and reset the index
        rdf = pd.concat(df_list).reset_index()
        self.rare_df = rdf.melt(id_vars=[rdf.columns[0], 'depth'], var_name='sample_id', value_name='counts')



    def saturation_curves(self):

        # Rarefy dataframe 
        self.rarefy_dataframe()
    
        # Calculate saturation curves based on cutoffs
        self.sat_curve_df = (self.rare_df.groupby(['sample_id', 'depth'])
                        .agg({'counts': [lambda x, c=c: (x > c).sum() for c in self.cutoffs]})
                        .reset_index())
        self.sat_curve_df.columns  = ['sample_id', 'depth'] + [f'>{c} reads' for c in self.cutoffs]
        self.sat_curve_df = self.sat_curve_df.melt(id_vars=['sample_id', 'depth'], var_name=['cutoffs'], value_name = ['num_features'])


    def plot_saturation_curve(self, sample_id: str):
        if self.sat_curve_df.empty:
            self.saturation_curves()
        fig = px.scatter(self.sat_curve_df[self.sat_curve_df.sample_id == sample_id].sort_values('depth'),
                         x = 'depth', y = 'num_features', color='cutoffs', height=600, width=800)



