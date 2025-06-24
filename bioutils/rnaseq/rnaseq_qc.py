import pandas as pd
from pathlib import Path
from typing import Union, List, Dict, Tuple, Optional
import numpy as np
from numpy.random import RandomState
import warnings
import plotly.express as px
import plotly.graph_objects as go


def validate_sample_data(sample_data: pd.DataFrame):
    """
    Make sure sample data dataframe is compatible with downstream use

    """
    if 'sample_id' not in sample_data.columns:
        warnings.warn("Required column 'sample_id' not found in sample data")
        return False
    return True


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
        not_mapped = df.iloc[-1, 0]
        df = df[:-1]
        to_keep = ['Genome', "%One_hit_one_genome", "%Multiple_hits_one_genome",
                   "%One_hit_multiple_genomes", "%Multiple_hits_multiple_genomes"]
        df = (df[to_keep].melt(id_vars=['Genome'], var_name='Hit', value_name='Percent')
              .assign(sample_id=file.stem.split("_screen")[0]))
        sample_id_expanded = df.sample_id.str.split("_", expand=True)
        df['orient'] = sample_id_expanded.iloc[:, -1]
        df['sample_id'] = sample_id_expanded[[0, 1]].agg('_'.join, axis=1)
        dfs.append(df)
    return pd.concat(dfs).merge(sample_data, on='sample_id', how='outer')


def parse_star_log_out():
    pass


def parse_salmon_log_out():
    pass
