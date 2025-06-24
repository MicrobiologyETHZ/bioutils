

def setup_notebook(project_root="."):

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    import seaborn as sns
    import plotly.express as px
    import yaml
    import plotly.graph_objects as go

    sns.set_context("notebook", font_scale=1.1)
    pd.set_option("display.max_columns", 100)
    pd.set_option("display.max_rows", 100)
    plt.rcParams["figure.figsize"] = (12, 8)
    plt.rcParams['savefig.dpi'] = 200
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['axes.titlesize'] = 20
    plt.rcParams['font.size'] = 16
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 8
    plt.rcParams['legend.fontsize'] = 14
    # True activates latex output in fonts!
    plt.rcParams['text.usetex'] = False
    pd.set_option('display.float_format', lambda x: '{:,.2f}'.format(x))
    from datetime import date
    today = date.today().strftime("%d-%m-%y")

    # Load config
    config_path = Path(project_root) / "config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    return config, Path, pd, np, plt, sns, px, go, today
