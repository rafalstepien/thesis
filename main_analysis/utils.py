import pandas as pd
from typing import Tuple
import matplotlib.pyplot as plt
import seaborn as sns


def get_percentage_of_missing_values(dataframe: pd.DataFrame) -> pd.Series:
    print(f"In the dataset there is {dataframe.shape[0]} rows")
    na_values = pd.DataFrame(dataframe.isna().sum(axis = 0), columns=['number_of_na'])
    return na_values / dataframe.shape[0] * 100


def encode_dataframe_to_plot_heatmap_of_missing_values(dataframe: pd.DataFrame) -> pd.DataFrame:
    missing_data_visualization = dataframe.replace(1, 0)
    missing_data_visualization = missing_data_visualization.replace(2, 0)
    missing_data_visualization = missing_data_visualization.replace(-1, 1)
    return missing_data_visualization.fillna(value=1)


def plot_and_save_heatmap(dataframe: pd.DataFrame, figsize: Tuple, colors: Tuple, filename: str) -> None:
    plt.figure(figsize = figsize)
    
    heatmap = sns.heatmap(dataframe, cbar=False, xticklabels=False, yticklabels=False, cmap=sns.color_palette(colors))
    heatmap.set_ylabel("Patient ID", fontsize = 20)
    heatmap.set_xlabel("SNP ID", fontsize = 20)
    
    figure = heatmap.get_figure()    
    figure.savefig(f'plots/{filename}.png', dpi=100)
    
    
    
def remove_by_variance_threshold(dataframe: pd.DataFrame, threshold: float) -> pd.DataFrame:
    variance_table = pd.DataFrame(dataframe.var(), columns=["variance"])
    ids_to_remove = variance_table[variance_table.variance < threshold].index.values
    return dataframe.drop(ids_to_remove, axis=1)
    
