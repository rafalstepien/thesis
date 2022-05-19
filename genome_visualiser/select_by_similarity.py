import pandas as pd


def _load_gwas_dataset(path):
    return pd.read_csv(path, delim_whitespace=True, usecols=["SNP", "P"])


dataset = _load_gwas_dataset("data/gwas_results.txt")
