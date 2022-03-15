import pandas as pd
import json
import dash_bio as dashbio
from dash import html
import numpy as np
from dash import Dash

from constants import SNP_CONFIG, LAYOUT_CONFIG


class GWASGenomeVisualizer(Dash):

    def __init__(self, gwas_file_path: str, chromosomes_file_path: str):
        self.gwas_file_path = gwas_file_path
        self.gwas_dataframe = None
        self.gwas_dataset_json = None
        self.chromosomes_file_path = chromosomes_file_path
        self.chromosomes_dataset = None
        self._load_gwas_dataset()
        self._load_chromosomes_data()
        self._parse_dataset()
        super().__init__(__name__)

    def run(self):
        self._add_circos_object()
        self.run_server(debug=True)

    def _load_gwas_dataset(self):
        self.gwas_dataframe = pd.read_csv(self.gwas_file_path, delim_whitespace=True, usecols=["SNP", "P"])

    def _load_chromosomes_data(self):
        with open(self.chromosomes_file_path, 'r') as file:
            self.chromosomes_dataset = json.load(file)

    def _parse_dataset(self):
        dataset = self.gwas_dataframe[self.gwas_dataframe["P"] < 0.05]["SNP"].str.split("_", expand=True).iloc[:, :2]
        dataset.columns = ["chromosome", "position"]
        dataset["value"] = np.random.uniform(0.001, 0.1, dataset.shape[0])
        dataset["block_id"] = dataset["chromosome"]
        dataset["position"] = dataset["position"].astype(int)
        self.gwas_dataset_json = [data for key, data in dataset.to_dict(orient="index").items()]

    def _add_circos_object(self):
        self.layout = html.Div([dashbio.Circos(
            layout=self.chromosomes_dataset["GRCh38"],
            config=LAYOUT_CONFIG,
            tracks=[
                {
                    "id": "snp",
                    "type": "SCATTER",
                    "data": self.gwas_dataset_json,
                    "config": SNP_CONFIG,
                },
            ],
            size=1920,
            enableDownloadSVG=True,
            enableZoomPan=True,
        )])


if __name__ == "__main__":
    app = GWASGenomeVisualizer("data/gwas_results.txt", "data/chromosome_lengths_GRCh38.json")
    app.run()
