import pandas as pd
import json
import dash_bio as dashbio
from dash import html
import numpy as np
from dash import Dash

from constants import SNP_CONFIG, LAYOUT_CONFIG


class GWASGenomeVisualizer(Dash):

    def __init__(self, snp_file_path: str, chromosomes_file_path: str):
        self.snp_file_path = snp_file_path
        self.snp_dataframe = None
        self.snp_dataset_json = None
        self.chromosomes_file_path = chromosomes_file_path
        self.chromosomes_dataset = None
        self._load_snp_data()
        self._load_chromosomes_data()
        self._parse_dataset()
        super().__init__(__name__)

    def run(self):
        self._add_circos_object()
        self.run_server(debug=True)

    def _load_snp_data(self):
        self.snp_dataframe = pd.read_csv(self.snp_file_path)

    def _load_chromosomes_data(self):
        with open(self.chromosomes_file_path, 'r') as file:
            self.chromosomes_dataset = json.load(file)

    def _parse_dataset(self):
        dataset = self.snp_dataframe.iloc[:, 0].str.split("_", expand=True).iloc[:, :2]
        dataset.columns = ["chromosome", "position"]
        dataset["value"] = np.random.uniform(0.1, 1, dataset.shape[0])
        dataset["block_id"] = dataset["chromosome"]
        dataset["position"] = dataset["position"].astype(int)
        self.snp_dataset_json = [data for key, data in dataset.to_dict(orient="index").items()]

    def _add_circos_object(self):
        self.layout = html.Div([dashbio.Circos(
            layout=self.chromosomes_dataset["GRCh38"],
            config=LAYOUT_CONFIG,
            tracks=[
                {
                    "id": "snp",
                    "type": "SCATTER",
                    "data": self.snp_dataset_json,
                    "config": SNP_CONFIG,
                },
            ],
            size=1920,
            enableDownloadSVG=True,
            enableZoomPan=True,
        )])


if __name__ == "__main__":
    app = GWASGenomeVisualizer(
        "data/data_after_removal_of_highly_correlated.csv",
        "data/chromosome_lengths_GRCh38.json"
    )
    # app = GWASGenomeVisualizer(
    #     "data/significant_features_no_chrX_without_seleciton.csv",
    #     "data/chromosome_lengths_GRCh38.json"
    # )
    app.run()
