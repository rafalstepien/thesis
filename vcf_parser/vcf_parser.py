import gzip
from constatnts import (
    VCF_FIRST_SAMPLE_ID_INDEX,
    SNP_ID_INDEX,
    INFO_INDEX,
    VCF_GENOTYPE_TO_GENOTYPE,
)
from typing import List, Iterable, Tuple, Union
import pandas as pd


class FileReader:
    """Class responsible for reading huge VCF files and returning a generator"""

    @staticmethod
    def load_vcf(filepath: str) -> Iterable:
        with gzip.open(filepath, "rb") as vcf:
            for row in vcf:
                yield row

    @staticmethod
    def load_gwas_dataset(filepath: str) -> pd.DataFrame:
        return pd.read_csv(filepath, delim_whitespace=True, usecols=["SNP", "P"])

    @staticmethod
    def load_variants_table(filepath: str) -> pd.DataFrame:
        with open(filepath, "r") as variants_table:
            for row in variants_table:
                yield row


class VCFHelper:

    @staticmethod
    def is_header_row(row: str) -> bool:
        return row.startswith("#") and not row.startswith("##")

    @staticmethod
    def is_data_row(row: str) -> bool:
        return not row.startswith("#")

    @staticmethod
    def append_row_to_file(data: str, output_path: str):
        if data:
            with open(output_path, "a") as f:
                f.write(data)


class VCFConverter(VCFHelper):
    def __init__(self, filepath: str, ouptut_path: str):
        self.filepath: str = filepath
        self.ouptut_path: str = ouptut_path
        self._genotype_position_in_info = None

    def convert_to_table(self):
        for row in FileReader.load_vcf(self.filepath):
            rows_as_string = self.generate_csv_row(row)
            self.append_row_to_file(rows_as_string)

    def generate_csv_row(self, row: bytes) -> Union[str, None]:
        decoded_row = row.decode("utf-8").strip()
        row_as_list = decoded_row.split("\t")

        if self.is_header_row(decoded_row):
            return self._generate_csv_header(row_as_list)

        elif self.is_data_row(decoded_row):
            return self._generate_csv_row(row_as_list)

    @staticmethod
    def _generate_csv_header(row: List) -> str:
        return "SNP_ID;" + ";".join(row[VCF_FIRST_SAMPLE_ID_INDEX:])

    def _generate_csv_row(self, row: List) -> str:
        if not self._genotype_position_in_info:
            self._genotype_position_in_info = row[INFO_INDEX].split(":").index("GT")
        return self._get_snp_id(row) + ";" + self._get_genotypes(row) + "\n"

    @staticmethod
    def _get_snp_id(row: List) -> str:
        return row[SNP_ID_INDEX]

    def _get_genotypes(self, row: List) -> str:
        genotypes_of_all_patients = []

        for genotype_of_patient in row[VCF_FIRST_SAMPLE_ID_INDEX:]:
            genotype_of_patient = self._extract_genotype_from_string(
                genotype_of_patient
            )
            genotype_encoded = VCF_GENOTYPE_TO_GENOTYPE.get(genotype_of_patient).value
            genotypes_of_all_patients.append(genotype_encoded)

        return ";".join(genotypes_of_all_patients)

    def _extract_genotype_from_string(self, genotype_string: str) -> str:
        return genotype_string.split(":")[self._genotype_position_in_info]


class VCFFilter:
    POSITION_OF_VARIANT_ID = 0
    NAME_OF_COLUMN_WITH_SNP_IDS = "SNP"

    def __init__(self, filepath: str):
        self.filepath = filepath

    def filter_by_list_of_variants(self, list_of_variants: pd.DataFrame, output_filepath: str):

        for row in FileReader.load_variants_table(self.filepath):
            if self._is_header(row):
                VCFHelper.append_row_to_file(row, output_filepath)

            else:
                snp_id = row.split(";")[VCFFilter.POSITION_OF_VARIANT_ID]
                if snp_id in list_of_variants[VCFFilter.NAME_OF_COLUMN_WITH_SNP_IDS].values:
                    VCFHelper.append_row_to_file(row, output_filepath)

    @staticmethod
    def _is_header(row):
        return row.startswith("SNP_ID")


if __name__ == "__main__":

    INPUT_TEST_FILE_PATH = "/home/rafalstepien/Magisterka/subset_multisample.vcf.gz"
    OUTPUT_TABLE_PATH_TEST = "/home/rafalstepien/Magisterka/code/vcf_parser/data/test_output_table.csv"

    INPUT_FILE_PATH = "/home/rafalstepien/Magisterka/multisample_20210716_just_snps.vcf.gz"
    OUTPUT_TABLE_PATH = "/home/rafalstepien/Magisterka/code/vcf_parser/data/output_table.csv"

    GWAS_INPUT_PATH = "/home/rafalstepien/Magisterka/code/genome_visualiser/data/gwas_results.txt"
    OUTPUT_FILTERED_PATH = "/home/rafalstepien/Magisterka/code/vcf_parser/data/output_table_filtered.csv"

    dataset = FileReader.load_gwas_dataset(GWAS_INPUT_PATH)
    VCFFilter(OUTPUT_TABLE_PATH).filter_by_list_of_variants(dataset, OUTPUT_FILTERED_PATH)
