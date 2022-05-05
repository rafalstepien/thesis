import os
import logging
import gzip
from constatnts import (
    VCF_FIRST_SAMPLE_ID_INDEX,
    SNP_ID_INDEX,
    INFO_INDEX,
    VCF_GENOTYPE_TO_GENOTYPE,
)
from typing import List, Iterable, Union


class VCFParser:
    def __init__(self, filepath: str, ouptut_path: str):
        self.filepath: str = filepath
        self.ouptut_path: str = ouptut_path
        self._genotype_position_in_info = None

    def convert_to_table(self):
        vcf_reader_generator = self.read_vcf_rows()
        for _id, row in enumerate(vcf_reader_generator):
            if _id % 10_000 == 0:
                logging.info(f"Saved line number {_id}")

            rows_as_string = self.generate_csv_row(row)

            self.append_row_to_file(rows_as_string)

    def read_vcf_rows(self) -> Iterable:
        with gzip.open(self.filepath, "rb") as vcf:
            for row in vcf:
                yield row

    def generate_csv_row(self, row: bytes) -> Union[str, None]:
        decoded_row = row.decode("utf-8").strip()
        row_as_list = decoded_row.split("\t")

        if self._is_header_row(decoded_row):
            return self._generate_csv_header(row_as_list)

        elif self._is_data_row(decoded_row):
            return self._generate_csv_row(row_as_list)

    @staticmethod
    def _generate_csv_header(row: List) -> str:
        return "SNP_ID;" + ";".join(row[VCF_FIRST_SAMPLE_ID_INDEX:])

    def _generate_csv_row(self, row: List) -> str:
        if not self._genotype_position_in_info:
            self._genotype_position_in_info = row[INFO_INDEX].split(":").index("GT")
        return self._get_snp_id(row) + ";" + self._get_genotypes(row) + "\n"

    def append_row_to_file(self, data: str):
        if data:
            with open(self.ouptut_path, "a") as f:
                f.write(data)

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

    @staticmethod
    def _is_header_row(row: str) -> bool:
        return row.startswith("#") and not row.startswith("##")

    @staticmethod
    def _is_data_row(row: str) -> bool:
        return not row.startswith("#")


if __name__ == "__main__":

    INPUT_TEST_FILE_PATH = "/home/rafalstepien/Magisterka/subset_multisample.vcf.gz"
    OUTPUT_TABLE_PATH_TEST = "/home/rafalstepien/Magisterka/code/vcf_parser/data/test_output_table.csv"

    INPUT_FILE_PATH = "/home/rafalstepien/Magisterka/multisample_20210716_just_snps.vcf.gz"
    OUTPUT_TABLE_PATH = "/home/rafalstepien/Magisterka/code/vcf_parser/data/output_table.csv"

    VCFParser(INPUT_FILE_PATH, OUTPUT_TABLE_PATH).convert_to_table()
