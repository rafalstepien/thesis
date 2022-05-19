from vcf_parser import FileReader
from typing import List, Iterable


def filter_by_p_value(path_to_genotypes: str, path_to_gwas_dataset: str) -> str:
    gwas_data = FileReader.load_gwas_dataset(path_to_gwas_dataset)
    gwas_data_significant_ids = gwas_data[gwas_data.P < 0.05]
    gwas_data_significant_ids_without_X_chromosome = gwas_data_significant_ids[~ gwas_data_significant_ids.SNP.str.startswith("chrX")].SNP.values
    genotypes_file = FileReader.load_variants_table(path_to_genotypes)
    return filter_genotypes_file(genotypes_file, gwas_data_significant_ids_without_X_chromosome)


def filter_genotypes_file(genotypes_generator: Iterable, snp_ids: List) -> str:
    data = ""
    for _id, row in enumerate(genotypes_generator):
        if not row.startswith("SNP"):
            chromosome_id, _, _ = row.partition(';')
            if _id % 10_000:
                print(_id)
            if chromosome_id in snp_ids:
                data += row
        else:
            data += row
    return data


if __name__ == "__main__":
    PATH_TO_GENOTYPES_FILE = "/home/rafalstepien/Magisterka/code/main_analysis/data/output_table_filtered_fixed.csv"
    PATH_TO_GWAS_FILE = "/home/rafalstepien/Magisterka/code/genome_visualiser/data/gwas_results.txt"

    data = filter_by_p_value(PATH_TO_GENOTYPES_FILE, PATH_TO_GWAS_FILE)
    pass
