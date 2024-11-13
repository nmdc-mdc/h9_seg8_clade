# -*- coding: UTF-8 -*-
'''Author:Kesheng Peng
e-mill:kesheng_peng_yx@163.com
Version:1.0
release-date:2023-10-19'''

import os
import argparse
import pandas as pd


def parse_args():
    """ Get command line parameters"""
    parser = argparse.ArgumentParser(description="SNP & annotation!")
    parser.add_argument("-w",
                        "--work_path",
                        required=True)
    args = parser.parse_args()

    return args


def merge_annotation_result(work_path):
    """ Merge annotation result to annotation_result_table.csv"""
    output_path = f"{work_path}/annotation_result_table.csv"
    if os.path.exists(output_path):
        os.remove(output_path)

    annotation_results = f"{work_path}/annotation_result"
    result = []
    for file_name in os.listdir(annotation_results):
        file_path = f"{annotation_results}/{file_name}"
        df = pd.read_csv(file_path, header=None,
                         dtype="str", keep_default_na=False,
                         low_memory=False)
        df.columns = ["segment_accession", "reference", "gene",
                      "ref_site", "ref_nucleicacid", "nucleicacid", "site",
                      "variant_type", "cds_name",
                      "ref_aminoacid_site", "ref_aminoacid", "aminoacid", "variant_amino_type"]
        result.append(df)

    df = pd.concat(result)
    df.to_csv(output_path, index=None, mode="a")


def merge_snp_result(work_path):
    """ Merge snp to snp_reuslt_table.csv """
    output_path = f"{work_path}/snp_reuslt_table.csv"
    if os.path.exists(output_path):
        os.remove(output_path)

    differ_site = f"{work_path}/differ_site"
    result = []
    for file_name in os.listdir(differ_site):
        file_path = f"{differ_site}/{file_name}"
        df = pd.read_csv(file_path,
                         dtype="str",
                         keep_default_na=False,
                         low_memory=False)
        mask = (df["ref"].str.contains("-") == False) & (df["query"].str.contains("-") == False)
        snp_df = df[mask]
        snp_df_differ = snp_df[snp_df["ref"] != snp_df["query"]]
        result.append(snp_df_differ)

    df = pd.concat(result)
    df.to_csv(output_path, mode="a", index=None)


if __name__ == "__main__":
    params = parse_args()
    work_path = params.work_path

    # Merge annotation result to annotation_result_table.csv
    merge_annotation_result(work_path)

    # Merge snp to snp_reuslt_table.csv
    merge_snp_result(work_path)
