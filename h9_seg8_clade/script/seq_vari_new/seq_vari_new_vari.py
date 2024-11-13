import json
import copy
import os.path
import shutil

import pandas as pd
import argparse


def get_input_vari(input_vari_path):
    # 读取变异表
    input_vari = pd.read_csv(input_vari_path,
                             dtype="str", keep_default_na=False)
    input_vari = input_vari[input_vari["variant_amino_type"] == "nonsynonymous SNV"]
    input_vari = input_vari.drop_duplicates(subset=[
        "cds_name", "ref_aminoacid_site", "ref_aminoacid", "aminoacid"
    ])
    input_vari_new_vari = {}
    groups = input_vari.groupby(["segment_accession", 'reference', 'cds_name'])
    for name, sub_df in groups:
        segment_accession = name[0]
        if segment_accession not in input_vari_new_vari:
            input_vari_new_vari[segment_accession] = {}
            cds_name = name[2]
        if cds_name not in input_vari_new_vari[segment_accession]:
            input_vari_new_vari[segment_accession][cds_name] = {
                "ref": name[1], "variants": []
            }
        amino_site = []
        amino_vari = []
        for row in sub_df.values.tolist():
            ref_aminoacid_site = row[9]
            ref_aminoacid = row[10]
            aminoacid = row[11]
            amino_site.append(int(ref_aminoacid_site))
            amino_vari.append([ref_aminoacid, ref_aminoacid_site, aminoacid])
        amino_site = sorted(amino_site)
        sort_amino_vari = [""] * len(amino_site)
        for vari in amino_vari:
            ref_aminoacid_site = int(vari[1])
            index = amino_site.index(ref_aminoacid_site)
            sort_amino_vari[index] = "".join(vari)

        input_vari_new_vari[segment_accession][cds_name]["variants"] = sort_amino_vari

    data = []
    for k, v in input_vari_new_vari.items():
        cds_name = sorted(list(v.keys()))
        vari = [""] * len(cds_name)
        for k1, v1 in v.items():
            ref = v1["ref"]
            variants = v1["variants"]
            index = cds_name.index(k1)
            vari[index] = k1 + "(" + " ".join(variants) + ")"

        row = [k, " ".join(vari)]
        data.append(row)

    return data


def get_new_vari(input_vari, node_parents_vari, node_vari, qc_filter_pass):
    input_vari_new_vari = []
    for row in input_vari:
        seq_id = row[0]
        seq_tree_node = qc_filter_pass[seq_id]["seq_tree_ref_name"]
        tree_node_vari = ""
        if seq_tree_node in node_vari:
            tree_node_vari = node_vari[seq_tree_node]
        tree_node_parents_vari = []
        if seq_tree_node in node_parents_vari:
            tree_node_parents_vari = node_parents_vari[seq_tree_node]
        if tree_node_vari != "":
            tree_node_cds = tree_node_vari.split("; ")
            for cds in tree_node_cds:
                cds_ls = cds.split(": ")
                cds_name = cds_ls[0]
                cds_vari = cds_ls[1]
                for vari in cds_vari.split(", "):
                    if cds_vari not in tree_node_parents_vari[cds_name]:
                        tree_node_parents_vari[cds_name].append(cds_vari)

        seq_vari = row[1]
        seq_cds = seq_vari.split("; ")
        seq_new_vari_dict = {}
        for cds in seq_cds:
            cds_ls = cds.split("(")
            cds_name = cds_ls[0]
            if cds_name not in seq_new_vari_dict:
                seq_new_vari_dict[cds_name] = []
            for vari in cds_ls[1].strip(")").split(" "):
                if vari not in tree_node_parents_vari[cds_name]:
                    seq_new_vari_dict[cds_name].append(vari)

        seq_new_vari = []
        cds_list = sorted(list(seq_new_vari_dict.keys()))
        for cds in cds_list:
            vari = seq_new_vari_dict[cds_name]
            if vari != []:
                seq_new_vari.append(cds + "(" + " ".join(vari) + ")")

        if seq_new_vari != []:
            seq_new_vari = "; ".join(seq_new_vari)
        else:
            seq_new_vari = ""

        text = row + [seq_new_vari]
        input_vari_new_vari.append(text)

    return input_vari_new_vari


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--input_isolate_name', required=True, help='Input isolate name')
    parser.add_argument('-w', '--work_dir', required=True, help='work_dir')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    params = parse_args()
    isolate_name = params.input_isolate_name
    work_dir = params.work_dir
    input_vari_path = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{isolate_name}/annotation_result_table.csv"
    input_vari = get_input_vari(input_vari_path)

    with open(f"{work_dir}/h9_seg8_clade/data/node_parents_vari.json") as file:
        node_parents_vari = json.load(file)

    node_varis_df = pd.read_csv(f"{work_dir}/h9_seg8_clade/data/seq_variations_table.csv", keep_default_na=False)
    node_vari = {x[0]: x[1] for x in node_varis_df.values.tolist()}

    with open(f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}/qc_filter_pass.json") as file:
        qc_filter_pass = json.load(file)

    input_vari_new_vari = get_new_vari(
        input_vari, node_parents_vari, node_vari, qc_filter_pass
    )

    df = pd.DataFrame(input_vari_new_vari, columns=["seq_id", "vari", "new_vari"])
    
    os.makedirs(f"{work_dir}/h9_seg8_clade/script/seq_vari_new/work_dir/{isolate_name}")
    df.to_csv(f"{work_dir}/h9_seg8_clade/script/seq_vari_new/work_dir/{isolate_name}/input_seq_vari_new_vari.csv", index=False)
