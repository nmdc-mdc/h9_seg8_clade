import os
import json
import shutil

import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
import argparse


def generate_table(qc_pass, qc_not_pass, input_seq_vari_new, seq_name_dict, meta_df):
    table = []
    for file_name, seq_info in seq_name_dict.items():
        name = seq_info["name"]
        oringin_name = seq_info["oringin_name"]
        qc_and_filter = ""
        if name in qc_pass:
            qc_and_filter = "PASS"
        if name in qc_not_pass:
            qc_and_filter = "FAIL"
        seq_tree_ref_name = "-"
        if name in qc_pass:
            seq_tree_ref_name = qc_pass[name]["seq_tree_ref_name"]
        segment = "-"
        clade_name = "-"
        if seq_tree_ref_name != "-":
            meta = meta_df[meta_df["seq_name"] == seq_tree_ref_name].values.tolist()[0]
            segment = meta[11]
            clade_name = meta[1]
        identity = "-"
        if name in qc_pass:
            identity = qc_pass[name]["seq_identity"]
        vari = ""
        new = ""
        for row in input_seq_vari_new.values.tolist():
            if row[0] == name:
                vari = row[1]
                new = row[2]
        text = [oringin_name, name, qc_and_filter, segment, clade_name,
                identity, vari, new, seq_tree_ref_name]
        table.append(text)

    return table


def input_seq_add_to_tree(tree_path, output_df, output_dir):
    seg_tree = {}
    output_dir = f"{output_dir}/tree"
    os.makedirs(output_dir, exist_ok=True)
    groups = output_df.groupby("segment")
    for seg, sub_df in groups:
        tree_file_path = tree_path[seg]
        tree = Phylo.read(tree_file_path, "newick")
        for row in sub_df.values.tolist():
            seq_in_tree_name = row[1]
            seq_tree_ref_name = row[8]
            node = tree.find_any(seq_tree_ref_name)
            # 创建新节点
            new_node = Clade(name=seq_in_tree_name)
            # 将新节点添加为原节点的子节点
            node.clades.append(new_node)
        output_path = f"{output_dir}/{seg}.nwk"
        Phylo.write(tree, output_path, 'newick')
        with open(output_path.replace(".nwk", ".txt"), "w") as f:
            f.write(str(tree))
        seg_tree[seg] = output_path

    return seg_tree


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--input_isolate_name', required=True, help='Input isolate name')
    parser.add_argument('-w', '--work_dir', required=True, help='work_dir')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    params = parse_args()
    isolate_name = params.input_isolate_name
    work_dir = params.work_dir
    output_dir = f"{work_dir}/h9_seg8_clade/script/output/{isolate_name}"
    os.makedirs(output_dir)
    with open(f"{work_dir}/h9_seg8_clade/script/parse_input/parse_output/{isolate_name}/seq_name_dict.json") as file:
        seq_name_dict = json.load(file)
    with open(f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}/qc_filter_pass.json") as file:
        qc_pass = json.load(file)
    with open(f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}/qc_filter_not_pass.json") as file:
        qc_not_pass = json.load(file)
    input_seq_vari_new = pd.read_csv(
        f"{work_dir}/h9_seg8_clade/script/seq_vari_new/work_dir/{isolate_name}/input_seq_vari_new_vari.csv",
        keep_default_na=False)
    meta_df = pd.read_csv(f"{work_dir}/h9_seg8_clade/data/all_meta_data_with_vari.csv", keep_default_na=False)
    meta_title = meta_df.columns.values.tolist()
    meta_dict = {}
    meta_groups = meta_df.groupby("seg_name")
    for seg_name, sub_df in meta_groups:
        meta_dict[seg_name] = [sub_df]

    output_table = generate_table(
        qc_pass, qc_not_pass, input_seq_vari_new, seq_name_dict, meta_df
    )
    title = ["inPutSeqName", "nameInTheTree", "qcAndFiler", "segment", "cladeName",
             "maximumIdenty", "variations", "variationsNotInAncestorsNodes",
             "seq_tree_ref_name"]
    output_df = pd.DataFrame(output_table, columns=title)
    output_df = output_df[output_df["segment"] != "-"]
    output_group = output_df.groupby("segment")
    for seg, sub_df in output_group:
        data = []
        for row in sub_df.values.tolist():
            meta_row = [""] * len(meta_title)
            meta_row[0] = row[1]
            meta_row[1] = row[4]
            meta_row[11] = row[3]
            meta_row[12] = row[6]
            meta_row[13] = row[7]
            data.append(meta_row)
        output_meta_df = pd.DataFrame(data, columns=meta_title)
        meta_dict[seg].append(output_meta_df)

    merge_meta_dict = {}
    for seq, meta_list in meta_dict.items():
        merge_meta_dict[seq] = pd.concat(meta_list, ignore_index=True)

    feature_json = {}
    for f in os.listdir(f"{work_dir}/h9_seg8_clade/data/feature_json"):
        f_seg_name = f.split(".")[0]
        f_path = f"{work_dir}/h9_seg8_clade/data/feature_json/{f}"
        feature_json[f_seg_name] = f_path

    tree_path = {x.split(".")[0]: f"{work_dir}/h9_seg8_clade/data/tree/{x}" for x in
                 os.listdir(f"{work_dir}/h9_seg8_clade/data/tree")}
    seg_tree = input_seq_add_to_tree(tree_path, output_df, output_dir)

    meta_out = f"{output_dir}/meta"
    os.makedirs(meta_out, exist_ok=True)
    for seq, df in merge_meta_dict.items():
        meta_out_path = f"{meta_out}/{seq}.csv"
        df.to_csv(meta_out_path, index=False)

    cladeRresultJson = {}
    for row in output_df.values.tolist():
        out_dir = f"{output_dir}/feature"
        os.makedirs(out_dir, exist_ok=True)
        feat_path = os.path.basename(feature_json[row[3]])
        featrue_path = f"{out_dir}/{feat_path}"

        source = feature_json[row[3]]
        destination = featrue_path

        # Check if source is a file or a directory, then copy accordingly
        if os.path.isfile(source):
            shutil.copy(source, destination)  # Copy file to destination
        elif os.path.isdir(source):
            shutil.copytree(source, destination, dirs_exist_ok=True)  # Copy directory recursively
        row_meta = merge_meta_dict[row[3]]
        meta_out_path = f"{meta_out}/{row[3]}.csv"

        row_dict = {
            "inPutSeqName": row[0],
            "nameInTheTree": row[1],
            "nameInTheTree": row[1],
            "qcAndFiler": row[2],
            "segment": row[3],
            "cladeName": row[4],
            "maximumIdenty": row[5],
            "variations": row[6],
            "variationsNotInAncestorsNodes": row[7],
            "featureVariations": featrue_path.replace("\\", "/"),
            "tree": seg_tree[row[3]].replace("\\", "/"),
            "meta": meta_out_path.replace("\\", "/"),
        }
        cladeRresultJson[row[1]] = row_dict

    with open(f"{output_dir}/cladeRresultJson.json", "w") as f:
        json.dump(cladeRresultJson, f, indent=4)
