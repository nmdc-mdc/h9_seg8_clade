import os
import json
import shutil
import argparse
import pandas as pd
import glob


def copy_files_recursively(source_pattern, destination_dir):
    # Ensure the destination directory exists
    os.makedirs(destination_dir, exist_ok=True)

    # Use glob to find all matching files in the source pattern
    files = glob.glob(source_pattern)

    for file_path in files:
        # If it's a file, copy it to the destination directory
        if os.path.isfile(file_path):
            shutil.copy(file_path, destination_dir)
        # If it's a directory, copy the entire directory recursively
        elif os.path.isdir(file_path):
            de_path = os.path.basename(file_path)
            destination_subdir = f"{destination_dir}/{de_path}"
            shutil.copytree(file_path, destination_subdir, dirs_exist_ok=True)


def get_input_seq(parse_output):
    files = os.listdir(parse_output)
    input_seqs = {}
    for f in files:
        f_path = f"{parse_output}/{f}"
        with open(f_path, "r") as f_obj:
            seq = f_obj.read().strip().split("\n")
        seq_name = seq[0].split(">")[1]
        seq_seq = "".join(seq[1:])
        seq_length = len(seq_seq)
        input_seqs[seq_name] = {
            "seq": seq[1:],
            "length": seq_length,
            "file_name": f.split(".")[0]
        }
    return input_seqs


def get_input_seq_ref_seq_filter(input_seqs, meta_df, ref_seq, blast_file):
    # 获取每条序列对上的最近序列名
    blast_result_files = os.listdir(blast_file)
    seq_tree_ref = {}
    for f in blast_result_files:
        df = pd.read_csv(
            f"{blast_file}/{f}",
            sep="\t"
        )
        first_row = df.iloc[0].tolist()
        identity = first_row[2]
        query_name = first_row[0]
        tree_ref_name = first_row[1]
        q_start = int(first_row[6])
        q_end = int(first_row[7])
        blastn_len = q_end - q_start + 1
        seq_tree_ref[query_name] = {
            "tree_ref_name": tree_ref_name,
            "identity": identity,
            "blastn_len": blastn_len
        }
    # 参考序列元数据表转换为字典
    meta_dict = {}
    for row in meta_df.values.tolist():
        meta_dict[row[0]] = row
    seq_qc_filter_pass = {}
    seq_qc_filter_not_pass = []
    for seq_name, seq_info in input_seqs.items():
        seq_length = seq_info["length"]
        seq_seq = seq_info["seq"]
        seq_tree_ref_name = seq_tree_ref[seq_name]["tree_ref_name"]
        seq_identity = seq_tree_ref[seq_name]["identity"]
        seq_blastn_len = seq_tree_ref[seq_name]["blastn_len"]
        seq_seg_name = meta_dict[seq_tree_ref_name][11]
        seq_ref_name = ""
        for k in ref_seq.keys():
            if k.endswith(seq_seg_name):
                seq_ref_name = k
        seq_ref_seq = ref_seq[seq_ref_name]
        seq_ref_percent = seq_length / len(seq_ref_seq)
        seq_blastn_base_percent = seq_blastn_len / seq_length
        judge = (
                (seq_identity >= 90) and
                (seq_blastn_base_percent >= 0.8) and
                (seq_ref_percent >= 0.8)
        )
        if judge:
            seq_qc_filter_pass[seq_name] = {
                "seq_seq": seq_seq,
                "seq_length": seq_length,
                "seq_tree_ref_name": seq_tree_ref_name,
                "seq_identity": seq_identity,
                "seq_blastn_len": seq_blastn_len,
                "seq_seg_name": seq_seg_name,
                "seq_ref_name": seq_ref_name,
                "seq_ref_seq": seq_ref_seq,
                "seq_ref_percent": seq_ref_percent,
                "seq_blastn_base_percent": seq_blastn_base_percent
            }
        else:
            seq_qc_filter_not_pass.append(seq_name)

    return seq_qc_filter_pass, seq_qc_filter_not_pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, help='Input file')
    parser.add_argument('-b', '--blast_file', required=True, help='Blast file')
    parser.add_argument('-n', '--isolate_name', required=True, help='Isolate name')
    parser.add_argument('-w', '--work_dir', required=True, help='Isolate name')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # 获取输入序列
    params = parse_args()
    input_file = params.input_file
    blast_file = params.blast_file
    isolate_name = params.isolate_name
    work_dir = params.work_dir
    input_seqs = get_input_seq(input_file)

    # 读取元数据表
    meta_df = pd.read_csv(f"{work_dir}/h9_seg8_clade/data/all_meta_data_with_vari.csv", dtype="str",
                          keep_default_na=False)
    # 读取参考序列
    with open(f"{work_dir}/h9_seg8_clade/data/reference_seq.fasta") as file:
        ref_seq = {x.split("\n")[0]: "".join(x.split("\n")[1:])
                   for x in file.read().split(">")[1:]}
    # 获取每条序列对上的最近序列和片段名并做过滤
    seq_qc_filter_pass, seq_qc_filter_not_pass = get_input_seq_ref_seq_filter(
        input_seqs, meta_df, ref_seq, blast_file
    )

    # 生成比对变异的文件
    isolate_name_path = f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}"
    if os.path.exists(isolate_name_path):
        shutil.rmtree(isolate_name_path)
    os.makedirs(isolate_name_path)
    query_path = f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}/querys"
    os.makedirs(query_path)
    for seq_name, seq_info in seq_qc_filter_pass.items():
        seq_seq = seq_info["seq_seq"]
        seq_ref = seq_info["seq_ref_name"]
        seq_file_name = f"{query_path}/{seq_ref}___querys.fasta"
        with open(seq_file_name, "a") as f:
            f.write(">" + seq_name + "\n" + "\n".join(seq_seq) + "\n")
    annotation_isolate = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{isolate_name}"
    if os.path.exists(annotation_isolate):
        shutil.rmtree(annotation_isolate)
    os.makedirs(annotation_isolate)
    annotation_isolate_input = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{isolate_name}/input"
    os.makedirs(annotation_isolate_input)
    annotation_isolate_input_query = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{isolate_name}/input/querys"
    os.makedirs(annotation_isolate_input_query)

    source_pattern = f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}/querys/*"
    destination_dir = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{isolate_name}/input/querys"

    copy_files_recursively(source_pattern, destination_dir)

    with open(f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}/qc_filter_pass.json", "w") as f:
        json.dump(seq_qc_filter_pass, f, indent=4)

    with open(f"{work_dir}/h9_seg8_clade/script/filter/filter_result/{isolate_name}/qc_filter_not_pass.json", "w") as f:
        json.dump(seq_qc_filter_not_pass, f, indent=4)
