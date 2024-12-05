#! /usr/bin/python3.9
import os
import argparse
import shutil
from multiprocessing import Pool
import json

import pandas as pd


def parse_input(work_dir, input_file, input_isolate_name, python):
    script = f"{work_dir}/h9_seg8_clade/script/parse_input/parse_input.py"
    cmd = (f"{python} {script} -i {input_file} -n {input_isolate_name} -w {work_dir}")

    os.system(cmd)
    print("\n----------------- parse_input finished! -----------------\n")


def blastn(work_dir, python, input_directory, output_directory, blastn_path):
    script = f"{work_dir}/h9_seg8_clade/script/blastn_filter/blastn.py"
    cmd = (f"{python} {script} -i {input_directory} -o {output_directory} -bn {blastn_path} -w {work_dir}")

    os.system(cmd)
    print("\n----------------- blastn finished! ------------------\n")


def qc_and_filter(work_dir, python, parse_fasta, blast_output, input_isolate_name):
    script = f"{work_dir}/h9_seg8_clade/script/filter/qc_and_filter.py"
    cmd = (
        f"{python} {script} -i {parse_fasta} -b {blast_output} -n {input_isolate_name} -w {work_dir}")

    os.system(cmd)
    print("\n----------------- qc_and_filter finished! ------------------\n")


def call_snp_and_annotation(work_dir, python, input_isolate_name,mafft_path):
    script = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/call_snp_and_annotation.py"
    cmd = (
        f"{python} {script} -i {work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{input_isolate_name}/input -n {input_isolate_name} -wp {work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{input_isolate_name} -p {python} -m 2 -w {work_dir} -mafft {mafft_path}")
    os.system(cmd)
    print("\n----------------- call_snp_and_annotation finished! ------------------\n")


def seq_vari_new_vari(work_dir, python, input_isolate_name):
    script = f"{work_dir}/h9_seg8_clade/script/seq_vari_new/seq_vari_new_vari.py"
    cmd = (
        f"{python} {script}  -n {input_isolate_name} -w {work_dir}"
    )
    os.system(cmd)
    print("\n----------------- seq_vari_new_vari finished! ------------------\n")


def generate_output(work_dir, python, input_isolate_name,output_file):
    script = f"{work_dir}/h9_seg8_clade/script/generate_output.py"
    cmd = f"{python} {script} -n {input_isolate_name} -w {work_dir} -o {output_file}"
    os.system(cmd)
    print("\n----------------- generate_output finished! ------------------\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, help='Input file')
    parser.add_argument('-n', '--input_isolate_name', required=True, help='Input isolate name')
    parser.add_argument('-o', '--output_directory', required=True, help='output directory')
    args = parser.parse_args()
    return args


def main(work_dir, input_file, python, input_isolate_name, blastn_path,mafft_path,output_path):
    input_file = input_file.strip("[").strip("]")
    input_isolate_name = input_isolate_name.strip("[").strip("]")
    parse_input(work_dir, input_file, input_isolate_name, python)
    parse_fasta = f"{work_dir}/h9_seg8_clade/script/parse_input/parse_output/{input_isolate_name}/seqs"
    blast_output = f"{work_dir}/h9_seg8_clade/script/blastn_filter/blast_output/{input_isolate_name}"
    blastn(work_dir, python, parse_fasta, blast_output, blastn_path)

    qc_and_filter(work_dir, python, parse_fasta, blast_output, input_isolate_name)

    call_snp_and_annotation(work_dir, python, input_isolate_name,mafft_path)

    seq_vari_new_vari(work_dir, python, input_isolate_name)

    generate_output(work_dir, python, input_isolate_name,output_path)


if __name__ == "__main__":
    # the path of the work dir
    work_dir = "/data7/sunxiuqiang/virus_database"
    # the path of the python
    python = "/data/homebackup/sunxiuqiang/tools/miniconda3/bin/python"
    # the path of the blastn
    blastn_p = "/data/homebackup/sunxiuqiang/tools/blast/bin/blastn"
    # the path of the mafft
    mafft = "/data/homebackup/sunxiuqiang/tools/mafft/bin/mafft"



    params = parse_args()
    # store the input fasta path into list
    input_file = params.input_file.split(",") if ',' in params.input_file else [
        params.input_file]
    # store the input isolate name into list
    input_isolate_name = params.input_isolate_name.split(
        ",") if ',' in params.input_isolate_name else [
        params.input_isolate_name]
    output_path = params.output_directory
    # store the modified name
    input_isolate_name_modify = []
    # 生成毒株名对应的seg序列的字典
    # modify the isolate name
    for name in input_isolate_name:
        input_isolate_name_modify.append(
            name.strip().replace(" ", "_").replace("/", "_").replace("(", "").replace(")", ""))
    isolate_segment_dic = dict(zip(input_isolate_name_modify, input_file))

    # remove old result file and make new one for each run
    parse_input_path = f"{work_dir}/h9_seg8_clade/script/parse_input/parse_output"
    if os.path.exists(parse_input_path):
        shutil.rmtree(parse_input_path)
    os.makedirs(parse_input_path)

    blastn_path = f"{work_dir}/h9_seg8_clade/script/blastn_filter/blast_output"

    if os.path.exists(blastn_path):
        shutil.rmtree(blastn_path)
    os.makedirs(blastn_path)

    filter_result_path = f"{work_dir}/h9_seg8_clade/script/filter/filter_result"
    if os.path.exists(filter_result_path):
        shutil.rmtree(filter_result_path)
    os.makedirs(filter_result_path)

    call_snp_work_path = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir"
    if os.path.exists(call_snp_work_path):
        shutil.rmtree(call_snp_work_path)
    os.makedirs(call_snp_work_path)

    seq_vari_work = f"{work_dir}/h9_seg8_clade/script/seq_vari_new/work_dir"
    if os.path.exists(seq_vari_work):
        shutil.rmtree(seq_vari_work)
    os.makedirs(seq_vari_work)


    python_list = []
    work_dir_list = []
    blastn_list = []
    mafft_list = []
    output_path_list = []
    # multiprocessing to each input fasta file and isolate name
    for i in range(len(input_file)):
        python_list.append(python)
        work_dir_list.append(work_dir)
        blastn_list.append(blastn_p)
        mafft_list.append(mafft)
        output_path_list.append(output_path)
    tuples = zip(work_dir_list, input_file, python_list, input_isolate_name_modify, blastn_list,mafft_list,output_path_list)
    with Pool(processes=5) as pool:
        pool.starmap(main, tuples)

    # store the output result into json
    output_json = {}
    # each input fasta file and isolate name will get an output file, read each output path
    output_lists = os.listdir(output_path)
    # read the file of genotype, use it to identify the genotype by clades of each segment
    df_genotype = pd.read_csv(f"{work_dir}/h9_seg8_clade/data/genotype.csv",
                              keep_default_na=False,
                              low_memory=False,
                              index_col=False)
    df_genotype = df_genotype.astype(str)
    output_json_path = f"{work_dir}/h9_seg8_clade/script/output/cladeRresultJson.json"
    output_no_pass_file = f"{work_dir}/h9_seg8_clade/script/output/no_pass_number.json"
    no_pass_file = {}

    filter_results = os.listdir(filter_result_path)
    for file_name in filter_results:
        file_dir = f"{filter_result_path}/{file_name}"
        no_pass_file_path = f"{file_dir}/qc_filter_not_pass.json"
        pass_file_path = f"{file_dir}/qc_filter_pass.json"
        if file_name == "none":
            with open(no_pass_file_path) as file:
                data = json.load(file)
                if not data:
                    no_pass_file["no_pass_number"] = ""
                    no_pass_file["pass"] = "all_pass"
                else:
                    no_pass_file["no_pass_number"] = data
                    no_pass_file["pass"] = "some_pass"
            with open(pass_file_path) as file:
                data = json.load(file)
                if not data:
                    no_pass_file["no_pass_number"] = ""
                    no_pass_file["pass"] = "no_pass"

        else:
            with open(no_pass_file_path) as file:
                data = json.load(file)
                if not data:
                    no_pass_file[file_name] = {
                        "no_pass_number": "",
                        "pass": "all_pass"
                    }
                else:
                    no_pass_file[file_name] = {
                        "no_pass_number": data,
                        "pass": "some_pass"
                    }

            with open(pass_file_path) as file:
                data = json.load(file)
                if not data:
                    no_pass_file[file_name] = {
                        "no_pass_number": "",
                        "pass": "no_pass"
                    }
    with open(output_no_pass_file, 'w') as json_file:
        json.dump(no_pass_file, json_file, indent=4)

    if os.path.exists(output_json_path):
        os.remove(output_json_path)
    for output_dir in output_lists:
        output_dir_path = f"{output_path}/{output_dir}"
        output_file = f"{output_dir_path}/cladeRresultJson.json"
        if os.path.exists(output_file):
            with open(output_file, 'r', encoding='utf-8') as file:
                json_str = file.read()
                data = json.loads(json_str)
                if output_dir == "none":
                    output_json = data

                else:
                    seg_clade = {}
                    for key, value in data.items():
                        seg_clade[value["segment"]] = value["cladeName"]
                    df = pd.DataFrame([seg_clade])
                    df = df.astype(str)

                    segs = df.columns.unique().tolist()
                    df_merge = pd.merge(df, df_genotype, on=segs, how="left")
                    df_merge["Genotype"] = df_merge["Genotype"].fillna('No genotype is matched')
                    genotype = df_merge["Genotype"].tolist()
                    if len(segs) == 8:
                        segment_complete = "complete"
                    else:
                        segment_complete = "partial"
                    output_json[output_dir] = {
                        "genotype": genotype,
                        "segment_complete": segment_complete,
                        "segment_count": len(segs),
                        "segment": segs,
                        "data": data
                    }
    with open(output_json_path, 'w') as json_file:
        json.dump(output_json, json_file, indent=4)
