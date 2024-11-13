import os
import argparse
import shutil


def multi_seq_blastn(input_directory, output_directory, blastn, work_dir):
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
    os.makedirs(output_directory, exist_ok=True)
    input_files = [f for f in os.listdir(input_directory) if f.endswith('.fasta')]
    database = f"{work_dir}/h9_seg8_clade/data/all_ref_seqs/all_ref_seqs"
    for input_file in input_files:
        input_file_name = input_file.split('.')[0]

        input_path = f"{input_directory}/{input_file}"
        out_file_path = f"{output_directory}/{input_file_name}.txt"

        cmd = (
            f"{blastn} -query {input_path} -db {database} -out {out_file_path} -outfmt 6"
        )

        os.system(cmd)


def parse_args():
    parser = argparse.ArgumentParser(description='blastn')
    parser.add_argument('-i', '--input_directory', type=str, required=True)
    parser.add_argument('-o', '--output_directory', type=str, required=True)
    parser.add_argument('-bn', '--blastn', type=str, required=True)
    parser.add_argument('-w', '--work_dir', type=str, required=True)
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    params = parse_args()
    multi_seq_blastn(params.input_directory, params.output_directory, params.blastn, params.work_dir)
