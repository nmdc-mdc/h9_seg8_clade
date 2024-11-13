import os
import re
import json
import shutil
import argparse
from multiprocessing import Pool


def parse_input(input_file, isolate_name,work_dir):
    """ Split sequence files and generate sequence's filenames. """
    output_dir = f"{work_dir}/h9_seg8_clade/script/parse_input/parse_output/{isolate_name}"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(output_dir)
    with open(input_file) as file:
        fastas = file.read().strip().split(">")[1:]
    seq_name_dict = {}
    seqs = {}
    for index, fasta in enumerate(fastas):
        lines = fasta.split("\n")
        # 最开始的序列名，用于最后展示到表格中
        oringin_name = lines[0]
        # 序列名，用于计算过程中的序列名，也是后面展示到进化树中的名字，不能有括号
        name = lines[0].split()[0].replace("(", ".").replace(")", ".")
        # 替换序列中特殊的字符，这样保证mafft能处理
        seq = []
        for line in lines[1:]:
            line = re.sub(r"[^ATCGU]", "U", line.upper())
            seq.append(line)
        seq = "\n".join(seq)
        name_seq = ">" + name + "\n" + seq
        # 每条序列对应一个文件用于后续比对的处理，序列名变为数字加下划线
        seq_file_name = "_" + str(index) + "_"
        seq_name_dict[seq_file_name] = {
            "oringin_name": oringin_name,
            "name": name
        }
        seqs[seq_file_name] = name_seq

    # 后续从这个文件中获取文件名和序列名之间的对应关系
    with open(f"{output_dir}/seq_name_dict.json", "w") as file:
        json.dump(seq_name_dict, file, indent=4)

    # 拆分后的序列文件
    seqs_path = f"{output_dir}/seqs"
    os.makedirs(seqs_path)
    for name, seq in seqs.items():
        with open(f"{seqs_path}/{name}.fasta", "w") as file:
            file.write(seq)


def parse_args():
    parser = argparse.ArgumentParser(description='Parse input file')
    parser.add_argument('-i', '--input_file', required=True, help='Input file')
    parser.add_argument('-n', '--isolate_name', required=True, help='Isolate name')
    parser.add_argument('-w', '--work_dir', required=True, help='Work dir path')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    input_file = args.input_file
    isolate_name = args.isolate_name
    work_dir = args.work_dir
    parse_input(input_file, isolate_name, work_dir)
