# -*- coding: UTF-8 -*-
'''Author:Kesheng Peng
e-mill:keseheng_peng_yx@163.com
Version:1.0
release-date:2023-10-19'''

import os
import argparse
import shutil


def parse_args():
    """ Get command line parameters """
    parser = argparse.ArgumentParser(description="call snp & annotation preparation of input data!")
    parser.add_argument("-i",
                        "--input_file_path",
                        required=True, help="Input folder path must be a full path!")
    parser.add_argument("-wp",
                        "--work_path",
                        required=True,
                        help="The path where the generated intermediate files and result files are stored")
    parser.add_argument("-p",
                        "--python",
                        required=True)
    parser.add_argument("-m",
                        "--multiprocessing_number",
                        help="The number of multiprocessing!")
    parser.add_argument("-g",
                        "--genetic_amino_code",
                        help="Nucleotide coding table!")
    parser.add_argument("-d",
                        "--delete",
                        )
    parser.add_argument("-n",
                        "--isolate_name",
                        )
    parser.add_argument("-w",
                        "--work_dir",
                        )
    parser.add_argument("-mafft",
                        "--mafft_path",
                        )
    args = parser.parse_args()
    return args


def run_program(input_file_path, work_path, multiprocessing_number, delete, current_path, python, isolate_name,
                work_dir,mafft_path):
    """ Run programs """
    script_path = f"{current_path}/scripts"

    def is_none(var):
        return isinstance(var, type(None))



    reference_file_path = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/input"
    praperation_path = f"{script_path}/input_praperation.py"

    input_praperation_cmd = python + " " + praperation_path + " -i " + reference_file_path + " -n " + isolate_name + " -wp " + work_path + " -w " + work_dir

    os.system(input_praperation_cmd)

    print("\nPreparation completed.\n")

    differ_path = f"{script_path}/align_differ_site_snp.py"
    query_ref_file_path = f"{work_path}/praperation"

    query_ref_file_path = f"{query_ref_file_path}/query_seqs"
    differ_cmd = python + " " + differ_path + " -p " + query_ref_file_path + " -w " + work_path + " -m " + multiprocessing_number + " -mafft " +mafft_path
    os.system(differ_cmd)

    print("Sequence alignment and identification of differential positions completed.\n")

    annotation_path = f"{script_path}/snp_annotation.py"
    differ_site_file_path = f"{work_path}/differ_site"
    query_ref_file_path = f"{work_path}/praperation"
    annotation_cmd = python + " " + annotation_path + " -p " + query_ref_file_path + " -w " + work_path + " -s " + differ_site_file_path + " -m " + multiprocessing_number
    os.system(annotation_cmd)

    print("Amino acid variation annotation completed.\n")

    merge_path = f"{script_path}/merger_result.py"
    merge_cmd = python + " " + merge_path + " -w " + work_path
    os.system(merge_cmd)

    print("Merging of results completed.\n")

    if is_none(delete) and delete == "delete":
        delete_path = f"{script_path}/delete_intermediate_files.py"
        delete_cmd = python + " " + delete_path + " -w " + work_path
        os.system(delete_cmd)

        print("Deletion of intermediate files completed.\n")


if __name__ == "__main__":


    params = parse_args()
    input_file_path = params.input_file_path
    isolate_name = params.isolate_name
    work_path = params.work_path
    work_dir = params.work_dir
    current_path = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation"
    multiprocessing_number = params.multiprocessing_number
    delete = params.delete
    python = params.python
    mafft_path = params.mafft_path
    run_program(input_file_path, work_path, multiprocessing_number, delete, current_path, python, isolate_name,
                work_dir,mafft_path)
