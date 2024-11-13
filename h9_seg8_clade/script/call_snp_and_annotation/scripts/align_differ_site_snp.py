# -*- coding: UTF-8 -*-
'''Author:Kesheng Peng
e-mill:keseheng_peng_yx@163.com
Version:1.0
release-date:2023-10-19'''

import pandas as pd
import os
import shutil
import argparse
import multiprocessing
from multiprocessing import Pool


def parse_args():
    """ Get command line paramiters """
    parser = argparse.ArgumentParser(description="call snp & annotation. align --> differ site --> snp !")
    parser.add_argument("-p",
                        "--query_ref_file_path",
                        required=True,
                        help="The paths for storing the query sequence and reference sequence separately!")
    parser.add_argument("-w",
                        "--work_path",
                        required=True,
                        help="The path where the generated intermediate files and result files are stored!")
    parser.add_argument("-m",
                        "--multiprocessing_number",
                        help="The number of multiprocessing!")
    parser.add_argument("-mafft",
                        "--mafft_path",
                        help="mafft path")
    args = parser.parse_args()

    return args


def align_seqs(input_file, result_file, mafft):
    """ Alignment of query sequences and reference sequences. """
    # Alignment of query sequences and reference sequences using the default options of MAFFT.
    command = f"{mafft} --quiet --thread -1 {input_file} > {result_file}"
    os.system(command)


def multiprocessing_align_seqs(query_ref_file_path, work_path, multiprocessing_number, mafft_path):
    """ Alignment of reference sequences and query sequences. """
    align_result_path = f"{work_path}/align_result"
    if os.path.exists(align_result_path):
        shutil.rmtree(align_result_path)
    os.mkdir(align_result_path)
    mafft_path_list = []
    query_ref_files = os.listdir(query_ref_file_path)
    input_file_paths = [f"{query_ref_file_path}/{x}" for x in query_ref_files]
    result_file_paths = [f"{align_result_path}/{x}" for x in query_ref_files]
    for i in range(len(result_file_paths)):
        mafft_path_list.append(mafft_path)
    tuples = zip(input_file_paths, result_file_paths, mafft_path_list)

    def is_none(var):
        return isinstance(var, type(None))

    if is_none(multiprocessing_number):
        multiprocessing_number = int(multiprocessing_number)
        with Pool(processes=multiprocessing_number) as pool:
            pool.starmap(align_seqs, tuples)
    else:
        multiprocessing_number = multiprocessing.cpu_count()
        with Pool(processes=multiprocessing_number) as pool:
            pool.starmap(align_seqs, tuples)


def get_differ_site(input_file, output_file):
    """ Align result to differ site """
    # Read the alignment result file of two sequences,
    # with the reference sequence following and the query sequence preceding.
    with open(input_file) as f:
        mafft_result = f.read().split(">")[1:]
    # The names and sequences of the two sequences,
    # with the first sequence being the query sequence
    # and the second sequence being the reference sequence
    two_seq_names = [x.split("\n")[0] for x in mafft_result]
    two_seqs = ["".join(x.split("\n")[1:]) for x in mafft_result]
    # Convert the sequences into lists,
    # with the first sequence being the query sequence
    # and the second sequence being the reference sequence.
    try:
        query_list = list(two_seqs[0])
        ref_list = list(two_seqs[1])
    except Exception as e:
        print(e)
        print(mafft_result)
    # Process the location information of query
    query_positions = []
    pos_query = 0
    for base in query_list:
        if base == '-' or base == ' ':
            query_positions.append(pos_query)
        else:
            pos_query += 1
            query_positions.append(pos_query)
    # Process the location information of reference
    ref_positions = []
    pos_ref = 0
    for base in ref_list:
        if base == '-' or base == ' ':
            ref_positions.append(pos_ref)
        else:
            pos_ref += 1
            ref_positions.append(pos_ref)

    reference_name = []
    for r in range(0, len(ref_positions)):
        reference_name.append(two_seq_names[1])
    query_name = []
    for r in range(0, len(ref_positions)):
        query_name.append(two_seq_names[0])
    # Convert the sequence and position information to a Pandas table
    df = pd.DataFrame({
        'ref_site': ref_positions,
        'ref': [x.upper() for x in ref_list],
        'query': [x.upper() for x in query_list],
        'query_site': query_positions,
        'ref_name': reference_name,
        'query_name': query_name
    })
    # df.to_csv("./differ_site/pre_duplicate_" + file_name + "_1.csv",mode="a",encoding="utf-8",index=None)
    # Get the position of the first
    # and last bases in the alignment length of the reference sequence and query sequence
    query_start = query_list.index(two_seqs[0].strip("-")[0])
    query_end = len(query_list) - query_list[::-1].index(two_seqs[0].strip("-")[-1]) - 1
    ref_start = ref_list.index(two_seqs[1].strip("-")[0])
    ref_end = len(ref_list) - ref_list[::-1].index(two_seqs[1].strip("-")[-1]) - 1
    # cut off the head and tail of the table
    # print(ref_start,query_start)
    # print(ref_end,query_end)
    cut_table_row_start = max(ref_start, query_start)
    cut_table_row_end = min(ref_end, query_end)
    df = df.loc[cut_table_row_start:cut_table_row_end]
    # df.to_csv("./differ_site/pre_duplicate_" + file_name + "_2.csv",mode="a",encoding="utf-8",index=None)
    # Remove the row in the table where the ref and query sequences are the same.
    df = df[df.iloc[:, 1] != df.iloc[:, 2]]
    valid_chars = ['A', 'G', 'C', 'T', '-']
    # Select the rows in which the values in any of the 'ref' and 'query' columns are not in the valid character list.
    invalid_rows = df.loc[~df['ref'].isin(valid_chars) | ~df['query'].isin(valid_chars)]
    # 删除这些行
    df.drop(invalid_rows.index, inplace=True)
    df.to_csv(output_file, mode="a", encoding="utf-8", index=None)


def multiprocessing_differ_site(work_path, multiprocessing_number):
    """ Align result to differ site """
    align_reuslt_path = f"{work_path}/align_result"
    differ_result_path = f"{work_path}/differ_site"
    if os.path.exists(differ_result_path):
        shutil.rmtree(differ_result_path)
    os.mkdir(differ_result_path)

    align_result_files = [f"{align_reuslt_path}/{x}" for x in os.listdir(align_reuslt_path)]

    differ_result_files = [f"{differ_result_path}/{x.split('.')[0]}.csv" for x in
                           os.listdir(align_reuslt_path)]

    tuples = zip(align_result_files, differ_result_files)

    def is_none(var):
        return isinstance(var, type(None))

    if is_none(multiprocessing_number):
        multiprocessing_number = int(multiprocessing_number)
        with Pool(processes=multiprocessing_number) as pool:
            pool.starmap(get_differ_site, tuples)
    else:
        multiprocessing_number = multiprocessing.cpu_count()
        with Pool(processes=multiprocessing_number) as pool:
            pool.starmap(get_differ_site, tuples)


if __name__ == "__main__":
    params = parse_args()
    query_ref_file_path = params.query_ref_file_path
    work_path = params.work_path
    mafft_path = params.mafft_path
    multiprocessing_number = params.multiprocessing_number

    # Alignment of reference sequences and query sequences.
    multiprocessing_align_seqs(query_ref_file_path, work_path, multiprocessing_number, mafft_path)

    # Align result to differ site
    multiprocessing_differ_site(work_path, multiprocessing_number)
