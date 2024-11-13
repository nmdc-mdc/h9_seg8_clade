# -*- coding: UTF-8 -*-
'''Author:Kesheng Peng
e-mill:kesheng_peng_yx@163.com
Version:1.0
release-date:2023-10-19'''


import os
import shutil
import argparse



def parse_args():
    """ Get command line parameters"""
    parser = argparse.ArgumentParser(description="SNP & annotation!")
    parser.add_argument(
                        "-w",
                        "--work_path",
                        required=True)
    args = parser.parse_args()

    return args


def delete_intermediate_files(work_path):
    """ Delete intermediate files """
    praperation_path = f"{work_path}/praperation"
    if os.path.exists(praperation_path):
        shutil.rmtree(praperation_path)

    differ_site_path = f"{work_path}/differ_site"
    if os.path.exists(differ_site_path):
        shutil.rmtree(differ_site_path)

    annotation_result_path = f"{work_path}/annotation_result"
    if os.path.exists(annotation_result_path):
        shutil.rmtree(annotation_result_path)

    
if __name__ == "__main__":
    params = parse_args()
    work_path = params.work_path

    delete_intermediate_files(work_path)
