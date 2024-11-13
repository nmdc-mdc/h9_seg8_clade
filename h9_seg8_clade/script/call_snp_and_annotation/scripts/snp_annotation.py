# -*- coding: UTF-8 -*-
'''Author:Kesheng Peng
e-mill:kesheng_peng_yx@163.com
Version:1.0
release-date:2023-10-19'''


import os
import json
import copy
import shutil
import argparse
import multiprocessing
import pandas as pd
from multiprocessing import Pool


def parse_args():
    """ Get command line parameters"""
    parser = argparse.ArgumentParser(description="SNP & annotation!")
    parser.add_argument("-p",
					"--praperation_path",
					 required=True)
    parser.add_argument("-s",
                        "--differ_site_file_path",
                        required=True)
    parser.add_argument("-w",
                        "--work_path",
                        required=True)
    parser.add_argument("-m",
                        "--multiprocessing_number",
                        help="The number of multiprocessing!")
    parser.add_argument("-g",
                        "--genetic_amino_code",
                        help="Nucleotide coding table!")
    args = parser.parse_args()

    return args


def get_annotation_meta(praperation_path):
    """ Reference sequence meta-information when obtaining SNP annotations """

    with open(f"{praperation_path}/non_cds_site.json") as file:
        non_cds_base_dict = json.load(file)

    with open(f"{praperation_path}/base_amino_site.json") as file:
        base_amino_dict = json.load(file)

    return non_cds_base_dict,base_amino_dict


def annotation(paths,metas):
    """ 
    Annotate the information of differential 
    loci using the reference sequence metadata dictionary.
    Read the difference sites, filter, 
    and extract the query sequence and reference sequence, 
    generate snps, insertions, and deletions, 
    and annotate them using the dictionary from step 2. 
    The annotation results are separately entered into annotation_ In the result folder"""

    differ_site_file = paths[0]
    annotation_result_file = paths[1]
    non_cds_base_dict = copy.deepcopy(metas[0])
    base_amino_dict = copy.deepcopy(metas[1])
    # genetic_amino_code_dict = metas[2]

    # print(differ_site_file,annotation_result_file,len(non_cds_base_dict),len(base_amino_dict))

    # Nucleotide coding table, 
    # data from https://www.mathworks.com/help/bioinfo/ref/revgeneticcode.html
    genetic_amino_code = {
        "GCT": "A","GCC": "A","GCA": "A","GCG": "A","CGT": "R","CGC": "R","CGA": "R","CGG": "R",
        "AGA": "R","AGG": "R","AAT": "N","AAC": "N","GAT": "D","GAC": "D","TGT": "C","TGC": "C",
        "CAA": "Q","CAG": "Q","GAA": "E","GAG": "E","GGT": "G","GGC": "G","GGA": "G","GGG": "G",
        "CAT": "H","CAC": "H","ATT": "I","ATC": "I","ATA": "I","TTA": "L","TTG": "L","CTT": "L",
        "CTC": "L","CTA": "L","CTG": "L","AAA": "K","AAG": "K","ATG": "M","TTT": "F","TTC": "F",
        "CCT": "P","CCC": "P","CCA": "P","CCG": "P","TCT": "S","TCC": "S","TCA": "S","TCG": "S",
        "AGT": "S","AGC": "S","ACT": "T","ACC": "T","ACA": "T","ACG": "T","TGG": "W","TAT": "Y",
        "TAC": "Y","GTT": "V","GTC": "V","GTA": "V","GTG": "V","TAA": "*","TAG": "*","TGA": "*"
    }

    def is_none(var):  
        return isinstance(var, type(None))
    # if is_none(genetic_amino_code_dict):
    #     genetic_amino_code = genetic_amino_code_dict

    # Read the differential site file, 
    # including the position and base of the reference sequence, 
    # query the base and position of the sequence, 
    # and the reference sequence name (ID number), 
    # query the sequence name (ID number).
    # The position of the reference sequence 
    # and the position of the query sequence are both arranged in ascending order.
    df = pd.read_csv(differ_site_file)
    if df.shape[0] != 0:
        variation_df = copy.deepcopy(df)
        first_row = list(df.iloc[0])
        reference_name = first_row[4]
        reference_dict = copy.deepcopy(base_amino_dict[reference_name])
        
        # Create a Boolean index to find out the rows that meet the conditions. 
        # The reference sequence and query sequence are not all "-" numbers. 
        # Some alignment results are "-" to indicate insertion or deletion.
        mask = (df["ref"].str.contains("-")==False) & (df["query"].str.contains("-")==False)               
        # -----------SNP-----------------
        # Put the body conditions into the table, snp_ df
        snp_df = df[mask] 
        # snp_ Convert to a dictionary, 
        # using the position of the reference sequence as the key, 
        # snp_ The df row is converted to a list as a value.                                                                                             
        snp_dict = {}                                                                                             
        for index,row in snp_df.iterrows():
            row = list(row)
            ref_site = row[0]
            snp_dict[str(ref_site)] = row
        # List of annotation results for snp
        all_snp_result = []                                                                                             
        for ref_key,ref_value in reference_dict.items():
            # ref_key is cds name
            # All base sites involved in the coding region.
            cds_base_site = []                                                                                          
            for ref_site in ref_value:
                cds_base_site.append(int(ref_site))
            cds_base_site = [str(x) for x in sorted(cds_base_site)]
            key_index = []                   
            # Used to read the position of the first base of 
            # each codon in the coding region in the reference sequence.                                                                       
            index_number = 0
            for index in range(0,int(len(cds_base_site)/3)):
                key_index.append(cds_base_site[index_number])
                index_number = index_number + 3
            for code_index in key_index:
                ref_code_dict = copy.deepcopy(ref_value[code_index])
                # Amino acid in reference sequence.
                pre_amino = ref_code_dict["amino"] 
                  
                # The name of the coding region, usually the name of the protein.                                                                   
                cds_name = ref_code_dict["cds_name"] 
                # Gene name.                                                                  
                gene = ref_code_dict["gene"] 
                # Position of reference sequence amino acid.                                                                           
                amino_site = ref_code_dict["amino_site"]                                                                
                exists_site_base = []
                # Whether the first base of the codon exists in the query variation.
                # if code_index in snp_dict:                                                                              
                #     query_row = snp_dict[code_index]
                #     exists_site_base.append(query_row)
                #     ref_code_dict["three_base"][0] = query_row[2]
                # Whether the second base of the codon exists in the query variation.
                if str(int(code_index)+1) in snp_dict:                                                                            
                    query_row = snp_dict[str(int(code_index)+1)]
                    exists_site_base.append(snp_dict[str(int(code_index)+1)])                                                     
                    ref_code_dict["three_base"][1] = query_row[2]
                # Whether the third base of the codon exists in the query mutation.
                if str(int(code_index)+2) in snp_dict:                                                                                                                                                       
                    query_row = snp_dict[str(int(code_index)+2)]
                    exists_site_base.append(snp_dict[str(int(code_index)+2)])
                    ref_code_dict["three_base"][2] = query_row[2]
                # If the query sequence base site exists in the reference codon dictionary.
                if exists_site_base != []:                                                                             
                    three_base = "".join(ref_code_dict["three_base"])
                    post_amino = genetic_amino_code[three_base]
                    # print(pre_amino,post_amino)
                    # If the amino acid after mutation is consistent with that before mutation.
                    if pre_amino == post_amino:                                                                         
                        if post_amino != "*":
                            annotaion_type = "synonymous SNV"
                            for row in exists_site_base:
                                all_text = [row[-1],row[4],gene] + row[0:4] + ["snp"] + [cds_name,amino_site,pre_amino,post_amino,"synonymous SNV"]
                                all_snp_result.append(all_text)
                    # If the amino acid after mutation is inconsistent with that before mutation.
                    else:                                                                                               
                        if pre_amino != "*":
                            if post_amino != "*":
                                annotaion_type = "nonsynonymous SNV"
                                for row in exists_site_base:
                                    all_text = [row[-1],row[4],gene] + row[0:4] + ["snp"] + [cds_name,amino_site,pre_amino,post_amino,"nonsynonymous SNV"]
                                    all_snp_result.append(all_text)
                            else:
                                anntation_type = "stopgain"
                                for row in exists_site_base:
                                    all_text = [row[-1],row[4],gene] + row[0:4] + ["snp"] + [cds_name,amino_site,pre_amino,post_amino,"stopgain"]
                                    all_snp_result.append(all_text)
                        else:
                            annotaion_type = "stoploss"
                            for row in exists_site_base:
                                    all_text = [row[-1],row[4],gene] + row[0:4] + ["snp"] + [cds_name,amino_site,pre_amino,post_amino,"stoploss"]
                                    all_snp_result.append(all_text)
        if all_snp_result != []:
            snp_result_df = pd.DataFrame(all_snp_result)
            snp_result_df.to_csv(annotation_result_file,mode="a",index=None,header=None) 

        # Put the rows where the reference sequence or query sequence base is "." into the table, indel_df.
        # Table containing base insertions and deletions.
        indel_df = df[~mask]
        # Base insertion dictionary, 
        # the key is the position of the reference sequence, 
        # and the value is a list of several base insertions at this position    
        insertion_dict = {}
        # The base deletion dictionary, 
        # the key is the position of the reference sequence, 
        # and the value is a list of the base deletions at this position.     
        deletion_dict = {}      
        for index,row in indel_df.iterrows():
            row = list(row)
            ref_site = row[0]
            ref_base = row[1]
            query_site = row[3]
            query_base = row[2]
            if ref_base=="-": 
                if int(ref_site) not in insertion_dict:
                    insertion_dict[int(ref_site)] = [row]
                elif int(ref_site) in insertion_dict:
                    insertion_dict[int(ref_site)].append(row)
            if query_base=="-":
                if int(query_site) not in deletion_dict:
                    deletion_dict[int(query_site)] = [row]
                elif int(query_site) in deletion_dict:
                    deletion_dict[int(query_site)].append(row)
        # -----Filter the Insertion and convert it to a dictionary. 
        # Insert 3 or 6 bases and leave the rest filtered out.------
        insert_filter_dict = {}
        for k,v in insertion_dict.items():
            # If the subsequent shift mutation needs to be considered,
            #  then this judgment should not be made 
            if len(v)==3 or len(v)==6:     
                ref_site = k
                ref_base = "-"
                query_site = ""
                query_bases = []
                ref = v[0][4]
                query = v[0][5]
                for q in v:
                    query_bases.append(q[2])
                query_bases = "".join(query_bases)
                insert_filter_dict[ref_site] = [ref_site,ref_base,query_bases,query_site,ref,query]

        # The annotation of insertion
        all_insertion_result = []
        for insertion_k,insertion_v in insert_filter_dict.items():
            for protein_name,pro_value in reference_dict.items():
                if str(insertion_k) in pro_value:
                    pro_list = pro_value[str(insertion_k)]
                    gene = pro_list['gene']
                    cds_name = pro_list['cds_name']
                    amino_site = pro_list['amino_site']
                    amino = pro_list['amino']
                    index_start = 0
                    insertion_base = insertion_v[2]
                    insertion_amino_list = []
                    for index in range(0,len(insertion_base)//3):
                        bases = insertion_base[index_start:index_start+3]
                        if bases in genetic_amino_code:
                            insertion_amino_list.append(genetic_amino_code[bases])
                            index_start = index_start + 3 
                    if insertion_amino_list != []:
                        insertion_amino = "".join(insertion_amino_list)
                        if "*" not in insertion_amino:
                            annotaion_type = "insertion"
                            all_text = ([insertion_v[-1],insertion_v[4],gene,insertion_v[0],insertion_v[1],insertion_v[2],insertion_v[3]] + 
                                        ["insertion",cds_name,amino_site,"-",insertion_amino,"insertion"])
                            all_insertion_result.append(all_text)
                        else:
                            annotaion_type = "stopgain"
                            all_text = ([insertion_v[-1],insertion_v[4],gene,insertion_v[0],insertion_v[1],insertion_v[2],insertion_v[3]] + 
                                        ["insertion",cds_name,amino_site,"-",insertion_amino,"stopgain"])
                            all_insertion_result.append(all_text)
        if all_insertion_result != []:
            insertion_result_df = pd.DataFrame(all_insertion_result)
            insertion_result_df.to_csv(annotation_result_file,mode="a",index=None,header=None)

        # Filter Deletion and convert it to a dictionary. 
        # Filter out those with 3 or 6 bases missing and keep the rest.
        deletion_filter_dict = {}
        for k_d,v_d in deletion_dict.items():
            # If there is a need to consider the subsequent shift mutation, 
            # then this judgment should not be made.
            if len(v_d)==3 or len(v_d)==6:      
                ref_site_d = v_d[0][0]
                ref_bases_d = []
                query_site_d = ""
                query_base_d = "-"
                ref_d = v_d[0][4]
                query_d = v_d[0][5]
                for q_d in v_d:
                    ref_bases_d.append(q_d[1])
                ref_bases_d = "".join(ref_bases_d)
                deletion_filter_dict[ref_site_d] = [ref_site_d,ref_bases_d,query_base_d,query_site_d,ref_d,query_d]
        # For the missing annotation, 
        # the leftmost position of the reference sequence that is missing is considered as the missing position.
        all_deletion_result = []
        for deletion_k,deletion_v in deletion_filter_dict.items():
            for protein_name,pro_value in reference_dict.items():
                if str(deletion_k) in pro_value:
                    pro_list = pro_value[str(deletion_k)]
                    gene = pro_list['gene']
                    cds_name = pro_list['cds_name']
                    amino_site = pro_list['amino_site']
                    amino = pro_list['amino']
                    index_start = 0
                    deletion_base = deletion_v[1]
                    deletion_amino = []
                    for index in range(0,len(deletion_base)//3):
                        bases = deletion_base[index_start:index_start+3]
                        if bases in genetic_amino_code:
                            deletion_amino.append(genetic_amino_code[bases])
                            index_start = index_start + 3 
                    if deletion_amino != []:
                        deletion_amino = "".join(deletion_amino)
                        if "*" not in deletion_amino:
                            annotaion_type = "deletion"
                            all_text = ([deletion_v[-1],deletion_v[4],gene,deletion_v[0],deletion_v[1],deletion_v[2],deletion_v[3]] + 
                                    ["deletion",cds_name,amino_site,deletion_amino,"-","deletion"])
                            all_deletion_result.append(all_text)
                        else:
                            annotaion_type = "stoploss"
                            all_text = ([deletion_v[-1],deletion_v[4],gene,deletion_v[0],deletion_v[1],deletion_v[2],deletion_v[3]] + 
                                ["deletion",cds_name,amino_site,deletion_amino,"-","stoploss"])
                            all_deletion_result.append(all_text)
        if all_deletion_result != []:
            deletion_result_df = pd.DataFrame(all_deletion_result)
            deletion_result_df.to_csv(annotation_result_file,mode="a",index=None,header=None)

        # Note that the variation in the non-coding region is mainly used to calculate the effectiveness of primers later.
        non_cds_base_dict = non_cds_base_dict
        non_cds_snp_result = []
        for index,row in variation_df.iterrows():
            row = list(row)
            reference = row[4]
            ref_site = row[0]
            ref_base = row[1]
            query_base = row[2]
            query_site = row[3]
            if reference in non_cds_base_dict:
                reference_53_dict = non_cds_base_dict[reference]
                if str(ref_site) in reference_53_dict:
                    ref_53_site = reference_53_dict[str(ref_site)][0]
                    # Only snps are retained here, regardless of insertions and deletions
                    if "-" not in ref_base and "-" not in query_base:
                        text = [row[-1]] + [reference,ref_53_site[1],ref_53_site[2],ref_base,query_base,query_site] + ["snp","","","",""] 
                        non_cds_snp_result.append(text)
                    # elif "-" in ref_base and "-" not in query_base:
                    #     text = [row[-1]] + [reference,ref_53_site[1],ref_53_site[2],ref_base,query_base] + ["insertion","","","",""] 
                    # elif "-" not in ref_base and "-" in query_base:
                    #     text = [row[-1]] + [reference,ref_53_site[1],ref_53_site[2],ref_base,query_base] + ["deletion","","","",""]
        if non_cds_snp_result != []:                
            non_cds_df = pd.DataFrame(non_cds_snp_result)
            non_cds_df.to_csv(annotation_result_file,mode="a",index=None,header=None)


def multiprocessing_annotation(work_path,
                               non_cds_base_dict,
                               base_amino_dict,
                               multiprocessing_number,
                            #    genetic_amino_code
                               ):
    """ Annotate the information of differential 
    loci using the reference sequence metadata dictionary """

    differ_site_path = f"{work_path}/differ_site"
    annotation_result_path = f"{work_path}/annotation_result"
    if os.path.exists(annotation_result_path):
        shutil.rmtree(annotation_result_path)
    os.mkdir(annotation_result_path)

    paths = [
            [f"{differ_site_path}/{x}",
             f"{annotation_result_path}/{x.split('.')[0]}.csv"]
             for x in os.listdir(differ_site_path)]
    metas = [[non_cds_base_dict,
              base_amino_dict,
            #   genetic_amino_code
              ]]*len(paths)

    tuples = zip(paths,metas)

    def is_none(var):  
        return isinstance(var, type(None))

    if is_none(multiprocessing_number):
        with Pool(processes=multiprocessing_number) as pool:
            pool.starmap(annotation,tuples)
    else:
        multiprocessing_number = multiprocessing.cpu_count()
        with Pool(processes=multiprocessing_number) as pool:
            pool.starmap(annotation,tuples)


if __name__ == "__main__":
    # Get command lines.
    params = parse_args()
    praperation_path = params.praperation_path
    differ_site_file_path = params.differ_site_file_path
    work_path = params.work_path
    multiprocessing_number = params.multiprocessing_number
    # genetic_amino_code = params.genetic_amino_code

    # Reference sequence meta-information when obtaining SNP annotations
    ref_meta = get_annotation_meta(praperation_path)
    non_cds_base_dict = ref_meta[0]
    base_amino_dict = ref_meta[1]

    # snp annotation
    multiprocessing_annotation(work_path,
                               non_cds_base_dict,
                               base_amino_dict,
                               multiprocessing_number,
                            #    genetic_amino_code
                               )