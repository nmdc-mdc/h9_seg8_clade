# -*- coding: UTF-8 -*-
'''Author:Kesheng Peng
e-mill:keseheng_peng_yx@163.com
Version:1.0
release-date:2023-10-19'''


import os
import re
import json
import shutil
import argparse
import pandas as pd 
from Bio.Seq import Seq
 

def parse_args():
    """ Get command line parameters """
    parser = argparse.ArgumentParser(description="call snp & annotation preparation of input data!")
    parser.add_argument("-i",
                        "--input_file_path",
                        required=True,help="Input folder path must be a full path!")
    parser.add_argument("-n",
                        "--isolate_name",
                        required=True, help="Isolate name!")
    parser.add_argument("-wp",
                        "--work_path",
                        required=True,
                        help="The path where the generated intermediate files and result files are stored")
    args = parser.parse_args()
    parser.add_argument("-w",
                        "--work_dir",
                        required=True,
                        help="The path where the generated intermediate files and result files are stored")
    args = parser.parse_args()
    return args

 
def praperation_directory(work_path):
    """ Create a directory to store the generated input files. """ 
    praperation_directory_path = f"{work_path}/praperation"
    print(praperation_directory_path)
    if os.path.exists(praperation_directory_path):
        shutil.rmtree(praperation_directory_path)
    os.mkdir(praperation_directory_path)

    return praperation_directory_path


def split_reference_files(input_file_path,work_path):
    """ Split the reference sequence into individual 
    files and generate a list of reference sequence names. """

    reference_file_path = f"{input_file_path}/reference_seq.fasta"
    references_directory_path = f"{praperation_directory_path}/references"
    if os.path.exists(references_directory_path):
        shutil.rmtree(references_directory_path)
    os.mkdir(references_directory_path)
    reference_list_path = f"{praperation_directory_path}/reference_list.txt"
    if os.path.exists(reference_list_path):
        os.remove(reference_list_path)

    with open(reference_file_path) as file:
        refs = file.read().strip().split(">")[1:]
    
    ref_names = []
    for ref in refs:
        ref_list = ref.split("\n")
        ref_name = ref_list[0]
        ref_names.append(ref_name)
        with open(f"{references_directory_path}/{ref_name}.fasta",mode="a") as file:
            file.write(">"+ref)

    with open(reference_list_path,mode="a") as file:
        file.write("\n".join(ref_names))
    

def gene_site_base(input_file_path):
    """ Read the meta information of the non-coding region and UTR 3 5`
      of the gene to convert the reference sequence into a dictionary """
    gene_region_file = f"{input_file_path}/gene_regions.txt"
    references_file = f"{input_file_path}/reference_seq.fasta"
    gene_region = pd.read_csv(gene_region_file,
                              dtype="str",
                              keep_default_na=False,
                              sep="\t")
    # gene_region_title = gene_region.columns.to_list() # ['reference', 'start_ends', 'gene/locus_tag']
    # print(gene_region_title)
    # Read the reference sequence file and convert it to a dictionary, 
    # where the key is the sequence name and the value is the sequence
    with open(references_file) as file:                   
        reference_seq = {x.split("\n")[0]:"".join(x.split("\n")[1:]) for x in file.read().strip().split(">")[1:]}
    # Generate a table of locations, bases, and corresponding gene names
    result = []
    for row in gene_region.values:
        row = list(row)
        reference = row[0]
        region = row[1].split(";")
        gene = row[2]
        regions = []
        for r in region:
            for x in r.split("-"):
                regions.append(int(x))
        regions = sorted(regions)
        start = regions[0]
        end = regions[1]
        # get sequence
        seq = reference_seq[reference][start-1:end]
        for base in seq:
            result_row = [reference,gene,start,base]
            result.append(result_row)
            start = start + 1
    title = ["reference","gene","ref_site","ref_nucleicacid"]
    df = pd.DataFrame(result,columns=title)

    return df


def base_amino_sites(input_file_path):
    """ Translate the coding region of the reference sequence, 
    match the base and amino acid positions accordingly, and create a dictionary. """
    cds_region_file = f"{input_file_path}/cds_regions.txt"
    references_file = f"{input_file_path}/reference_seq.fasta"
    cds_meta = pd.read_csv(cds_region_file,
                             dtype="str",
                             keep_default_na=False,
                             sep="\t"
                             )
    # cds_region_title = cds_region.columns.tolist() # ['reference', 'cds_name', 'start_ends', 'gene/locus_tag']
    # print(cds_region_title)

    # 参考序列名作为键,序列作为值                                                     
    with open(references_file) as file:                   
        ref_seqs = {x.split("\n")[0]:"".join(x.split("\n")[1:]) for x in file.read().strip().split(">")[1:]}

    # Store the correspondence between the bases and amino acids in the coding region.
    nflu_component_site_list = []
    # Store the reference amino acid sequence.                                       
    ref_protein_dict = {}  
    # The corresponding amino acid for each base at every position on the reference sequence,
    #  where some positions have a single correspondence while others have multiple correspondences
    base_amino_sites_dict = {}                                             
    for m in cds_meta.values:
        m = list(m)
        # Reference sequence name in the metadata.
        ref = m[0]  
        if ref not in base_amino_sites_dict:
            base_amino_sites_dict[ref] = {}  
        # Protein name encoded by the CDS in the metadata.                                                  
        protein = m[1]
        # The start and end positions of the CDS region in the metadata.                                               
        start_ends = []                                                 
        # The processing of splicing together multiple coding regions to form a protein.
        if ";" in m[2]:                                                 
            for region in m[2].split(";"):
                start_ends.append(region.split("-"))
        else:                                                       
            start_ends.append(m[2].split("-"))
        # The gene name in the metadata.
        gene = m[3]  
        # To obtain the coding region sequence, 
        # which represents the complete nucleotide sequence encoding a protein, 
        # in Python indexing starts from 0, so the starting position needs to be reduced by 1.                                                   
        cds_base = ""                                                 
        for base_region in start_ends:
            cds_base = cds_base+ref_seqs[ref][int(base_region[0])-1:int(base_region[1])]
        # Translate the entire coding sequence of the protein into an amino acid sequence.
        protein_seq = Seq(cds_base).translate()                    
        ref_protein_dict[ref+"___"+gene+"___"+protein] = protein_seq
        # The protein amino acid sequence output here starts with 'M' and ends with '*'.
        # print(ref,protein,protein_seq)                           

        """ Only the amino acid sequence ending with a stop codon 
        and without '*' in other parts is retained. """
        if protein_seq.endswith("*"):
            if protein_seq.startswith("M"):
                """ Put bases, base positions, amino acids, 
                and amino acid positions into separate lists, 
                all of which have the same length. """
                # The list of base positions in the coding sequence (CDS) region,
                #  when the coding region is in segments, 
                # the base positions are not continuous
                ref_base_site_list = []                                                      
                for base_region in start_ends:
                    for site in range(int(base_region[0]),int(base_region[1])+1):
                        ref_base_site_list.append(site)
                # The list of bases in the coding sequence (CDS) region, 
                # convert the entire protein-coding base sequence into a list.
                ref_base_list = list(cds_base) 
                # The list of amino acid residue positions in the coding sequence (CDS) region, 
                # each position stored three times. This list has the same length as the base position list.         
                ref_gene_amino_site_list = []                                           
                for index in [x for x in range(1,len(protein_seq)+1)]:
                    for number in range(0,3):
                        ref_gene_amino_site_list.append(index)
                # The list of encoded amino acids in the coding sequence (CDS) region, 
                # with each amino acid stored three times. 
                # This list has the same length as the base position list
                ref_gene_amino_list = []                                                   
                for amino in  protein_seq:
                    for number in range(0,3):
                        ref_gene_amino_list.append(amino)

                for index,value in enumerate(ref_base_site_list):
                    result_row = [ref,gene,ref_base_site_list[index],
                                  ref_base_list[index],protein,
                                  ref_gene_amino_site_list[index],
                                  ref_gene_amino_list[index]]
                    nflu_component_site_list.append(result_row)

                """ Add the list of corresponding situations for each position 
                on the reference sequence to the 'base_amino_sites_dict' dictionary, 
                for the annotation of nucleotide variations. """
                # The indices of four identical lists with the same length: 
                # (base, base position, amino acid, amino acid position).
                index_start = 0    
                for index in range(0,len(ref_base_site_list)//3):
                    # "Retrieve the positions of the three bases in a codon.
                    ref_base_site = ref_base_site_list[index_start:index_start+3] 
                    # Retrieve the three bases of a codon.              
                    ref_base = ref_base_list[index_start:index_start+3]  
                    # Retrieve the positions of the amino acids corresponding to the three bases in a codon.                       
                    ref_gene_amino_site = ref_gene_amino_site_list[index_start:index_start+3]   
                    # Retrieve the amino acids corresponding to the three bases in a codon.
                    ref_gene_amino = ref_gene_amino_list[index_start:index_start+3]             
                    """ Retrieve the positions of the reference sequence in this encoded protein. """
                    for site in ref_base_site:
                        site_dict = {
                                    # The position of a base within a codon, 
                                    # for ease of annotation, is denoted as 0, 1, or 2.
                                    "index_three": 0,   
                                    # The three bases of a codon.                                       
                                    "three_base": "", 
                                    # The position of an amino acid encoded by a codon in a protein sequence.                                          
                                    "amino_site": 0,      
                                    # The amino acid encoded by a codon.                                      
                                    "amino": "",    
                                    # The corresponding coding region's name (protein).                                 
                                    "cds_name": "",  
                                    # The corresponding gene's name.                                           
                                    "gene": ""}                                                 
                        index_three = ref_base_site.index(site)
                        site_dict["index_three"] = index_three
                        site_dict["three_base"] = ref_base
                        site_dict["amino_site"] = ref_gene_amino_site[index_three]
                        site_dict["amino"] = ref_gene_amino[index_three]
                        site_dict["cds_name"] = protein
                        site_dict["gene"] = gene
                        if protein not in base_amino_sites_dict[ref]:
                            base_amino_sites_dict[ref][protein] = {}
                        # Using the reference sequence name followed by
                        #  the base position as the key in the dictionary, linked with "___".
                        key_refernce_site = site      
                        if key_refernce_site not in base_amino_sites_dict:                      
                            base_amino_sites_dict[ref][protein][key_refernce_site] = site_dict                                                             
                    index_start = index_start + 3

    title = ["reference","gene","ref_site","ref_nucleicacid","protein","ref_aminoacid_site","ref_aminoacid"]
    amino_aite_df = pd.DataFrame(nflu_component_site_list,columns=title)

    return amino_aite_df,ref_protein_dict,base_amino_sites_dict


def non_cds_site(gene_site_base_df,amino_site_df):
    """ Generate a dictionary for non coding regions to be used for SNP annotation. """
    gene_site_dict = {}
    for row in gene_site_base_df.values:
        row = list(row)
        key = tuple([row[0],row[2],row[1]])
        gene_site_dict[key] = row

    amino_base_site_dict = {}
    for row in amino_site_df.values:
        row = list(row)
        key = tuple([row[0],row[2],row[1]])
        amino_base_site_dict[key] = row
    
    non_cds_site_row = []
    for ref_site in gene_site_dict:
        if ref_site not in amino_base_site_dict:
            row = gene_site_dict[ref_site]
            non_cds_site_row.append(row)
    result = non_cds_site_row
    title = ["reference","gene","ref_site","ref_nucleicacid"]
    df = pd.DataFrame(result,columns=title)
    
    non_cds_site_dict = {}
    for row in df.values:
        row = list(row)
        reference = row[0]
        ref_site = row[2]
        if reference not in non_cds_site_dict:
            non_cds_site_dict[reference] = {}
            if ref_site not in non_cds_site_dict[reference]:
                non_cds_site_dict[reference][ref_site] = [row]
            else:
                non_cds_site_dict[reference][ref_site].append(row)
        else:
            if ref_site not in non_cds_site_dict[reference]:
                non_cds_site_dict[reference][ref_site] = [row]
            else:
                non_cds_site_dict[reference][ref_site].append(row)

    return non_cds_site_dict


def split_query_add_reference(input_file_path,praperation_directory_path,isolate_name,work_dir):
    """ Split the query sequence file, 
    with each query sequence followed by an associated reference sequence. 
    The query sequence comes first, followed by the reference sequence. """
    input_querys_path = f"{work_dir}/h9_seg8_clade/script/call_snp_and_annotation/work_dir/{isolate_name}/input"
    querys_files = f"{input_querys_path}/querys"
    query_seqs_path = f"{praperation_directory_path}/query_seqs"
    references_path = f"{praperation_directory_path}/references"

    # Create a split query sequence file.
    if os.path.exists(query_seqs_path):
        shutil.rmtree(query_seqs_path)
    os.mkdir(query_seqs_path)

    # Open the list of reference sequence file names.
    with open(f"{praperation_directory_path}/reference_list.txt") as file:
        reference_content = file.read().strip().split("\n")
        reference_list = list(set(reference_content))
    if len(reference_list) < len(reference_content):
        print(""" The entered list of reference sequence names contains duplicates.
               Please remove the duplicates and input again! """)

    split_query_number = 1
    if len(reference_list) == len(reference_content):
        for i in reference_list:
            if os.path.exists(f"{querys_files}/{i}___querys.fasta"):
                # Read the reference sequences.
                with open(f"{references_path}/{i}.fasta") as file:
                    reference = file.read().strip()
                # Read the query sequences.
                with open(f"{querys_files}/{i}___querys.fasta") as file:
                    querys = file.read().strip().split(">")[1:]
                # Write the query sequences and reference sequences into a separate file.
                for query in querys:
                    query_list = query.split("\n")
                    query_name = query_list[0]
                    query_seq = []
                    for x in query_list[1:]:
                        # If there are any characters other than 
                        # ATCGU in the DNA sequence, replace them with '-'.
                        x = re.sub(r"[^ATCGU]","U",x.upper())             						
                        query_seq.append(x)
                    query_seq = "\n".join(query_seq)
                    query_name_seq = ">"+query_name+"\n"+query_seq
                    query_reference = query_name_seq.strip()+"\n"+reference.strip()
                    if query_reference.strip() == "":
                        print(query_name)
                    if ">" in query_reference:
                        # "Use 'split_query_number' as the filename for the split file,
                        #  instead of the sequence name."
                        with open(f"{query_seqs_path}/{str(split_query_number)}.fasta",mode="a") as file:
                            file.write(query_reference)
                        split_query_number = split_query_number+1


if __name__ == "__main__":
    params = parse_args()
    input_file_path = params.input_file_path
    work_path = params.work_path
    isolate_name = params.isolate_name
    work_dir = params.work_dir
    # Create a folder to store praperation files.
    praperation_directory_path = praperation_directory(work_path)

    # Split the reference sequence file.
    split_reference_files(input_file_path,praperation_directory_path)

    # Label the positions of the reference sequence with gene labels.
    gene_site_base_df = gene_site_base(input_file_path)

    # Generate a dictionary for coding regions to be used for SNP annotation.
    base_amino_sites_list = base_amino_sites(input_file_path)
    amino_site_df = base_amino_sites_list[0]
    ref_protein_dict = base_amino_sites_list[1]
    base_amino_sites_dict = base_amino_sites_list[2]
    base_amino_sites_path = f"{praperation_directory_path}/base_amino_site.json"
    with open(base_amino_sites_path,mode="a") as file:
        json.dump(base_amino_sites_dict,file,indent=4)

    # Generate a dictionary for non coding regions to be used for SNP annotation.
    non_cds_site_dict = non_cds_site(gene_site_base_df,amino_site_df)
    non_cds_site_json_path = f"{praperation_directory_path}/non_cds_site.json"
    with open(non_cds_site_json_path,mode="a") as file:
        print(1)
        json.dump(non_cds_site_dict,file,indent=4)

    split_query_add_reference(input_file_path,praperation_directory_path,isolate_name,work_dir)

    print("The preparation of files for SNP calling and annotation is complete!")