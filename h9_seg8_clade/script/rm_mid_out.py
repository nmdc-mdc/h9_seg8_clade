import os 


cmd = "rm -rf /h9_seg8_clade/parse_input/parse_output"
os.system(cmd)

cmd = "rm -rf /h9_seg8_clade/blastn_filter/blast_output"
os.system(cmd)

f = "/h9_seg8_clade/filter"
cmd = f"rm -rf {f}/qc_filter_not_pass.json {f}/qc_filter_pass.json {f}/querys"
os.system(cmd)

c = "/h9_seg8_clade/call_snp_and_annotation/"
cmd = ("rm -rf " + c + "align_result " + c + 
		"annotation_result " + c + "annotation_result_table.csv " + c + 
		"differ_site " + c + "praperation " + c + "snp_reuslt_table.csv " + 
		"/h9_seg8_clade/call_snp_and_annotation/input/querys/*")
os.system(cmd)

cmd = "rm -rf /h9_seg8_clade/seq_vari_new/input_seq_vari_new_vari.csv"
os.system(cmd)

cmd = "rm -rf /h9_seg8_clade/output"
os.system(cmd)
