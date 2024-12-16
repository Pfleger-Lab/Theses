#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 08:56:31 2023

@author: Joshua Dietrich
"""
import sys
import pandas as pd
import os
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import csv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast import NCBIXML
import time



start = time.time()
#%%  Get some statistics from raw reads
# Import the fastq as a list of seq records
# only for the desired length
filepath = os.getcwd()
filename = sys.argv[1]

#%%  Get stats from all reads
all_reads = [s for s in SeqIO.parse(os.path.join(filepath,filename),'fastq')]
num_of_reads = len(all_reads)
read_lengths = [len(r.seq) for r in all_reads]

print(sum(read_lengths)/1000,'kb of total reads')
print(len(all_reads),'total reads generated')

#%%  Get histogram from raw reads
binwidth = 10
figure(figsize=(7, 3), dpi=80)
plt.hist(read_lengths,bins=range(0, max(read_lengths) + binwidth, binwidth))
plt.xlim(0,12000)
#plt.ylim(0,5000)
plt.ylabel('Number of reads')
plt.xlabel('Read length (bases)')
plt.title('All reads')
plt.savefig(filename[:-6]+'_all_reads.png', bbox_inches='tight')
# plt.show()

#%%  Check for empty backbone

# Make database for blasting against
make_db_path = 'makeblastdb'

bkb_db_input_file= "pJD071.fa"
bkb_cline = NcbimakeblastdbCommandline(cmd=make_db_path,
                                  dbtype = "nucl",
                                  input_file = bkb_db_input_file,
                                  out="bkb_blast_db")
bkb_cline()
blast_n_path = 'blastn'

# bkb is 4151 bp
min_len = 3951
max_len = 4351
bkb_reads = [s for s in all_reads if min_len<len(s.seq)<max_len]

bkb_result = r"JGI_bkb_BLAST_output.xml"
E_VALUE_THRESH = 1e-6
bkb_out_file = [] # list of results to write to final output file
bkb_read_count = 0

for idx,curr_read_rec in enumerate(bkb_reads):
    curr_read_length = len(curr_read_rec.seq) 
    
    SeqIO.write([curr_read_rec],'temp_query_fasta.fasta','fasta')
    
    bkb_cline = NcbiblastnCommandline(cmd = blast_n_path,
                          query='temp_query_fasta.fasta',
                          db= 'bkb_blast_db', 
                          strand='both',
                          evalue=E_VALUE_THRESH,
                          num_alignments = 3,
                          max_hsps = 1,
                          out=bkb_result, 
                              outfmt=5)
    bkb_cline()
    
    ### Parse BLAST results
    with open(bkb_result,'r') as result_handle:
        for record in NCBIXML.parse(result_handle):
            bkb_out_file.append("Query:\t"+str(record.query.split()[0])+"\tlength=\t"+str(curr_read_length))
            if record.alignments: #skip queries with no match
                curr_bitscores = []
                curr_hits = []
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            hit_ID = align.hit_def.split()[0]
                            bkb_out_file.append('\t'+str(hit_ID))
                            bkb_out_file.append("\t\tAlignment length:\t"+str(hsp.align_length))
                            curr_bitscores.append(hsp.bits)
                            curr_hits.append(hit_ID)
                            bkb_out_file.append ("\t\tBit score:\t"+str(hsp.bits))
                    if len(curr_bitscores)>1:
                        scores_copy = curr_bitscores.copy()
                        scores_copy.remove(max(curr_bitscores))
                        score_diff = max(curr_bitscores)-max(scores_copy)
                        best_hit = curr_hits[curr_bitscores.index(max(curr_bitscores))]
                    else:
                        score_diff =501 # if only one hit found just use it
                        bkb_read_count += 1
            else:
                bkb_out_file.append("\tNo alignments found.")
    
### Save the BLAST output text files
out_file_name = 'bkb_'+filename[:-6]+'_BLAST_results.txt'
with open(out_file_name, 'w') as text_file:
    for line in bkb_out_file:
        text_file.write("%s\n" % line)
print(bkb_read_count,'reads aligned to the backbone')    
    

#%%  Get stats from just reads in range
min_len = 6000
max_len = 12000
reads = [s for s in all_reads if min_len<len(s.seq)<max_len]
num_of_reads = len(reads)
read_lengths = [len(r.seq) for r in reads]

print(sum([len(r.seq) for r in reads])/1000,'kb of reads between',min_len,'and',max_len)		 
print(len(reads),'reads generated between',min_len,'and',max_len)

#%%  Get histogram from raw reads in range
binwidth = 10
figure(figsize=(7, 3), dpi=80)
plt.hist(read_lengths,bins=range(0, max(read_lengths) + binwidth, binwidth))
plt.xlim(min_len,max_len)
plt.ylim(0)
plt.ylabel('Number of reads')
plt.xlabel('Read length (bases)')
plt.title(str(min_len)+'-'+str(max_len))
plt.savefig(filename[:-6]+'_'+str(min_len)+'-'+str(max_len)+'.png', bbox_inches='tight')
# plt.show()

#%% Make databases for blasting against
# make_db_path = r"C:\Program Files\NCBI\blast-2.13.0+\bin\makeblastdb.exe"
make_db_path = 'makeblastdb'

ALS_db_input_file= "JGI_lib_2.0_ALS_TUs.fa"
ALS_cline = NcbimakeblastdbCommandline(cmd=make_db_path,
                                  dbtype = "nucl",
                                  input_file = ALS_db_input_file,
                                  out="JGI_lib_2.0_ALS_TUs")

KARI_db_input_file= "JGI_lib_2.0_KARI_TUs.fa"
KARI_cline = NcbimakeblastdbCommandline(cmd=make_db_path,
                                  dbtype = "nucl",
                                  input_file = KARI_db_input_file,
                                  out="JGI_lib_2.0_KARI_TUs")

DHAD_db_input_file= "JGI_lib_2.0_DHAD_TUs.fa"
DHAD_cline = NcbimakeblastdbCommandline(cmd=make_db_path,
                                  dbtype = "nucl",
                                  input_file = DHAD_db_input_file,
                                  out="JGI_lib_2.0_DHAD_TUs")
ALS_cline()
KARI_cline()
DHAD_cline()


#%% # Import the library of sequences to check against
lib_file = os.path.join(os.getcwd(),'JGI_lib_2.0_TUs.csv')
lib_db = pd.read_csv(lib_file,index_col=0)

#%% # Just ALS, KARI, and DHADs
ALS_list = [list(a) for a in zip(list(lib_db.index[:72]),list(lib_db['Sequence'][:72]))]
KARI_list = [list(a) for a in zip(list(lib_db.index[72:165]),list(lib_db['Sequence'][72:165]))]
DHAD_list = [list(a) for a in zip(list(lib_db.index[165:240]),list(lib_db['Sequence'][165:240]))]
TU_lists = [ALS_list,KARI_list,DHAD_list]
TU_names = ['ALS','KARI','DHAD']

ALS_lengths = [len(i[1]) for i in ALS_list]
KARI_lengths = [len(i[1]) for i in KARI_list]
DHAD_lengths = [len(i[1]) for i in DHAD_list]
lengths = [ALS_lengths,KARI_lengths,DHAD_lengths]
print('Smallest TU for:')
for t,l in zip(TU_names,[min(i) for i in lengths]):
    print(t,l)
    
#%% # BLAST and find full cassettes
ALS_output = [] # list of results to write to final output file
KARI_output = []
DHAD_output = []
out_files = [ALS_output,KARI_output,DHAD_output] 
E_VALUE_THRESH = 1e-6
ALS_result = r"JGI_TU_2.0_ALS_BLAST_output.xml"
KARI_result = r"JGI_TU_2.0_KARI_BLAST_output.xml"
DHAD_result = r"JGI_TU_2.0_DHAD_BLAST_output.xml"

# blast_n_path = r"C:\Program Files\NCBI\blast-2.13.0+\bin\blastn.exe"
blast_n_path = 'blastn'

read_count_dict = dict()
read_dict = dict()
TU_count_dicts = []
curr_progress = -1

for TU_list in TU_lists:
    TU_count_dicts.append(dict.fromkeys([i[0] for i in TU_list],0))

for idx,curr_read_rec in enumerate(reads):
    curr_read_length = len(curr_read_rec.seq) 
    
    SeqIO.write([curr_read_rec],'temp_query_fasta.fasta','fasta')
    results = [ALS_result, KARI_result, DHAD_result] # list of BLAST output files
    
    ALS_cline = NcbiblastnCommandline(cmd = blast_n_path,
                          query='temp_query_fasta.fasta',
                          db= 'JGI_lib_2.0_ALS_TUs', 
                          strand='both',
                          evalue=E_VALUE_THRESH,
                          num_alignments = 3,
                          max_hsps = 1,
                          out=ALS_result, 
                              outfmt=5)

    KARI_cline = NcbiblastnCommandline(cmd = blast_n_path,
                          query='temp_query_fasta.fasta',
                          db= 'JGI_lib_2.0_KARI_TUs', 
                          strand='both',
                          evalue=E_VALUE_THRESH,
                          num_alignments = 3,
                          max_hsps = 1,
                          out=KARI_result, 
                          outfmt=5)     

    DHAD_cline = NcbiblastnCommandline(cmd = blast_n_path,
                          query='temp_query_fasta.fasta',
                          db= 'JGI_lib_2.0_DHAD_TUs', 
                          strand='both',
                          evalue=E_VALUE_THRESH,
                          num_alignments = 3,
                          max_hsps = 1,
                          out=DHAD_result, 
                          outfmt=5)
    ALS_cline()
    KARI_cline()
    DHAD_cline()
    
    ### Parse BLAST results
#         s_len = 100 # length to print
    # idx indicates ALS (0), KARI (1), or DHAD (2)
    cassette_match = [] # list of three strings with the ID of the best match for each TU
    for idx,(result,out_file) in enumerate(zip(results,out_files)):
        with open(result,'r') as result_handle:
            for record in NCBIXML.parse(result_handle):
                out_file.append("Query:\t"+str(record.query.split()[0])+"\tlength=\t"+str(curr_read_length))
                if record.alignments: #skip queries with no match
                    curr_bitscores = []
                    curr_hits = []
                    for align in record.alignments:
                        for hsp in align.hsps:
                            if hsp.expect < E_VALUE_THRESH:
                                hit_ID = align.hit_def.split()[0]
                                out_file.append('\t'+str(hit_ID))
                                out_file.append("\t\tAlignment length:\t"+str(hsp.align_length))
                                out_file.append('\t\t'+str(list(lib_db.loc[hit_ID][:5])))
                                curr_bitscores.append(hsp.bits)
                                curr_hits.append(hit_ID)
                                out_file.append ("\t\tBit score:\t"+str(hsp.bits))
                    if len(curr_bitscores)>1:
                        scores_copy = curr_bitscores.copy()
                        scores_copy.remove(max(curr_bitscores))
                        score_diff = max(curr_bitscores)-max(scores_copy)
                        best_hit = curr_hits[curr_bitscores.index(max(curr_bitscores))]
                    else:
                        score_diff =501 # if only one hit found just use it
                        best_hit = hit_ID
                    if max(curr_bitscores) > 1800 and score_diff > 500:
                        cassette_match.append(best_hit)
                        # TU_count_dicts[idx][best_hit] += 1 # always count TUs found
                    else:
                        cassette_match.append('NoUniqueMatch')
                else:
                    out_file.append("\tNo alignments found.")
                    cassette_match.append('NoUniqueMatch')
    cassette_ID = '_'.join(cassette_match)
    if 'NoUniqueMatch' not in cassette_ID:
        for idx,TU in enumerate(cassette_match):
            TU_count_dicts[idx][TU] += 1
        if cassette_ID in read_count_dict:
            read_count_dict[cassette_ID] += 1
        else:
            read_count_dict[cassette_ID] = 1
        
        read_dict[curr_read_rec.id] = cassette_match
            

#%% # Save the BLAST output text files
for result,out_file,TU in zip(results,out_files,TU_names):
    out_file_name = TU+'_'+filename[:-6]+'_full_cassette_BLAST_results.txt'
    with open(out_file_name, 'w') as text_file:
        for line in out_file:
            text_file.write("%s\n" % line)
            
#%% # Save the .csv with the unique cassette counts
with open(filename[:-6]+'_read_counts.csv', 'w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['Cassette','# of reads'])
    w.writerows(read_count_dict.items())
    
#%% # Save the .csv with the TU counts
for TU_count_dict,TU in zip(TU_count_dicts,['ALS','KARI','DHAD']):
    TUs_with_no_reads = [k for k,i in TU_count_dict.items() if i == 0]
    print(TU,"TUs with no reads:",TUs_with_no_reads)
    with open(filename[:-6]+TU+'_read_counts.csv', 'w',newline='') as f:
        w = csv.writer(f)
        w.writerow(['TU','# of reads'])
        w.writerows(TU_count_dict.items())
    
#%% # Save the .csv with the reads and their matches
with open(filename[:-6]+'_read_matches.csv', 'w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['Read ID','ALS','KARI','DHAD'])
    for key in sorted(read_dict.keys()):
        w.writerow([key] + read_dict[key])
		
end = time.time()
print('Time to run: ',(end - start)/60,'minutes')