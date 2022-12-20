#!/usr/bin/env python
# coding: utf-8

import sys
archivo=sys.argv[1]

# 1. Counts the abundance of each possible 3-mers (histogram)

#import packages
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import os
import re

def histogram(dicc):
    keys = dicc.keys()
    vals = dicc.values()
    plt.bar(range(len(keys)), vals,align='center', width=0.3)
    plt.title("3-mer distribution")
    plt.ylabel("3-mer count")
    plt.xticks(range(len(keys)), keys,rotation="vertical",fontsize=7)
    plt.savefig("3mer_histogram.png")#save histogram to png
    
def kmer_count(seq,k):
    dicc={}
    for i in range(len(seq)):
        kmer = seq[i:i+k]#store in kmer the positions DNA[0:2], DNA[1:3], DNA[2:4],...
        if len(kmer)==3:
            if kmer not in dicc:
                dicc[kmer]=0 #if kmer not in dicc -> create new key with value = 0
            dicc[kmer] += 1 #there are one key for each kmer and value is the count
    return dicc

string=""
for record in SeqIO.parse(archivo, archivo.split(".")[::-1][0]):#we extract the sequences from the file
    string=string+str(record.seq)#append all of them to a string
histogram(kmer_count(string,3))#compute 3-mer count and plot it


# 2. Splits the input in FASTA/FASTQ into files of only 1 sequence each

def splitfile(input_file):
    f=open(input_file,"rt")
    line_no=0
    counter=1
    if input_file.split(".")[-1]=="fastq":#for fastq files
        for line in f:
            if line_no%4==0:#open a new output file for each tag
                output_file=open(str(counter)+"_"+"".join(input_file.split(".")[0:-1])+".fastq","wt")
                counter+=1
                output_file.write(line.strip("\n")+"\n")#write tag
            elif line_no%4==1:
                output_file.write(line.strip("\n")+"\n")#write sequence
            elif line_no%4==2:
                output_file.write(line.strip("\n")+"\n")#write "+"
            elif line_no%4==3:#calidades
                output_file.write(line.strip("\n")+"\n")#write qualities
                output_file.close()
            line_no+=1
    elif input_file.split(".")[-1]=="fasta":#similar process for fasta files
        for line in f:
            if line.startswith(">"):
                output_file=open(str(counter)+"_"+"".join(input_file.split(".")[0:-1])+".fasta","wt")
                counter+=1
                output_file.write(line.strip("\n")+"\n")
            else:
                output_file.write(line.strip("\n")+"\n")
                output_file.close()
            line_no+=1
    f.close()
splitfile(archivo)


# 3. Hard-trims (from the right) all sequences from all files 20nt

cwd = os.getcwd()#get directory

# Get .txt files
for f_name in os.listdir(cwd):#iterate over every file in the directory
    match = re.search(r'(^\d+_)',f_name)#if it starts with numbers + _ it is a splitted file
    if match:
        output_file=open("".join(f_name.split(".")[0:-1])+"_trimmed"+"."+f_name.split(".")[1],"w")#open a new file for every splitted file
        for record in SeqIO.parse(f_name,f_name.split(".")[1]):#get the sequence
            cut_record = record[:-20]  # remove 20 nt from the write
            SeqIO.write(cut_record, output_file, f_name.split(".")[1])#re-write the record in the output file
        output_file.close()
        os.remove(f_name)#delete untrimmed files

# 4. Using the genome mapping tool BWA and the reference genome of the Scaromice Cerevisiae (any strain will do), aligns each of the files producing its corresponding SAM file

files = os.listdir('.') #store every file in directory in a list
os.system("bwa index s_cerevisiae.fna")#create index for bwa alignment
os.system("samtools faidx s_cerevisiae.fna")
for i in files:
    if "_trimmed" in i:#for every trimmed file
        os.system("bwa mem s_cerevisiae.fna %s > %s"%(i,str(i)+"_sam.sam"))#create the sam file

# 5. Merges all SAM files ignoring headers

files = os.listdir(".")
out_file_sam_merge=open("sam_merged.sam","wt")
for i in files:
    if i.endswith(".sam"):#for every .sam file in the directory
        l=open(i,"rt")
        for line in l:
            if not line.startswith("@"):#skip header (starts with @)
                out_file_sam_merge.write(line)#write remaining lines in the merged file
        l.close()
out_file_sam_merge.close()

# 6. Sorts the SAM file by chromosome and position

#lines = open("sam_merged.sam", 'r').readlines()
#output = open("sam_merged_sorted.sam", 'w')
#for line in sorted(lines, key=lambda line: line.split("\t")[2]):#sort files by 3 field when splitting by \t --> chromosome
 #   output.write(line)
#output.close()
os.system("cat sam_merged.sam | sort -k3,3 -k4n,4 | sed -r '/^\s*$/d' > sam_merged_sorted.sam")
# 7. Computes how many reads have been aligned

h=open("sam_merged_sorted.sam","rt")#the number of lines in the sorted file is the number of reads alligned
count=0
for line in h:
    count+=1
h.close()
print("#################\n"+str(count)+" reads aligned.\n#################")


