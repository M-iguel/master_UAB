import sys #import packages
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Process fastx files.') #call function argparse
parser.add_argument('--input', type=str, help="Introduce an output file in fastx format.",required=True) #add arguments and specify type and help instructions (--help)
parser.add_argument('--output', type=str, help="Introduce an output file in fastx format.",required=True)
parser.add_argument('--operation', type=str, choices=['rc', 'trim', 'remove-adaptor'], help="Choose an operation to perform.",required=True)
parser.add_argument('--trim-left', type=int, help="Number of bases trimmed from the left", required=False) #some arguments are not compulsory
parser.add_argument('--trim-right', type=int, help="Number of bases trimmed from the right", required=False)
parser.add_argument('--adaptor', type=str, help="Introduce an adaptor to remove from the sequences.", required=False)
args = parser.parse_args() #save arguments into objects with attributes
variables = vars(args) #store values in a dictionary
if variables["input"].split(".")[::-1][0] == variables["output"].split(".")[::-1][0]: #if input and output files share extension
    print("Operation: "+variables["operation"]) #print operation
else: #if not
    print("[File Extension Error]: Please make sure your input and output files have the same extension (fastq/fasta).")
    sys.exit() #print Error and exit program

def remove_adap(seq, adap): #define function to remove adaptors
    adap_no = 0
    yes=""
    if seq.startswith(adap):  #if sequence starst with adaptor
        seq2 = seq[len(adap):]  #adaptor is removed
        adap_no += 1 #keep track of how many adaptors we delete
        yes="y" # save label to know if the adaptor was removed
    else: #if there is no adaptor the sequence remains the same
        seq2 = seq
        yes="n"
    return seq2, adap_no,yes

def trim(sec, x, y): #define trim function given the sequence and the number of bases to delete from each side (x and y)
    seq1 = sec[x:len(sec) - y]  #we retrieve the sequence contained between x and len(seq)-y
    return seq1

def rev_comp(sequence): #define function to write the reverse complement with a dictionary
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N"
    }
    return "".join(complement[sequence[i]] for i in range(len(sequence) - 1, -1, -1))

def summary(seqq): #define function to compute % of each base
    base_proces = len(seqq) #total number of bases processed
    a_per = (seqq.count("A") / len(seqq)) * 100 #individual base count
    c_per = (seqq.count("C") / len(seqq)) * 100
    g_per = (seqq.count("G") / len(seqq)) * 100
    t_per = (seqq.count("T") / len(seqq)) * 100
    n_per = (seqq.count("N") / len(seqq)) * 100
    result = str(base_proces) + " bases processed " + "(" + str(round(a_per, 2)) + "% A, " + str(
        round(c_per, 2)) + "% C, " + str(round(g_per, 2)) + "% G, " + str(round(t_per, 2)) + "% T, " + str(
        round(n_per, 2)) + "% N)"
    return result #return the string containing information

def summary2(archivo): #define function to compute the summary of the whole file
    string=""
    for record in SeqIO.parse(archivo, archivo.split(".")[::-1][0]): #extract all sequences with SeqIO library
        string=string+str(record.seq) #append every sequence in a single string
    return summary(string) #call function summary function over that string

def write_seq(line,operacion): #define operation to write the modified DNA sequence depending on the chosen operation
    trimmed = ""
    if operacion == "rc": #reverse complement
        sequ=rev_comp(line.strip("\n"))  #line will be the sequence line extracted from the input file
    elif operacion == "trim":
        trimmed = trimmed + line[:variables["trim_left"]] + line[(len(line) - variables["trim_right"]):] #store all trimmed sequences in a single string
        sequ=trim(line.strip("\n"), variables["trim_left"], variables["trim_right"])
    elif operacion == "remove-adaptor":
        sequ=remove_adap(line.strip("\n"), variables["adaptor"])[0] #argument 0 of remove_adap function is the modified sequence
    return sequ,trimmed

def write_qual(line,operacion,adap_no,yes): #define function to write qualities in fastq formats
    qual=""
    if operacion == "rc":
        qual=line.strip("\n")[::-1] #reverse quality line in case of reverse-complementing the DNA sequence
    elif operacion == "trim":
        qual=trim(line.strip("\n"), variables["trim_left"], variables["trim_right"]) #also trim the quality sequence
    elif operacion == "remove-adaptor":
        if yes=="y": #if the adaptor was removed from the begining
            qual = line.strip("\n")[len(variables["adaptor"]):]  # delete from the begining len(adaptor) positions
        else:
            qual=line.strip("\n") #dont modify the sequence if there is no adaptor
    return qual

def imprimir(line_no,input_file,trimmed,adap_no,operacion): #define imprimir function to display the summaries
    print(str(line_no // 4) + " reads processed.")
    print(summary2(input_file)+".")
    if operacion == "trim":
        print(summary(trimmed) + " and trimmed.") #compute the summary of the trimmed sequences
    if operacion == "remove-adaptor":
        print(str(adap_no) + " adaptor(s) removed") #print number of adaptors removed

def write_out(operacion, input_file, output_file): #define function to write the output file
    trimmed=""
    yes=""
    adap_no=0
    f = open(input_file, "rt")  # open file
    out = open(output_file, "wt")  # create output file
    line_no = 0
    if input_file.split(".")[::-1][0] == "fastq": #if the input is fastq
        for line in f: #for every line in the file
            if line_no % 4 == 0:  # first line is tag
                out.write(line)  # write the tag in the new file
            elif line_no % 4 == 1: # second line is the sequence
                out.write(write_seq(line,operacion)[0]+"\n") #call function that returns the modified sequence (line) as its first argument ([0])
                if operacion=="trim": # if operation was trim
                    trimmed=trimmed+write_seq(line,operacion)[1] #save trimmed sequences (second argument of the function trim ([1]))
                elif operacion=="remove-adaptor":
                    adap_no += remove_adap(line.strip("\n"), variables["adaptor"])[1] # save number of adaptors (second argument of function remove_adap)
                    yes=remove_adap(line.strip("\n"), variables["adaptor"])[2]
            elif line_no % 4 == 2: # third line is the "+"
                out.write(line.strip("\n")+"\n")  # write the separator
            elif line_no % 4 == 3:  # fourth line is the quality
                out.write(write_qual(line,operacion,adap_no,yes)+"\n") #call the quality-writing function
            line_no += 1 # continue iterating
        imprimir(line_no,input_file,trimmed,adap_no,operacion) # display summaries in the terminal
    elif input_file.split(".")[::-1][0] == "fasta": # if the input file is in fasta format
        for line in f:
            if line.startswith(">"):
                out.write(line) # keep the tag
            else:
                out.write(write_seq(line,operacion)[0]+"\n") # write modified sequence
                if operacion=="trim":
                    trimmed=trimmed+write_seq(line,operacion)[1]
                elif operacion=="remove-adaptor":
                    adap_no += remove_adap(line.strip("\n"), variables["adaptor"])[1]
            line_no += 1
        imprimir(line_no*2,input_file,trimmed,adap_no,operacion) # display summaries -> note line_no*2 because fasta files do not have separator and qualities
    f.close() # close files
    out.close()

write_out(variables["operation"], variables["input"], variables["output"]) # call function providing operation, input file and output file from the arguments dictionary
