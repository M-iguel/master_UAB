#!/usr/bin/env python
# coding: utf-8

import sys
pattern=sys.argv[1]
text=sys.argv[2]

#we get the blosum62 matrix from the Bio.Align package
from Bio.Align import substitution_matrices
blosum_62=substitution_matrices.load("BLOSUM62")#get blosum62 matrix


#following a dynamic programming algorithm we define the functions

def edit_distance_dp(pattern,text):
    for h in range(1,len(text)+1):
        for v in range(1,len(pattern)+1):
            pair=(pattern[v-1], text[h-1])#we get the tuple of the compared positions (e.g. (N,D))
            if pair not in blosum_62:#if that pair is not in blosum_62
                blosum_62[pair]=blosum_62[tuple(reversed(pair))]#we subsitute it by the reversed (if (N,D) is not in the matrix, we can use (D,N))
                #because blosum62 is a lower triangular matrix
    #init        
    dp_matrix = [[0 for _ in range(len(text)+1)] for _ in range(len(pattern)+1)]#define the distance matrix
    for v in range(len(pattern)+1):
        dp_matrix[v][0] = v*(-2)#assign values to the first column (multiply by -2 to set the penalty for deletion to -2)
    for h in range(len(text)+1):
        dp_matrix[0][h] = h*(-4)#assign values to the first row (multiply by -4 to set the penalty for insertion to -4)
    # Compute DP Matrix
    for h in range(1,len(text)+1):
        for v in range(1,len(pattern)+1):
            dp_matrix[v][h] = max (#penalizations are negative, thus we want the maximum
                dp_matrix[v-1][h-1] + blosum_62[(pattern[v-1], text[h-1])],#the penalty for subtitution depends on blosum62
                dp_matrix[v][h-1] -4,#penalty for insertion is -4
                dp_matrix[v-1][h] -2)#penalty for deletion is -2
    return dp_matrix

def backtrace_matrix(pattern,text,dp_matrix):#CIGAR
    v = len(pattern)
    h = len(text)
    cigar = []
    while v > 0 and h > 0:
        if dp_matrix[v][h] == dp_matrix[v-1][h] + 1:
            v -= 1
            cigar.insert(0,"D")
        elif dp_matrix[v][h] == dp_matrix[v][h-1] + 1:
            h -= 1
            cigar.insert(0,"I")
        else:
            v -= 1
            h -= 1
            if pattern[v] == text[h]:
                cigar.insert(0,"M")
            else:
                cigar.insert(0,"X")
    if v > 0:
        for _ in range(v): cigar.insert(0,"D")
    if h > 0:
        for _ in range(h): cigar.insert(0,"I")
    return cigar

def print_matrix(dp_matrix):
    for y in dp_matrix:
        for x in y:
            print("%3d" % (x),end='')
        print("\n",end='')

def pretty_print_alignment(pattern,text,cigar):#pretty print the alignment
    (pattern_txt,i) = ("",0)
    operation_txt = ""
    (text_txt,j) = ("",0)
    for op in cigar:
        if op == "M":
            pattern_txt += pattern[i]
            i += 1
            operation_txt += "|"
            text_txt += text[j]
            j += 1
        elif op == "X":
            pattern_txt += pattern[i]
            i += 1
            operation_txt += " "
            text_txt += text[j]
            j += 1
        elif op == "I":
            pattern_txt += " "
            operation_txt += " "
            text_txt += text[j]
            j += 1
        elif op == "D":
            pattern_txt += pattern[i]
            i += 1
            operation_txt += " "
            text_txt += " "
    print(pattern_txt)
    print(operation_txt)
    print(text_txt)

#call functions
print("Distance matrix:")
print_matrix(edit_distance_dp(pattern,text))
print("\nAlignment:")
pretty_print_alignment(pattern,text,backtrace_matrix(pattern,text,edit_distance_dp(pattern,text)))
