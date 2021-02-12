#!/usr/bin/env python
 
"""
Extract only the sequences from the 
genome that are NOT in the provided GTF
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"
 
from Bio import SeqIO
from Bio.Seq import Seq
import sys

tsv = open(sys.argv[1], "r")
readlines = tsv.readlines()
tsv.close()

gtf = open(sys.argv[2], "r")
readgtf = gtf.readlines()
gtf.close()

genomefile = open(sys.argv[3], "rU")
genomeSeqIOdict = SeqIO.to_dict(SeqIO.parse(genomefile, "fasta"))
genomefile.close()

#feature_dict = dict()
my_set = set()

### Make dictionary containing every nucleotide index that was identified in tRNA genome mmseqs blast
print("Building dictionary of nucleotide indices identified using mmseqs")
for line in readlines:
    linesplit = line.split("\t")
    chrom = linesplit[0]
    start = linesplit[1]
    stop = linesplit[2].strip()
    start_to_stop = chrom + "_" + start + "-"  + stop
    my_set.add(start_to_stop)

### Add GTF coordinates to the feature dict
print("Adding to dictionary using nucleotide indices from GTF")
for line in readgtf:
    linesplit = line.split("\t")
    chrom = linesplit[0]
    start = linesplit[3]
    stop = linesplit[4]
    start_to_stop = chrom + "_" + start + "-"  + stop
    my_set.add(start_to_stop)

### Function to replace nucleotide given sequence and nucleotide index
def replace_nucleotide_with_N(my_sequence, my_start_index, my_stop_index):
    my_start_index = int(my_start_index)
    my_stop_index = int(my_stop_index)
    string_len = my_stop_index - my_start_index + 1
    N_string = ''.join([char*string_len for char in "N"]) # Repeat N n times
    modified_sequence = my_sequence[:my_start_index - 1] + N_string + my_sequence[my_stop_index:] # -1 correcting for pythonic counting
    return modified_sequence

for coordinates in my_set:
    linesplit = coordinates.split("_")
    chrom = linesplit[0]
    start = linesplit[1].split("-")[0]
    stop = linesplit[1].split("-")[1]
    #print(chrom, start, stop)
    print("Replacing nucleotides in sequence %s: %s:%s" % (chrom, start, stop))
    sequence = genomeSeqIOdict[chrom].seq  #Defines chromosome and coordinates
    #print(sequence)
    ### Loop over indices and replace with Ns
    sequence = replace_nucleotide_with_N(sequence, start, stop)
    #print(sequence)
    genomeSeqIOdict[chrom].seq = sequence # Replace dict sequence with new masked one
    
print("Writing SeqIO dict to file...")
with open(sys.argv[4], "w") as handle:
    SeqIO.write(genomeSeqIOdict.values(), handle, 'fasta')
print("Done")
