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

feature_dict = dict()

### Make dictionary containing every nucleotide index that was identified in tRNA genome mmseqs blast
for line in readlines:
    linesplit = line.split("\t")
    chrom = linesplit[0]
    start = linesplit[1]
    stop = linesplit[2]
    start_to_stop = range(int(start), int(stop) + 1)
    if chrom in feature_dict.keys():
        dict_values = feature_dict[chrom]
        new_values = list(set(dict_values + start_to_stop))
        feature_dict[chrom] = (new_values)
    else:
        feature_dict[chrom] = (start_to_stop)

### Add GTF coordinates to the feature dict
for line in readgtf:
    linesplit = line.split("\t")
    chrom = linesplit[0]
    start = linesplit[3]
    stop = linesplit[4]
    start_to_stop = range(int(start), int(stop) + 1)
    if chrom in feature_dict.keys():
        dict_values = feature_dict[chrom]
        new_values = list(set(dict_values + start_to_stop))
        feature_dict[chrom] = (new_values)
    else:
        feature_dict[chrom] = (start_to_stop)

### Function to replace nucleotide given sequence and nucleotide index
def replace_nucleotide_with_N(my_sequence, my_index):
    modified_sequence = my_sequence[:my_index - 1] + "N" + my_sequence[my_index:] # -1 correcting for pythonic counting
    return modified_sequence

for k,v in feature_dict.iteritems():
    print("Replacing nucleotides in sequence %s" % (k))
    sequence = genomeSeqIOdict[k].seq  #Defines chromosome and coordinates
    #print(sequence)
    ### Loop over indices and replace with Ns
    for index in v:
        sequence = replace_nucleotide_with_N(sequence, index)
        #print(sequence)
    #print(sequence)
    genomeSeqIOdict[k].seq = sequence # Replace dict sequence with new masked one
    
print("Writing SeqIO dict to file...")
with open(sys.argv[4], "w") as handle:
    SeqIO.write(genomeSeqIOdict.values(), handle, 'fasta')
print("Done")
