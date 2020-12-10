#!/usr/bin/env python
 
"""
Extract sequences from genome using GTF
and write new GTF with correct coordinates
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
from Bio import SeqIO
from Bio.Seq import Seq
import sys

gtf = open(sys.argv[1], "r")
newfile = open(sys.argv[3], "w")
newgtf = open(sys.argv[4], "w")
readlines = gtf.readlines()
gtf.close()

genomefile = open(sys.argv[2], "rU")
genomeSeqIOdict = SeqIO.to_dict(SeqIO.parse(genomefile, "fasta"))
genomefile.close()

geneid_set = set()

for line in readlines:
    if line.startswith("#"):
        pass
    else:
        splitline = line.strip().split("\t")
        featureID = splitline[8].split('"')
        if featureID[1] in geneid_set:
            """ 
            If the gene ID has already been added to this set, 
            do not add it again and do not write it to the GTF 
            and FASTA files
            """
            pass
        else:
            geneid_set.add(featureID[1])
            left_coord = int(splitline[3]) - 50 # -50 to get 50 nt upsteam
            right_coord = int(splitline[4]) + 50 # +50 to get 50 nt downstream
            gtf = (featureID[1], splitline[1], splitline[2], str(1), str((right_coord + 1) - left_coord), splitline[5], splitline[6], splitline[7], splitline[8]) # +1 to account for python indexing
            newgtf.write("\t".join(gtf) + "\n")
            sequence = genomeSeqIOdict[splitline[0]].seq[left_coord:right_coord]  #Defines chromosome and coordinates
            if splitline[6] == "+":
                newfile.write(">" + featureID[1] + "\t" + splitline[8] + " (+)\n" + str(sequence) + "\n")	
            elif splitline[6] == "-":
                revSequence = sequence.reverse_complement()
                newfile.write(">" + featureID[1] + "\t" + splitline[8] + " (-)\n" + str(revSequence) + "\n")

newfile.close()
newgtf.close()
