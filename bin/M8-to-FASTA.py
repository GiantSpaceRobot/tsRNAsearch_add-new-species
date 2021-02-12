#!/usr/bin/env python
 
"""
Extract sequences from genome using M8
and write new GTF with correct coordinates
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
from Bio import SeqIO
from Bio.Seq import Seq
import sys

gff = open(sys.argv[1], "r")
newfile = open(sys.argv[3], "w")
newgff = open(sys.argv[4], "w")
readlines = gff.readlines()
gff.close()

genomefile = open(sys.argv[2], "rU")
genomeSeqIOdict = SeqIO.to_dict(SeqIO.parse(genomefile, "fasta"))
genomefile.close()

geneid_set = set()

count = 0
for line in readlines:
    count = count + 1
    splitline = line.strip().split("\t")
    gene_name = splitline[0]
    if "chr" in gene_name:
        gene_name = gene_name.replace("chr", "c")
    if ".trna" in gene_name:
        gene_name = gene_name.replace(".trna", "t")
    if ".tRNA" in gene_name:
        gene_name = gene_name.replace(".tRNA", "t")
    if gene_name.startswith("ENS"):
        gene_name = gene_name.replace("ENS", "NCRNA")
    lookalike_name = "lookalike" + str(count) + "_" + gene_name
    attributes="gene_id \"" + str(lookalike_name) + "\"; transcript_id \"" + str(lookalike_name) + "\";" 
    chrom = splitline[1]
    start = splitline[8]
    stop = splitline[9]
    query_start = splitline[6]
    query_stop = splitline[7]
    if int(query_start) < int(query_stop):
        # Extract + strand
        strand = "+"
    else:
        # Extract minus strand
        strand = "-"
    gff = (chrom, "mmseqs", "exon", start, stop, ".", strand, ".", attributes)
    newgff.write("\t".join(gff) + "\n")
    sequence = genomeSeqIOdict[chrom].seq[int(start):int(stop)]  #Defines chromosome and coordinates
    if strand == "+":
        newfile.write(">" + lookalike_name + "\n" + str(sequence) + "\n")	
    elif strand == "-":
        revSequence = sequence.reverse_complement()
        newfile.write(">" + lookalike_name + "\n" + str(revSequence) + "\n")
newfile.close()
newgff.close()
