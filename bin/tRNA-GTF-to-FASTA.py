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

count = 0
for line in readlines:
    if line.startswith("#"):
        pass
    else:
        splitline = line.strip().split("\t")
        featureID = splitline[8].split('"')
        geneName = featureID[1]
        if featureID[1] in geneid_set:
            """ 
            If the gene ID has already been added to this set, 
            do not add it again and do not write it to the GTF 
            and FASTA files
            """
            continue#pass
        else:
            geneDescription = splitline[8]
            if "Mt_tRNA" in line:
                count = count + 1
                geneName = splitline[8].split("gene_name ")[1].split('"')[1]
                geneName = geneName.replace("-", "_")
                geneName = "chrMT.tRNA" + str(count) + "-" + geneName # Create new name for mitochondrial tRNAs
                #geneDescription = 'gene_id "' + geneName + '"' + " ".join(splitline[8].split('"', 2)[2:]) # Replace gene_id with new one
                geneDescription = 'gene_id "' + geneName + '"; gene_version "1"; gene_name "' + geneName + '"; gene_source "insdc"; gene_biotype "Mt_tRNA";' # Replace gene_id with new one
            geneid_set.add(featureID[1])
            try:
                chromosome = genomeSeqIOdict[splitline[0]].id
            except: # If chromosome name from GTF not in genome file, skip. 
                print("WARNING: Couldn't find Chromosome %s in genome file, skipping the following feature..." % splitline[0])
                print(splitline)
                continue 
            chromosomeLen = len(genomeSeqIOdict[splitline[0]].seq)
            left_coord = int(splitline[3]) - 50 # -50 to get 50 nt upsteam
            right_coord = int(splitline[4]) + 50 # +50 to get 50 nt downstream
            ### Extract sequence from chromosome:
            if (left_coord > 0) and (right_coord < chromosomeLen):
                sequence = genomeSeqIOdict[chromosome].seq[left_coord:right_coord]  #Defines chromosome and coordinates
            ### If the left coordinate falls out of bounds of chromosome, add NNNs to fill space
            if left_coord < 0:
                ambiguousNucleotideNumber = abs(left_coord) + 1 #  + 1 to correct for passing 0
                ambigNucString = ''.join([char*ambiguousNucleotideNumber for char in "N"]) # Replicate N to fill in neg space
                sequence = ambigNucString + genomeSeqIOdict[chromosome].seq[0:right_coord]  #Defines chromosome and coordinates
            ### If the right coordinate falls out of bounds of chromosomes, add NNNs to fill space
            if right_coord > chromosomeLen:
                ambiguousNucleotideNumber = (right_coord - chromosomeLen) # Get no. of NNNs we need
                ambigNucString = ''.join([char*ambiguousNucleotideNumber for char in "N"]) # Replicate N to fill in neg space
                sequence = genomeSeqIOdict[chromosome].seq[left_coord:chromosomeLen] + ambigNucString # Get sequence and add NNNs
            gtf = (geneName, splitline[1], splitline[2], str(1), str((right_coord + 1) - left_coord), splitline[5], splitline[6], splitline[7], geneDescription) # +1 to account for python indexing
            newgtf.write("\t".join(gtf) + "\n")
            if splitline[6] == "+":
                newfile.write(">" + geneName + "\t" + splitline[8] + " (+)\n" + str(sequence) + "\n")	
            elif splitline[6] == "-":
                revSequence = sequence.reverse_complement()
                newfile.write(">" + geneName + "\t" + splitline[8] + " (-)\n" + str(revSequence) + "\n")

newfile.close()
newgtf.close()
