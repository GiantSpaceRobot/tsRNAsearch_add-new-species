from Bio import SeqIO
import sys

if len(sys.argv) < 4:
    print("Usage: python FASTA-to-GTF.py GRCm38_tRNA GRCm38_tRNAs.fa Relative.fa Relative.gtf")
    exit(1)

feature_type = str(sys.argv[1])
fasta_sequences = SeqIO.parse(open(sys.argv[2]),'fasta')
out_file = open(sys.argv[3], "w")
gtf_file = open(sys.argv[4], "w")
for fasta in fasta_sequences:
    name, sequence, description = fasta.id, str(fasta.seq), str(fasta.description)
    namesplit = name.split("-")
    tRNAname = namesplit[1]
    codon = namesplit[2]
    tRNA = (description.split("(")[1]).split(")")[0] # Get tRNA name, e.g. 'tRNAscan SE chr13.trna34-AlaAGC'
    #tRNA = "chr" + tRNA.split("chr")[1] # Remove tRNAscan string if present
    tRNA = tRNA.split(" ")[-1] # Remove tRNAscan string if present, e.g. 'chr13.trna34-AlaAGC'
    tRNA = tRNA.split("-")[0] # Remove suffix, e.g. 'chr13.trna34'
    tRNA = '%s-%s%s' % (tRNA, tRNAname, codon) # Rebuild suffix, e.g. 'chr13.trna34-AlaAGC'   # Necessary for awkward/inconsistent annotations
    long_line = 'gene_id "%s"; transcript_id "%s";' % (tRNA, tRNA)
    #strand = (description.split("(")[2]).split(")")[0]
    #upper = 0
    lower = 0
    count = 0
    total = ""
    #exon1 = ""
    #intron = ""
    if "(+)" in description: # Determine strand that tRNA is on
        strand = "+"
        sequence = sequence + "CCA"
    else:
        strand = "-"
        sequence = sequence + "CCA" # GtRNAdb tRNA sequences are all in the same orientation
        #sequence = "ACC" + sequence
        
    #if strand == "+":
    #    sequence = sequence + "CCA"
    #elif strand == "-":
    #    sequence = "ACC" + sequence
    #else:
    #    print ("Error: ", tRNA, strand)
    
    if sum(1 for c in sequence if c.islower()) > 5: # If the sum total of lowercase in a sequence is greater than 5, check for intron
        for i in sequence:
            total = total + i
            if (i.isupper()):
                lower = 0
            elif (i.islower()):
                lower = lower + 1
            if int(lower) == 5:   # If there are 5 or more lowercase in a row, it is an intron
                exon1 = total[:-(lower)]
            elif int(lower) > 5:
                letterIndex = sequence[count]
                if (sequence[count + 1]).isupper():
                    exon2 = sequence[(count + 1):]
                    break
            count = count + 1
        gtf_file.write(tRNA + "\t" + feature_type + "\texon\t1" + "\t" + str(len(exon1)) + "\t1000" + "\t" + strand + "\t.\t" + long_line + "\n")
        gtf_file.write(tRNA + "\t" + feature_type + "\texon\t" + str(len(sequence) - len(exon2)) + "\t" + str(len(sequence)) + "\t1000" + "\t" + strand + "\t.\t" + long_line + "\n")
    else:
        gtf_file.write(tRNA + "\t" + feature_type + "\texon\t1" + "\t" + str(len(sequence)) + "\t1000" + "\t" + strand + "\t.\t" + long_line + "\n")
    
    out_file.write(">" + tRNA + "\t" + description + "\n" + sequence + "\n")

out_file.close()
gtf_file.close()


