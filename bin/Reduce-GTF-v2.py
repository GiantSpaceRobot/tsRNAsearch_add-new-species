###
# Remove GTF duplicates
###
import sys

### Sequence file
M = open(sys.argv[1], "r")
ReaderM = M.readlines()
M.close()
my_set = set()

for line in ReaderM:
    if ">" in line:
        my_line = line.strip(">").split()[0]
        my_set.add(my_line)
print("Genes in FASTA input: %s" % len(my_set))

### GTF output
outputGTF = open(sys.argv[3], "w")

### GTF input
gtf = open(sys.argv[2], "r")
gtfread = gtf.readlines()
gtf.close()
print("Lines in GTF to start with: %s" % len(gtfread))

### Create new GTF by matching genes in set
for i in gtfread:
    gene = i.split('"')[1] 
    if gene in my_set:
        outputGTF.write(i)

outputGTF.close()
