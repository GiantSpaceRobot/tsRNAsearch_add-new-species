###
# Remove GTF duplicates
###
import sys

### Sequence file
M = open(sys.argv[1], "r")
ReaderM = M.readlines()
M.close()

### GTF file
with open(sys.argv[2]) as file:
    lines = [line.strip() for line in file]

### GTF output
outputGTF = open(sys.argv[3], "w")

for eachline in ReaderM:
    StrippedM = eachline.strip()
    SplitM = StrippedM.split()
    SplitString1 = SplitM[0]
    if str(SplitString1).startswith(">"):
        tRNA = SplitString1[1:]
        for line in lines:
            if tRNA in line:
                #if tRNA == "chr9.tRNA4-ArgTCT":
                #    print "Yes, %s and %s" % (tRNA, line)
                newline = line.replace("\tgene\t", "\texon\t")
                outputGTF.write(newline + "\n")
                #print "Yes, %s and %s" % (tRNA, line)

outputGTF.close()
