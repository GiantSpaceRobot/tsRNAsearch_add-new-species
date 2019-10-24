#!/usr/bin/env python
 
"""
Combine tRNAs and remove introns to create new GTF
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys

gff = open(sys.argv[1], "r")   # GTF
newgff = open(sys.argv[2], "w") # New GTF
readlines = gff.readlines()
gff.close()

count = 0
for line in readlines:
    count = count + 1
    splitline = line.strip().split(" ")
    if (count % 2) == 0:   # Count is even (we are on a even number file line)
        newgff.write(splitline[3] + "\n")
    else: # Count is odd, we are on an odd line of the file
        newgff.write(splitline[0] + "\t" + splitline[4] + "\t")

newgff.close()
