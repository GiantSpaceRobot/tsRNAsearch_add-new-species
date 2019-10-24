#!/usr/bin/env python
 
"""
Write tRNA lengths to file
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys

gff = open(sys.argv[1], "r")
readlines = gff.readlines()
gff.close()

my_list = list()
my_dict = dict()

for line in readlines:
    splitline = line.strip().split("\t")
    feature_end = int(splitline[4]) 
    feature_id = (splitline[0].split('-')[1])[:3]
    if feature_id not in my_dict: # Add all tRNA lengths to a dictionary (e.g. Gly: 78,75,76)
        my_dict[feature_id] = [feature_end]
    else:
        my_dict[feature_id].append(feature_end)

for k,v in my_dict.iteritems():
    keyval = (k + "\t" + str(max(v)) + "\n")
    my_list.append(keyval)
my_list.sort()

with open(sys.argv[2], 'w') as f:
    for item in my_list:
        f.write(item)
