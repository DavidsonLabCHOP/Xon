#!/usr/bin/env python

import sys
import os
import argparse
import itertools
import re

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Input a intropolis dataset.')
parser.add_argument('-l', '--lineNumber', required=True, help='Input the line number you wish to print.')
args = parser.parse_args()

num = int(args.lineNumber)

with open(args.input, "r") as fp:
    for i, line in enumerate(fp):
        if i == num:
            # line of interest.
            split1 = re.split(r'\t+', str(line.rstrip()))
            chromosome = split1[0]
            startPos = split1[1]
            endPos = split1[2]
            countsList = re.split(',', split1[7])
            for i in range(0, len(countsList)): 
                countsList[i] = int(countsList[i]) 
            sumCounts = sum(countsList)
            numDatasets = len(countsList)
            print(str(num), ",", chromosome, ":", startPos, "-", endPos, ",", numDatasets, ",", sumCounts, sep='')
            #print(str(line))
        elif i > num:
            break