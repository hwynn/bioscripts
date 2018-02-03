#!/usr/bin/python3
#this is to show how biopython can work with clustalw. 
#It turns out it doesn't handle clustalw very well at all. 
import os
import math, string, sys
from Bio.Align.Applications import ClustalwCommandline
print("Yes this is running");
help(ClustalwCommandline)
#cline = ClustalwCommandline("clustalw2", infile="TIGR00056.SEED") #does nothing
#print(cline); #also does nothing
