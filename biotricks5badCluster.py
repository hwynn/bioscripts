#!/usr/bin/python3

#we tried to make a cluster with this file. But the library is a machine learning one, and we don't want that.

import os
import math, string, sys
from Bio import AlignIO
#alignment = AlignIO.read("TIGR00056.sto", "stockholm")
#print(alignment)

from Bio.Seq import Seq #Seq()
from Bio.Alphabet import IUPAC #IUPAC.protein
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


my_prot1 = Seq("GMVLGLQGYVVL", IUPAC.protein)
my_prot2 = Seq("GMVLGLQGYLVL", IUPAC.protein)
my_prot3 = Seq("GMVFTIQVAREF", IUPAC.protein)
my_prot4 = Seq("GMIFTIQVAREF", IUPAC.protein)
my_prot5 = Seq("GGVIALQTYSTF", IUPAC.protein)
my_prot6 = Seq("GVVIAYQSAVQL", IUPAC.protein)
my_prot7 = Seq("GFAVALQGALQL", IUPAC.protein)

prots = [my_prot1, my_prot2, my_prot3, my_prot4, my_prot5, my_prot6, my_prot7]


#from Bio.SubsMat import MatrixInfo as matlist 

#matrix = matlist.blosum62 
#gap_open = -11
#gap_extend = -1

#scoreLists = [];
#for i in range(1,len(prots)):
	#sl = [];
	#for j in range(0,i):
		#for a in pairwise2.align.globaldx(prots[j], prots[i], matrix):
			#al1,al2, score, begin, end = a
		#sl.append(score);
	#scoreLists.append(sl);

#for item in scoreLists:
	#print(item);

import sys, getopt
from Bio import SeqIO,pairwise2
import Bio.SubsMat.MatrixInfo as matrices
import sklearn.cluster as cluster

scores = [[0 for i in range(len(prots))] for j in range(len(prots))]
for i in range(0,len(prots)):
	for j in range(0,len(prots)):
			a = pairwise2.align.globalds(prots[i],prots[j],matrix,gap_open,gap_extend)
			(s1,s2,score,start,end) = a[0]
			scores[i][j] = score

for item in scores:
	print(item);
