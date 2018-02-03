#!/usr/bin/python3

#this file shows us how to run pairwise2 alignments, display them, and return their scores

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


for a in pairwise2.align.localms("ACCGT", "ACGC",1,-1,-2,-2):
      al1,al2, score, begin, end = a
print(score);


for a in pairwise2.align.localms(my_prot1, my_prot2,1,-1,-2,-2):
      al1,al2, score, begin, end = a
print(score);

from Bio.SubsMat import MatrixInfo as matlist 

matrix = matlist.blosum62 
gap_open = -11
gap_extend = -1

for a in pairwise2.align.globaldx(my_prot1, my_prot2, matrix):
	
	print(format_alignment(*a)) 

for a in pairwise2.align.globaldx(my_prot1, my_prot2, matrix):
	al1,al2, score, begin, end = a	
print(score);
