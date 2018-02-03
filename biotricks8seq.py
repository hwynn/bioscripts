#!/usr/bin/python3

#This file is to show us how to create MultipleSeqAlignment objects and what we can do with them. 
#the goal is to create an output file that hmmbuild can use. 

import os
import math, string, sys
from Bio import AlignIO
#alignment = AlignIO.read("TIGR00056.sto", "stockholm")
#print(alignment)

from Bio.Seq import Seq #Seq()
from Bio.Alphabet import IUPAC #IUPAC.protein
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from Bio.Align import MultipleSeqAlignment

itemList = [];

for alignment2 in AlignIO.parse("TIGR00056_Mot(copy).sto", "stockholm"):
	print("Alignment of length %i" % alignment2.get_alignment_length())
	print("type of: ", type(alignment2));
	itemList.append(alignment2);
print("we have ", len(itemList), " alignments");

for whatever in itemList[0]:
	#print("type of whatever's inside alignment2: ", type(whatever));#<class 'Bio.SeqRecord.SeqRecord'>
	#print(whatever.seq, "   ", type(whatever.seq)); #GMVFTIQVAREF     <class 'Bio.Seq.Seq'>
	print(whatever.seq, whatever.id, whatever.features);
	#we don't get many features from this. maybe because it's just from our badly made stockholm file
	#but I'd love to be able to slice a sequence and have its location attributes adjusted
	#http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:SeqRecord-slicing
print(type(AlignIO.parse("TIGR00056_Mot(copy).sto", "stockholm"))); #<class 'generator'>


#http://biopython.org/DIST/docs/api/Bio.Align.MultipleSeqAlignment-class.html#add_sequence
#http://biopython.org/DIST/docs/api/Bio.AlignIO-module.html
##http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:SeqRecord-slicing

#my_prot1 = Seq("GMVLGLQGYVVL", IUPAC.protein)
#my_prot2 = Seq("GMVLGLQGYLVL", IUPAC.protein)
#my_prot3 = Seq("GMVFTIQVAREF", IUPAC.protein)
#my_prot4 = Seq("GMIFTIQVAREF", IUPAC.protein)
#my_prot5 = Seq("GGVIALQTYSTF", IUPAC.protein)
#my_prot6 = Seq("GVVIAYQSAVQL", IUPAC.protein)
#my_prot7 = Seq("GFAVALQGALQL", IUPAC.protein)

##my_prot1 = Seq("AAAAAAAAAAAA", IUPAC.protein)
##my_prot2 = Seq("AAAAAAAAAAAA", IUPAC.protein)
##my_prot3 = Seq("CCCCCCCCCCCC", IUPAC.protein)
##my_prot4 = Seq("AAAAAAAAAAAA", IUPAC.protein)
##my_prot5 = Seq("CCCCCCCCCCCC", IUPAC.protein)
##my_prot6 = Seq("CCCCCCCCCCCC", IUPAC.protein)
##my_prot7 = Seq("AAAAAAAAAAAA", IUPAC.protein)

#prots = [my_prot1, my_prot2, my_prot3, my_prot4, my_prot5, my_prot6, my_prot7]
#print("type(prots): ", type(prots));
##type(prots):  <class 'list'>

#print("Our initial motif block")
#for i in prots:
	#print(i);
	#print("type(prots[i]): ", type(i));
	##type(prots[i]):  <class 'Bio.Seq.Seq'>


from Bio import AlignIO

alignment = AlignIO.read(open("TIGR00056_Mot.sto"), "stockholm")
#print("type(alignment): ", type(alignment));
#type(alignment):  <class 'Bio.Align.MultipleSeqAlignment'>
#print("Alignment length %i" % alignment.get_alignment_length())
#for record in alignment :
    #print(record.seq + " " + record.id)
    #print("type(record.seq): ", type(record.seq));
    #type(record.seq):  <class 'Bio.Seq.Seq'>
    #print("type(record.id): ", type(record.id));
    #type(record.id):  <class 'str'>
    


