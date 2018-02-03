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


for alignment2 in AlignIO.parse("TIGR00056_Mot(copy).sto", "stockholm"):
	print("Alignment of length %i" % alignment2.get_alignment_length())
	print("type of: ", type(alignment2));
	for whatever in alignment2:
		print("type of whatever's inside alignment2: ", type(whatever));
		print(whatever.seq, "   ", type(whatever.seq));
	#sprint(alignment2.seq)
print(type(AlignIO.parse("TIGR00056_Mot(copy).sto", "stockholm")));


my_prot1 = Seq("GMVLGLQGYVVL", IUPAC.protein)
my_prot2 = Seq("GMVLGLQGYLVL", IUPAC.protein)
my_prot3 = Seq("GMVFTIQVAREF", IUPAC.protein)
my_prot4 = Seq("GMIFTIQVAREF", IUPAC.protein)
my_prot5 = Seq("GGVIALQTYSTF", IUPAC.protein)
my_prot6 = Seq("GVVIAYQSAVQL", IUPAC.protein)
my_prot7 = Seq("GFAVALQGALQL", IUPAC.protein)

#my_prot1 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot2 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot3 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot4 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot5 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot6 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot7 = Seq("AAAAAAAAAAAA", IUPAC.protein)

prots = [my_prot1, my_prot2, my_prot3, my_prot4, my_prot5, my_prot6, my_prot7]
print("type(prots): ", type(prots));
#type(prots):  <class 'list'>

print("Our initial motif block")
for i in prots:
	print(i);
	print("type(prots[i]): ", type(i));
	#type(prots[i]):  <class 'Bio.Seq.Seq'>


from Bio import AlignIO

alignment = AlignIO.read(open("TIGR00056_Mot.sto"), "stockholm")
#print("type(alignment): ", type(alignment));
#type(alignment):  <class 'Bio.Align.MultipleSeqAlignment'>
#print("Alignment length %i" % alignment.get_alignment_length())
for record in alignment :
    print(record.seq + " " + record.id)
    #print("type(record.seq): ", type(record.seq));
    #type(record.seq):  <class 'Bio.Seq.Seq'>
    #print("type(record.id): ", type(record.id));
    #type(record.id):  <class 'str'>
    
input_handle = open("TIGR00056_Mot.sto", "rU")
#print("type(input_handle): ", type(input_handle));
#type(input_handle):  <class '_io.TextIOWrapper'>
output_handle = open("testnew.sto", "w")
#print("type(output_handle): ", type(output_handle));
#type(output_handle):  <class '_io.TextIOWrapper'>

alignments = AlignIO.parse(input_handle, "stockholm")
#print("type(alignments): ", type(alignments));
#type(alignments):  <class 'generator'>
AlignIO.write(alignments, output_handle, "stockholm")

output_handle.close()
input_handle.close()


#hmmbuild testnew.hmm /home/sunlight/Documents/bioinformatics/multiple-alignment/bioscripts/testnew.sto
#holy crap. That command worked! We just made a hmm profile!
#we now want to make a {#type(alignments):  <class 'generator'>} out of {type(prots[i]):  <class 'Bio.Seq.Seq'>}

#hmmsearch globins4.hmm uniprot sprot.fasta > globins4.out
#hmmsearch /home/sunlight/Documents/bioinformatics/multiple-alignment/bioscripts/testnew.hmm /home/sunlight/Documents/bioinformatics/multiple-alignment/bioscripts/TIGR00056_Mot.sto > testnew.out
#wow. That also worked! But we searched this against a multiple sequence alignment, not against the individual sequences themselves. 
