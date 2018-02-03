#!/usr/bin/python3

#this file shows how we can create a matrix of pairwise scores and turn the matrix into a cluster

import os
import math, string, sys
from Bio import AlignIO
#alignment = AlignIO.read("TIGR00056.sto", "stockholm")
#print(alignment)

from Bio.Seq import Seq #Seq()
from Bio.Alphabet import IUPAC #IUPAC.protein
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


my_prot1 = Seq("AAAAAAAACCAA", IUPAC.protein)
my_prot2 = Seq("AAAAAAAAAAAA", IUPAC.protein)
my_prot3 = Seq("AAAAAAAAACCA", IUPAC.protein)
my_prot4 = Seq("AAAAAAAAAAAA", IUPAC.protein)
my_prot5 = Seq("CCACCCCCCCCC", IUPAC.protein)
my_prot6 = Seq("CCCCCAAACCCC", IUPAC.protein)
my_prot7 = Seq("CCCCCCACCCCC", IUPAC.protein)
my_prot8 = Seq("CCCCCCACCCCC", IUPAC.protein)

#my_prot1 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot2 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot3 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot4 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot5 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot6 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot7 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot8 = Seq("CCCCCCCCCCCC", IUPAC.protein)

prots = [my_prot1, my_prot2, my_prot3, my_prot4, my_prot5, my_prot6, my_prot7, my_prot8]

print("Our initial motif block")
for i in prots:
	print(i);

from Bio.SubsMat import MatrixInfo as matlist 

matrix = matlist.blosum62 
gap_open = -11
gap_extend = -1

scoreLists = [];
#for i in range(1,len(prots)):
	#sl = [];
	#for j in range(0,i):
		#for a in pairwise2.align.globaldx(prots[j], prots[i], matrix):
			#al1,al2, score, begin, end = a
		#sl.append(score);
	#scoreLists.append(sl);
for i in range(len(prots)):
	sl = [];
	for j in range(len(prots)):
		for a in pairwise2.align.globaldx(prots[j], prots[i], matrix):
			al1,al2, score, begin, end = a
		sl.append(score);
	scoreLists.append(sl);

#print([[int(x) for x in item] for item in scoreLists]);
	
import numpy as np
#mymat = np.matrix([[1, 2], [3, 4]])
mymat = np.matrix([[int(x) for x in item] for item in scoreLists])
#mymat = np.matrix([[56.0], [26.0, 26.0], [25.0, 25.0, 58.0], [28.0, 28.0, 26.0, 25.0], [28.0, 28.0, 26.0, 25.0, 29.0], [33.0, 33.0, 26.0, 26.0, 25.0, 33.0]]);
print("ourScoreMatrix:");
print(mymat);

from scipy.cluster import hierarchy

#from scipy.cluster import hierarchy
#from scipy.cluster.hierarchy import fclusterdata
#chang depth from 2 to 1 makes us only calculate 1 cluster
myclusy=hierarchy.fclusterdata(mymat, 1.0,depth=4, method='single'); #single,complete,average,weighted,centroid,median,ward
print("fclusterdata(ourScoreMatrix,10.0): ", myclusy);
print("total clusters entries: ",len(myclusy));
print("largest result: ",max(myclusy));

#this is a really bad function to find the cluster with the most entries (cluster with largest entries in case of tie)
def mostCommon(lst):
	counts = [];
	items = [];
	for i in range(len(lst)):
		if lst[i] not in items:
			items.append(lst[i]);
			counts.append(1);
		else:
			counts[items.index(lst[i])] = counts[items.index(lst[i])]+1;
	#if we have a tie for most most_common
	mostcount = max(counts);
	if(counts.count(mostcount)>1):
		maxes = [];
		for j in range(len(items)):
			if counts[j]==mostcount:
				maxes.append(items[j]);
		return max(maxes);
	else: 
		return items[counts.index(mostcount)];
		
print("highest of most common: ",mostCommon(myclusy));
#return a list of indexes that match our most common result
def getCommon(lst, most):
	items = [];
	for i in range(len(lst)):
		if lst[i]==most:
			items.append(i);
	return items;
print(getCommon(myclusy, mostCommon(myclusy)));
clusi = getCommon(myclusy, mostCommon(myclusy));
for ind in clusi:
	print(prots[ind]);

#https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
#from scipy.cluster.hierarchy import linkage
#myclus1=hierarchy.linkage(mymat);
#print("linkage(ourScoreMatrix):");
#print(myclus1);

##from scipy.cluster.hierarchy import dendrogram
#plt.figure(figsize=(25, 10))
#plt.title('Hierarchical Clustering Dendrogram')
#plt.xlabel('sample index')
#plt.ylabel('distance')

#dendrogram(
    #myclus1,
    #leaf_rotation=90.,  # rotates the x axis labels
    #leaf_font_size=8.,  # font size for the x axis labels
#)
#plt.show()