#!/usr/bin/python3

#this file also makes clusters from a matrix of pairwise scores. But it also extracts the entries

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

#my_prot1 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot2 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot3 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot4 = Seq("AAAAAAAAAAAA", IUPAC.protein)
#my_prot5 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot6 = Seq("CCCCCCCCCCCC", IUPAC.protein)
#my_prot7 = Seq("AAAAAAAAAAAA", IUPAC.protein)

prots = [my_prot1, my_prot2, my_prot3, my_prot4, my_prot5, my_prot6, my_prot7]

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
		sl.append((score/len(prots[0])*20)); #score normalizing. Went from (score) to (score/len[0]) to normalize. But that made everything in one cluster
		#so multiplying by 20 readusted it so the clustering function would seperate them again. 
		#actually length is the same for all of these since we used a global align
		#maybe this would only work for a different pairing system that aligns only similar segments. find way to do this and get segment length
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

from random import randint

def getWeightedCluster(p_seg, p_clst): 
	"""This function takes a cluster result from hierarchy.fclusterdata, calculates a weighted random cluster to choose, then returns a list of all items in that cluster, and another list of all the items that aren't.
	p_seg: parameter list of sequence segments in the initial motif block. Each item in this list will be put into one of the two lists we'll be returning. 
	p_clst: parameter cluster list. will be something like [2 2 1 1 2 2 2] this is a list of entries and the numbers are which cluster they are assigned to"""
	f_included=[];
	f_excluded=[];
	print("random number: ",randint(0, len(p_clst)-1), " cluster result: ", p_clst[randint(0, len(p_clst)-1)]);
	f_ind=p_clst[randint(0, len(p_clst)-1)]; #this is the weighted randomly chosen cluster. Now we we'll create lists of entries included and entries excluded in it
	for i in range(len(p_clst)):
		if(p_clst[i]==f_ind): #is this entry in our chosen cluster?
			f_included.append(p_seg[i]); #if so, include that entry's sequence segment in our included list
		else:
			f_excluded.append(p_seg[i]); #otherwise put the sequence segment in the exluded list
	return (f_included, f_excluded);

#return a list of indexes that match our most common result
def getCommon(lst, most):
	items = [];
	for i in range(len(lst)):
		if lst[i]==most:
			items.append(i);
	return items;

print(getCommon(myclusy, mostCommon(myclusy)));
print(getWeightedCluster(prots, myclusy));
clusi = getCommon(myclusy, mostCommon(myclusy));
for ind in clusi:
	print(prots[ind]);




#don't pick, the largest one,
#but choose a random one with the results weighted on probabilities.
#so don't choose the clusterA with 8 entries instead of the clusterB with 2 entries
#choose the clusterA 80% of the time and clusterB 20% of the time.
#so score doesn't in our choice. Just random choice weighted by size.


#build profile and use it to search the other entries then use that to shift the alignment

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