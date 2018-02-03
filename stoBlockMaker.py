import re
import copy
import math, string, sys
from Bio.Seq import Seq #Seq()
from Bio.Alphabet import IUPAC #IUPAC.protein
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import Bio.Alphabet
from Bio import AlignIO
#rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT", IUPAC.protein), id="1JOY", name="EnvZ", description="Homodimeric domain of EnvZ from E. coli", annotations= );

def pullApart(p_seqName):
	seqNameRegex = re.compile(r'([A-Z0-9|_]+)/(\d+)-(\d+)');
	mo = seqNameRegex.search(p_seqName);
	f_seqName = str(mo.group(1));
	f_seqBeg = int(mo.group(2));
	f_seqEnd = int(mo.group(3));
	return(f_seqName, f_seqBeg, f_seqEnd);

def seqRecordSlice(p_seqRec, p_snip1, p_snip2): 
	"""Takes a SeqRecord and where to slice it, and it returns a new SeqRecord with the start and end values properly adjusted. 
	This assumes we have the start, end, and original items in the annotations.
	p_snip1 and p_snip2 are the ints for a string slice [p_snip1:p_snip2] on the aligned sequence string
	This could be done with one absurdly long statement:
	return SeqRecord(Seq(str(p_seqRec.seq)[p_snip1:p_snip2], SingleLetterAlphabet()), id=(p_seqRec.id), description=(p_seqRec.description), annotations={'start': (p_seqRec.annotations)['start']+len(str(str(p_seqRec.seq)[0:p_snip1]).replace("-","")), 'end': (p_seqRec.annotations)['end']-len(str(str(p_seqRec.seq)[p_snip2:len(str(p_seqRec.seq))]).replace("-","")), 'original': str(p_seqRec.seq)[p_snip1:p_snip2].replace("-","")};"""
	crust1 = str(p_seqRec.seq)[0:p_snip1]; #this is everything to the left of the slice we're taking from the aligned sequence string
	crust2 = str(p_seqRec.seq)[p_snip2:len(str(p_seqRec.seq))]; #this is everything to the right of the slice we're taking from the aligned sequence string
	f_c1 = len(str(crust1).replace("-","")); #this number is how much we need to add to the start# of the sequence
	f_c2 = len(str(crust2).replace("-","")); #this number is how much we need to subtract from the end# of the sequence
	#f_seqSeg = Seq(str(p_seqRec.seq)[p_snip1:p_snip2], SingleLetterAlphabet());
	f_seqSeg = Seq(str(p_seqRec.seq)[p_snip1:p_snip2]);
	f_segRec = SeqRecord(f_seqSeg, id=(p_seqRec.id), description=("a sequence segment in an initial motif"), annotations={'start': (p_seqRec.annotations)['start']+f_c1, 'end': (p_seqRec.annotations)['end']-f_c2, 'original': str(f_seqSeg).replace("-","")});
	return f_segRec;

def avg(p_list):
	return sum(p_list)/len(p_list);

#helper function for findBlocks
def motPull(numlist, t1, t2, p_index, strtBarrier):
	#default strtBarrier is 0 because every list starts at 0 and we can't access elements previous to it (if any somehow existed)
	curAverage = numlist[p_index];
	start = p_index; #includes
	f_end = p_index; #includes
	lAvr = 0; #stores the average score of a left expansion
	rAvr = 0; #stores the average score of a right expansion
	#we stop if we cannot expand selection
	while((start>strtBarrier) or (f_end<(len(numlist)-1))):
		#avg(numlist[start:f_end+1]) #average of current selection
		#avg(numlist[start:f_end+2]) #average of selection expanded+1 right
		#avg(numlist[start-1:f_end+1]) #average of selection expanded+1 left

		#is expansion of motif in both directions possible?
		if((start>strtBarrier) and (f_end<(len(numlist)-1))):
			#is expansion in both directions tolerated?
			if(avg(numlist[start-1:f_end+1])<t2 and avg(numlist[start:f_end+2])<t2):
				if(numlist[start-1]>=numlist[f_end+1]):
					f_end = f_end + 1;
				else:
					start = start - 1;
			#is expanding left tolerated?
			elif(avg(numlist[start-1:f_end+1])<t2):
				start = start - 1;
			#is expanding right tolerated?
			elif(avg(numlist[start:f_end+2])<t2):
				f_end = f_end + 1;
			else:
				break;
		#is expanding left possible?
		elif(start>strtBarrier):
			#is expanding left tolerated?
			if(avg(numlist[start-1:f_end+1])<t2):
				start = start - 1;
			else:
				break;
		#is expanding right possible?
		elif(f_end<(len(numlist)-1)):
			#is expanding right tolerated?
			if(avg(numlist[start:f_end+2])<t2):
				f_end = f_end + 1;
			else:
				break;
		#is expansion impossible?
		elif(start==strtBarrier and f_end==(len(numlist)-1)):
			break;
		else:
			print("unknown logic error");
		
	#now we need to shave the high entropy ends off the the motif block
	#shave left side
	while(numlist[start]>t1):
		start= start+1;
	#shave right side
	while(numlist[f_end]>t1):
		f_end= f_end-1;
	
		#return the valuepair of start value and length of the motif
	return[start,(f_end-start+1)]

def findBlocks(scores, t1, t2, minsize): #t1 is lower best results threshold, t2 is higher and more lenient
	motList = [];
	ignore = 0
	firstAvail = 0 #after each successful motif result, we build a barrier to prevent new motifs overlapping with previous ones
	tempmot = []; #used to temporary store valuepairs to put into motList
	for i in range(len(scores)):
		if(ignore<=0):
			if(scores[i] < t1):
				tempmot = motPull(scores, t1, t2, i, firstAvail);
				if(tempmot[1]>=minsize):
					motList.append(tempmot); #if the motif is long enough to keep, it's added to our list
					firstAvail = tempmot[1]+tempmot[0] #should be the last index of found motif+1, start+length = last index + 1
					ignore = firstAvail-i; #since the motif expanded ahead, we want to skip to where it ended
				tempmot = [];
		else:
			ignore = ignore-1; #ticks down until we reach a spot we haven't checked yet
	return motList;

#possible help: http://blog.dkbza.org/2007/05/scanning-data-for-entropy-anomalies.html
#http://pythonfiddle.com/shannon-entropy-calculation/
# Stolen from Ero Carrera
# http://blog.dkbza.org/2007/05/scanning-data-for-entropy-anomalies.html
def range_bytes (): return range(256)
def range_printable(): return (ord(c) for c in string.printable)
def H(data, iterator=range_bytes):
	if not data:
		return 0;
	entropy = 0;
	#print(data,end='');
	for x in iterator():
		p_x = float(data.count(chr(x)))/len(data)
		if p_x > 0:
			if(chr(x)!='-'): #I made this modification
				entropy += - p_x*math.log(p_x, 2)
			else:
				entropy += - p_x*math.log((float(1)/len(data)), 2)
				#print(chr(x),end=" ");
	return entropy

def getEntropy(xList):
	"""this calculates the columnwise entropy of the multiple alignment
	it returns a list of the entropy values of each column
	"""
	#for item in [(''.join([xList[j][i] for j in range(len(xList))])) for i in range(len(xList[0]))]:
		#print(item, ": ",H(item, range_printable), sep='');
		#print([H(item, range_printable) for item in [(''.join([xList[j][i] for j in range(len(xList))])) for i in range(len(xList[0]))]]);
	return([H(item, range_printable) for item in [(''.join([xList[j][i] for j in range(len(xList))])) for i in range(len(xList[0]))]]);

def getBlocks(source, p_e1=1.2, p_e2=1.5, p_e3=10):
	sequenceList = [str(i_rec.seq) for i_rec in source];
	places = findBlocks(getEntropy(sequenceList), p_e1, p_e2, p_e3);
	motList = [];
	currentMot = [];
	currentSeq = "";
	for mot in places: #this and the stuff below it needs heavy modification
		for sequence in source:
			currentSeq = seqRecordSlice(sequence, mot[0], (mot[0]+mot[1]));
			currentMot.append(currentSeq);
		motList.append(currentMot);
		currentMot = [];
	return copy.deepcopy(motList);

from Bio.SubsMat import MatrixInfo as matlist
import numpy as np
from scipy.cluster import hierarchy
from random import randint
def clusteringStep(p_prots):
	matrix = matlist.blosum62 
	gap_open = -11
	gap_extend = -1
	scoreLists = [];
	
	for i in range(len(p_prots)): #I think this creates our score matrix.
		sl = [];
		for j in range(len(p_prots)):
			for a in pairwise2.align.globaldx(p_prots[j], p_prots[i], matrix):
				al1,al2, score, begin, end = a
			sl.append(score);
		scoreLists.append(sl);
	
	mymat = np.matrix([[int(x) for x in item] for item in scoreLists])
	#print("ourScoreMatrix:");
	#print(mymat);
	
	f_clusteringResults=hierarchy.fclusterdata(mymat, 1.0,depth=4, method='single'); #single,complete,average,weighted,centroid,median,ward
	print("fclusterdata(ourScoreMatrix,10.0): ", f_clusteringResults);
	#print("total clusters entries: ",len(f_clusteringResults));
	return f_clusteringResults;
	
	
def choosenCLists(p_clst):
	"""return a list of indexes that match our chosen result, and a list of indexes that don't.
	p_clst: parameter cluster list. will be something like [2 2 1 1 2 2 2] this is a list of entries and the numbers are which cluster they are assigned to"""
	
	f_ind=p_clst[randint(0, len(p_clst)-1)]; #this is the weighted randomly chosen cluster.
	f_included = [];
	f_excluded = [];
	for i in range(len(p_clst)):
		if p_clst[i]==f_ind:
			f_included.append(i);
		else:
			f_excluded.append(i);
	return (f_included, f_excluded);

def ourMain(originFile="TIGR00056.sto"):
	recList = [];
	for thing in AlignIO.parse(originFile, "stockholm"):
		for record in thing: #this extracts records from multiple sequence alignments from a stockholm file.
			recList.append(SeqRecord(record.seq, id=pullApart(record.id)[0], description="An aligned sequence", annotations={'start': pullApart(record.id)[1], 'end': pullApart(record.id)[2], 'original': str(record.seq).replace("-","")}));

	#entScores = getEntropy(sequenceList);

	#slsize = extractAligned(originFile)[2]; #this accounts for name lengths for sequence lines
	cores = getBlocks(recList); #this is a list of motifs. each motif has a list of sequence segments
	"""for biip in cores: #show off the adjust sequence segmment records
		for item in biip:
			print(len(item));
			print(item.id);
			print("end, start:" ,(item.annotations)['end'],(item.annotations)['start']);
			print(item.seq)"""
	prots = [i_seq.seq for i_seq in cores[0]];		
	#print("Our initial motif block")
	#for i in prots:
	#	print(i);

	clusteringResults = clusteringStep(prots);

	#multiple assignment trick
	included_entries, excluded_entries = choosenCLists(clusteringResults);
	print("included_entries", included_entries); 
	for ind in included_entries:
		print(prots[ind]);
	#these will be used to create an hmm profile
	from Bio.Align import MultipleSeqAlignment
	littleChoosen = MultipleSeqAlignment([]);
	protRecs = [i_seq for i_seq in cores[0]];
	for i_ind in included_entries:
		littleChoosen.append(protRecs[i_ind]); #filling MultipleSeqAlignment with included sequence segment records
	#for i_rec in littleChoosen:
		#print(i_rec.seq);
	print(littleChoosen.format("stockholm"));
	fo = open("chosenMotif.sto", "w");
	fo.write(littleChoosen.format("stockholm"));
	fo.close();
	
	print("excluded_entries", excluded_entries);
	for ind in excluded_entries:
		print(prots[ind]);
	#these will be the full sequences the profile searches against
	bigSoap = MultipleSeqAlignment([]);
	for i_ind in excluded_entries:
		bigSoap.append(recList[i_ind]); #filling MultipleSeqAlignment with excluded full sequence records
	fo = open("searchUs.sto", "w");
	fo.write(bigSoap.format("stockholm"));
	fo.close();

	#hmmbuild chosenProfile.hmm /home/sunlight/Documents/bioinformatics/multiple-alignment/bioscripts/chosenMotif.sto
	#hmmsearch /home/sunlight/Documents/bioinformatics/multiple-alignment/bioscripts/chosenProfile.hmm /home/sunlight/Documents/bioinformatics/multiple-alignment/bioscripts/searchUs.sto > testResults.out

ourMain("TIGR00056.sto");