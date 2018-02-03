import re
import copy
import math, string, sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
#rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT", IUPAC.protein), id="1JOY", name="EnvZ", description="Homodimeric domain of EnvZ from E. coli", annotations= );


def pullApart(p_seqName):
	seqNameRegex = re.compile(r'([A-Z0-9|_]+)/(\d+)-(\d+)');
	mo = seqNameRegex.search(p_seqName);
	f_seqName = str(mo.group(1));
	f_seqBeg = int(mo.group(2));
	f_seqEnd = int(mo.group(3));
	return(f_seqName, f_seqBeg, f_seqEnd);

from Bio import SeqIO

#for record in SeqIO.parse("TIGR00056.sto", "stockholm"):
    #print(type(record));
    #print(record.id);
    #print(record.features);
    #print(record.seq);
##for both SeqIO and AlignIO. features is empty

from Bio import AlignIO
#this extracts a multiple sequence alignment from a stockholm file.
"""recList = [];
for thing in AlignIO.parse("TIGR00056.sto", "stockholm"):
	for record in thing:
		recList.append(SeqRecord(record.seq, id=pullApart(record.id)[0], description="An aligned sequence", annotations={'start': pullApart(record.id)[1], 'end': pullApart(record.id)[2], 'original': str(record.seq).replace("-","")}));"""

#we can make this with a list comprehension to look super cool. 
recList = [[SeqRecord(record.seq, id=pullApart(record.id)[0], description="An aligned sequence", annotations={'start': pullApart(record.id)[1], 'end': pullApart(record.id)[2], 'original': str(record.seq).replace("-","")}) for record in thing] for thing in AlignIO.parse("TIGR00056.sto", "stockholm")][0];
"""for item in recList:
	print(item.id);
	print("end-start:" ,(item.annotations)['end']-(item.annotations)['start']);
	print("sequence length:",len((item.annotations)['original']))
	#this shows us that seq.length = (seq.end - seq.start)+1"""

for item in recList:
	print(item.id);
	print("end-start:" ,(item.annotations)['end']-(item.annotations)['start']);
	print("sequence length:",len((item.annotations)['original']))
	print("aligned length:",len(item.seq))



def seqRecordSlice(p_seqRec, p_snip1, snip2): 
	"""Takes a SeqRecord and where to slice it, and it returns a new SeqRecord with the start and end values properly adjusted. 
	This assumes we have the start, end, and original items in the annotations.
	snip1 and snip2 are the ints for a string slice [snip1:snip2] on the aligned sequence string
	This could be done with one absurdly long statement:
	return SeqRecord(Seq(str(p_seqRec.seq)[snip1:snip2], SingleLetterAlphabet()), id=(p_seqRec.id), description=(p_seqRec.description), annotations={'start': (p_seqRec.annotations)['start']+len(str(str(p_seqRec.seq)[0:p_snip1]).replace("-","")), 'end': (p_seqRec.annotations)['end']-len(str(str(p_seqRec.seq)[snip2:len(str(p_seqRec.seq))]).replace("-","")), 'original': str(p_seqRec.seq)[snip1:snip2].replace("-","")};"""
	crust1 = str(p_seqRec.seq)[0:p_snip1]; #this is everything to the left of the slice we're taking from the aligned sequence string
	crust2 = str(p_seqRec.seq)[snip2:len(str(p_seqRec.seq))]; #this is everything to the right of the slice we're taking from the aligned sequence string
	f_c1 = len(str(crust1).replace("-","")); #this number is how much we need to add to the start# of the sequence
	f_c2 = len(str(crust2).replace("-","")); #this number is how much we need to subtract from the end# of the sequence
	f_seqSeg = Seq(str(p_seqRec.seq)[snip1:snip2], SingleLetterAlphabet());
	f_segRec = SeqRecord(f_seqSeg, id=(p_seqRec.id), description=(p_seqRec.description), annotations={'start': (p_seqRec.annotations)['start']+f_c1, 'end': (p_seqRec.annotations)['end']-f_c2, 'original': str(f_seqSeg).replace("-","")});
	return f_segRec;
