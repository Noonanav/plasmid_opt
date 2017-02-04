import json
import regex
from Bio import SeqIO
from Bio import SearchIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

with open('enz_list.json') as data_file:
	enz_list = json.load(data_file)

#OPEN SEQUENCE FILE
for record in SeqIO.parse("pCYAko.fasta","fasta"):
    plasmid = record

seq = plasmid.seq
comseq = seq.complement()

full_seq = seq + comseq
full_seq_str = str(full_seq).upper()

# TARGET CODE CHANGE

init_site = []
target_site = []

for i in enz_list:
	if 'Arthrospira platensis' in i[0]:
		if '?' not in i[2]:
        		init_site.append(i[2])

for i in init_site:
	i = i.replace('R', '[GA]')
	i = i.replace('Y', '[CT]')
	i = i.replace('M', '[AC]')
	i = i.replace('K', '[GT]')
	i = i.replace('S', '[GC]')
	i = i.replace('W', '[AT]')
	i = i.replace('B', '[GCT]')
	i = i.replace('D', '[AGT]')
	i = i.replace('H', '[ACT]')
	i = i.replace('V', '[ACG]')
	i = i.replace('N', '[ACGT]')
	i = i.replace('^', '')
	target_site.append(i.replace("\r", ""))

match_str = "(" + ")|(".join(target_site) + ")"

cut_site = regex.finditer(match_str, full_seq_str)

for i in cut_site:
	print i.group(0)
