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
target_enz =[]
target_site = []
cut = []
rm_site = []
meth_sites = []
re_sites = []
nick_sites = []
homing = []


for i in enz_list:
	if 'Arthrospira platensis' in i[0]:
		if '?' not in i[2]:
        		init_site.append(i)

for i in init_site:
	i[2] = i[2].replace('R', '[GA]')
	i[2] = i[2].replace('Y', '[CT]')
	i[2] = i[2].replace('M', '[AC]')
	i[2] = i[2].replace('K', '[GT]')
	i[2] = i[2].replace('S', '[GC]')
	i[2] = i[2].replace('W', '[AT]')
	i[2] = i[2].replace('B', '[GCT]')
	i[2] = i[2].replace('D', '[AGT]')
	i[2] = i[2].replace('H', '[ACT]')
	i[2] = i[2].replace('V', '[ACG]')
	i[2] = i[2].replace('N', '[ACGT]')
	i[2] = i[2].replace('^', '')
	i[2] = i[2].replace('\r', '')
	target_enz.append(i)

for i in target_enz:
	cut_site = regex.finditer(i[2], full_seq_str, overlapped=True)
	for c in cut_site:
		target_site.append(c)

for t in target_site:
	for i in target_enz:
		if i[2] in str(t.re):
			rm_hit = [i[1], i[3], t.group(0), t.start()]
			rm_site.append(rm_hit)

for rm in rm_site:
	if 'methyltransferase' in rm[1]:
		meth_sites.append(rm)
	if 'restriction' in rm[1]:
		re_sites.append(rm)
	if 'Nicking' in rm[1]:
		nick_sites.append(rm)
	if 'Homing' in rm[1]:
		homing.append(rm)



# print len(target_site)
# print len(rm_site)



# for t in target_enz:
# 	if t[2] in i in target_site:
# 		print str(i.re)
# # print target_site
#
# # print enz_list
#
# for i in target_site:
# 	print i.re
