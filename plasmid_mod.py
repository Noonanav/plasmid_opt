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

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

with open('enz_list.json') as data_file:
	enz_list = json.load(data_file)

#OPEN SEQUENCE FILE
for record in SeqIO.parse("pCYAko.fasta","fasta"):
    plasmid = record

seq = plasmid.seq
comseq = seq.complement()

full_seq = seq + comseq
full_seq_str = str(full_seq).upper()

# LOCATE ORFs
stop_codon_list = []
stop_codon = regex.finditer(r'(?!TGG)(T[AG][AG])', full_seq_str, overlapped=True)
for match in stop_codon:
    stop_codon_list.append(match)

frame1 = []
frame2 = []
frame3 = []

plasmid_range = range(len(full_seq_str) / 3)
for p in plasmid_range:
    frame1.append(0 + (p * 3))
    frame2.append(1 + (p * 3))
    frame3.append(2 + (p * 3))

f1_codon = []
f2_codon = []
f3_codon = []

for stop in stop_codon_list:
    if stop.start() in frame1:
        f1_codon.append(stop)
    if stop.start() in frame2:
        f2_codon.append(stop)
    if stop.start() in frame3:
        f3_codon.append(stop)

orf_list = []

for i in range(len(f1_codon) - 1):
    orf = (f1_codon[i+1].start() - f1_codon[i].start())
    if orf >= 300:
        o = [f1_codon[i],orf]
        orf_list.append(o)

for i in range(len(f2_codon) - 1):
    orf = (f2_codon[i+1].start() - f2_codon[i].start())
    if orf >= 300:
        o = [f2_codon[i],orf]
        orf_list.append(o)

for i in range(len(f3_codon) - 1):
    orf = (f3_codon[i+1].start() - f3_codon[i].start())
    if orf >= 300:
        o = [f3_codon[i],orf]
        orf_list.append(o)

orfs = []
for orf in orf_list:
    o = full_seq_str[orf[0].end() : (orf[0].end() + orf[1])]
    o_seq = Seq(o)
    orfs1 = [orf[0].end()+1, orf[1], o_seq, o_seq.translate(table=11)]
    orfs.append(orfs1)

# SeqFeature FROM ORF
orf_features = []
for orf in orfs:
	orf_feature = SeqFeature(FeatureLocation(orf[0], orf[0] + orf[1]), type='ORF')
	orf_features.append(orf_feature)

plasmid.features = orf_features

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
rm_loc = []


for i in enz_list:
	if 'Synechocystis species PCC 6803' in i[0]:
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
			rm_hit = [i[1], i[3], t.group(0), t.start(), t.end()]
			rm_site.append(rm_hit)

for rm in rm_site:
	rm_loc.append(rm[3])
rm_loc1 = list(set(rm_loc))

# RM SITE IN ORF
rm_orf = []
rm_nonc = []
orf_ranges = []


for orf in orfs:
	for i in range((orf[0]-1), (orf[0]-1)+orf[1]):
		orf_ranges.append(i)
		for rm in rm_loc1:
			if rm == i:
				rm_o = [rm, orf]
				rm_orf.append(rm_o)

# RM SITE IN NON-CODING PLASMID
orf_bp = list(set(orf_ranges))
non_coding = list(range(len(full_seq_str) - 1))
for x in orf_bp:
	non_coding.remove(x)
for nc in non_coding:
	for rm in rm_loc1:
		if rm == nc:
			rm_nonc.append(rm)
# print rm_nonc
# print rm_orf
# print len(rm_orf)
# print len(rm_nonc)
# print type(non_coding)

# IN-FRAME RM SEQUENCE/LOCATION
site = []
for n in rm_orf:
	x = (n[0])-(n[1][0])
	for i in range(3):
		if (x + i) % 3 == 0:
			seq = n[1][2][(x-(3-i)):((x-(3-i))+15)]
			aa = n[1][3][((x-(3-i))/3):((x-(3-i))/3)+5]
			for rm in rm_site:
				if n[0] == rm[3]:
					site1 = [rm[0], rm[1], rm[2], seq, n[0]-(3-i)]
					site.append(site1)
			# print seq.translate()
			# print aa
# print site
# for s in site:
# 	print [s[0], s[2], s[3].translate(), s[4]]
# 	print s[3].complement()
# 	print (s[4]-len(plasmid))
	# if x % 3 == 0:
	# 	print 'fuck ya'
	# else:
	# 	print 'fuck no'
# for meth in meth_loc1:
# 	print meth

# for orf in orfs:
# 	print orf

# for re in re_loc:
# 	print re
# print len(re_loc)
# print len(re_loc1)
#
# print len(meth_loc)
# print len(meth_loc1)


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
