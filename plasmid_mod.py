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
re_loc = []


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
			rm_hit = [i[1], i[3], t.group(0), t.span()]
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

re_loc = []
meth_loc = []
nick_loc = []
hom_loc = []

for re in re_sites:
	re_loc.append(re[3])
re_loc1 = list(set(re_loc))
for met in meth_sites:
	meth_loc.append(met[3])
meth_loc1 = list(set(meth_loc))
for n in nick_sites:
	nick_loc.append(n[3])
nick_loc1 = list(set(nick_loc))
for h in homing:
	hom_loc.append(h[3])
hom_loc1 = list(set(hom_loc))

# RM SITE IN ORF
re_orf = []
meth_orf = []
nick_orf = []
hom_orf = []

for orf in orfs:
	for i in range((orf[0]-1), (orf[0]-1)+orf[1]):
		for re in re_loc1:
			if re[0] == i:
				re_o = [re, orf]
				re_orf.append(re_o)
		for meth in meth_loc1:
			if meth[0] == i:
				meth_o = [meth, orf]
				meth_orf.append(meth_o)
		for nick in nick_loc1:
			if nick[0] == i:
				n_o = [nick, orf]
				nick_orf.append(n_o)
		for h in hom_loc1:
			if h[0] == i:
				h_o = [h, orf]
				hom_orf.append(h_o)

# for n in re_orf:
# 	print n
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
