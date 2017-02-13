import json
import regex
import string
import random
from Bio import SeqIO
from Bio import SearchIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment


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
target_list = []


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
	i[2] = i[2].replace('(-11/-8)', '')
	i[2] = i[2].replace('^', '')
	i[2] = i[2].replace('\r', '')
	target_enz.append(i)
	target_list.append(i[2])

for i in target_enz:
	cut_site = regex.finditer(i[2], full_seq_str, overlapped=True)
	for c in cut_site:
		target_site.append(c)

for t in target_site:
	for i in target_enz:
		if i[2] in str(t.re):
			rm_hit = [i[1], i[3], t.group(0), t.start(), t.end()]
			# rm_hit = [Enzyme name, Enzyme type, target sequence, target start, target end]
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
				# rm_orf = [rm target site, ORF]

# RM SITE IN NON-CODING PLASMID
orf_bp = list(set(orf_ranges))
non_coding = list(range(len(full_seq_str) - 1))
for x in orf_bp:
	non_coding.remove(x)
for nc in non_coding:
	for rm in rm_loc1:
		if rm == nc:
			rm_nonc.append(rm)

# IN-FRAME RM SEQUENCE/LOCATION
site = []
for n in rm_orf:
	x = (n[0])-(n[1][0])
	for i in range(1,4):
		if (x + i) % 3 == 0:
			seq = n[1][2][(x-(3-i)):((x-(3-i))+15)]
			aa = n[1][3][((x-(3-i))/3):((x-(3-i))/3)+5]
			for rm in rm_site:
				if n[0] == rm[3]:
					site1 = [rm[0], rm[1], rm[2], seq, n[0]-(3-i), 3-i]
					# site1 = [Enzyme name, Enzyme type, Recognition bp, In-frame bp, In-frame target loc]
					site.append(site1)

# SEQUENCE MODIFICATION
loc = []
ran_ins = []
mod_inserts = []


for s in site:
	if s[4] not in loc:
		random_seq = ''.join(random.SystemRandom().choice('ACTG') for _ in range(len(s[2])))
		random_ins = str(s[3]).replace(s[2], random_seq)
		x = regex.search("(" + ")|(".join(target_list) + ")", random_ins, overlapped=True)
		ins_perm = []
		while Seq(random_ins).translate() != s[3].translate() or x != None:
			random_seq = ''.join(random.SystemRandom().choice('ACTG') for _ in range(len(s[2])))
			ins_perm.append(random_seq)
			random_ins = str(s[3]).replace(s[2], random_seq)
			x = regex.search("(" + ")|(".join(target_list) + ")", random_ins, overlapped=True)
			if len(list(set(ins_perm))) == 4**len(s[2]):
				mod_insert = [str(s[3]), s[3], s[2], s[4]]
				mod_inserts.append(mod_insert)
				loc.append(s[4])
				break
		else:
			mod_insert = [random_ins, s[3], s[2], s[4]]
			mod_inserts.append(mod_insert)
			loc.append(s[4])

unmod_seq = full_seq_str

for mod in mod_inserts:
	full_seq_str = full_seq_str.replace(full_seq_str[mod[3]-1:mod[3]+len(mod[0])-1], mod[0])

# alignment = pairwise2.align.globalxx(Seq(unmod_seq), Seq(full_seq_str))
# print(pairwise2.format_alignment(*alignment[0]))


# PLASMID DIAGRAM
# gd_diagram = GenomeDiagram.Diagram("pCYAko")
# gd_track_for_features = gd_diagram.new_track(1, name="ORFs")
# gd_feature_set = gd_track_for_features.new_set()
#
# for orfs in plasmid.features:
#     if orfs.type == "ORF":
#         color = colors.blue
#     else:
#         color = colors.lightblue
#     gd_feature_set.add_feature(orfs, color=color, label=True)
#
# gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4',
#                 fragments=4, start=0, end=len(full_seq_str))
# gd_diagram.write("plasmid_linear.pdf", "PDF")

# gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
#                 start=0, end=len(plasmid), circle_core=0.7)
# gd_diagram.write("plasmid_circular.pdf", "PDF")
