# REQUIRED INPUTS
# Plasmid sequence file (line 15)
# Entrez.email (line 35)
# Genome ID (line 38)
# Assembly ID

import re
import regex
import json
import os
from Bio import SeqIO
from Bio import SearchIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC
from Bio.Blast import NCBIWWW

# PLASMID INPUT
for record in SeqIO.parse("pCYAko.fasta","fasta"):
    plasmid = record

seq = plasmid.seq
comseq = seq.complement()

full_seq = seq + comseq
full_seq_str = str(full_seq).upper()

# 6 BP PLASMID SEGMENTS

substrings = {}
seq_len = len(full_seq_str)

for i in range(3,11):
    if i not in substrings:
        substrings[i] = []
    for j in range(0,seq_len-i-1):
        substrings[i].append(full_seq_str[j:j+i])

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

# ESEARCH HOST GENOME IDLIST

# Entrez.email = "avery.noonan@mail.utoronto.ca"

# gen_search = Entrez.esearch(db="nucleotide", term="GCF_000009725.1", report="full")
# genome_search = Entrez.read(gen_search)
#
# open('./genome_IdList.json', 'w').write(json.dumps(genome_search['IdList']))

# EFETCH HOST GENOME SEQUENCES
#
# sequence_list = []
# gene_list = []
#
with open('genome_IdList.json') as data_file:
	genome_IdList = json.load(data_file)



# GENEBANK FILES
subdir = "genbank_dir"
# try:
#     os.mkdir(subdir)
# except Exception:
#     pass
# #
# for id in genome_IdList:
#     filename = id
#     if not os.path.isfile(filename):
#         net_handle = Entrez.efetch(db="nucleotide", id=id, rettype="gbwithparts", retmode="text")
#         out_handle = open(os.path.join(subdir, filename), "w")
#         out_handle.write(net_handle.read())
#         out_handle.close()
#         net_handle.close()

# for id in genome_IdList:
#     nuc_fetch = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
#     nucleotide_fetch = nuc_fetch.read()
#     sequence_list.append(nucleotide_fetch)
#
# for id in genome_IdList:
#     gen_fetch = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
#     genome_fetch = gen_fetch.read()
#     gene_list.append(genome_fetch)
#
# open('./sequence_list.json', 'w').write(json.dumps(sequence_list))
# open('./gene_list.json', 'w').write(json.dumps(gene_list))

# DICTIONARY - GENOME_IDLIST : [GENE, NUCLEOTIDE]

with open('sequence_list.json') as data_file:
	sequence_list = json.load(data_file)
with open('gene_list.json') as data_file:
	gene_list = json.load(data_file)

pair_list = []

for num in range(5):
    seq_pair = [gene_list[num], sequence_list[num]]
    pair_list.append(seq_pair)

genome_dict = dict(zip(genome_IdList, pair_list))

# GENOME PARSING

gene_dicts = []
cds_list = []
AA_seqs = []

for id in genome_IdList:
    for gene in SeqIO.parse(os.path.join(subdir, id), "genbank"):
        genes = gene.features
        for g in genes:
            t = g.type
            if 'CDS' in t:
                cds_list.append(g)
            else:
                pass

for cds in cds_list:
    gene_dicts.append(cds.qualifiers)

# print gene_dicts[50]
# for trans in gene_dicts:
#     if 'product' in trans:
#         # print(trans['product'])
#     else:
#         pass

for trans in gene_dicts:
    if 'translation' in trans:
        AA_seqs.append(trans['translation'])
    else:
        pass


# BLAST GENOME

# for aa in AA_seqs:
#     blast_result = NCBIWWW.qblast("blastp", "nr", aa)

# X BP GENOME SEGMENTS
# gen_bp_ol = []
#
# full_gen = Seq("", generic_dna)
# for sequence in sequence_list:
#     full_gen += sequence
#     full_gen_str = str(full_gen)
#
# stop = False
# for i in range(3,11):
#     # print "substrings of length " + str(i)
#     for ol in substrings[i]:
#         if ol not in full_gen_str:
#             stop = True
#             gen_bp_ol.append(ol)
#             # print ol
#     if stop:
#         break
#
# open('./gen_bp_overlap.json', 'w').write(json.dumps(gen_bp_ol))

with open('gen_bp_overlap.json') as data_file:
	gen_bp_ol = json.load(data_file)


#
# rseq = seq.reverse_complement()
# # print seq.count("G")
# print type(full_seq_str)
# print rseq
