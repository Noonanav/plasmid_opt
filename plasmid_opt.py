# REQUIRED INPUTS
# Plasmid sequence file (line 15)
# Entrez.email (line 35)
# Genome ID (line 38)
# Assembly ID

import re
import regex
import json
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC

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


# ESEARCH HOST GENOME IDLIST

# Entrez.email = "avery.noonan@mail.utoronto.ca"

# gen_search = Entrez.esearch(db="nucleotide", term="GCF_000009725.1", report="full")
# genome_search = Entrez.read(gen_search)
#
# open('./genome_IdList.json', 'w').write(json.dumps(genome_search['IdList']))

# EFETCH HOST GENOME SEQUENCES
#
sequence_list = []
# gene_list = []
#
with open('genome_IdList.json') as data_file:
	genome_IdList = json.load(data_file)
#
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

# print(gene_list[0])

# DICTIONARY - GENOME_IDLIST : [GENE, NUCLEOTIDE]

with open('sequence_list.json') as data_file:
	sequence_list = json.load(data_file)
with open('gene_list.json') as data_file:
	gene_list = json.load(data_file)
#
# pair_list = []
#
# for num in range(5):
#     seq_pair = [gene_list[num], sequence_list[num]]
#     pair_list.append(seq_pair)
#
# genome_dict = dict(zip(genome_IdList, pair_list))
# # print(genome_dict[genome_IdList[0]][1])
# # print(len(genome_dict))

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
