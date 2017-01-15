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
full_seq_str = str(full_seq)

# 6 BP PLASMID SEGMENTS
psixbp = []

sixbpp = regex.findall(r'([a-zA-Z]{6})', full_seq_str, overlapped=True)
psixbp.extend(sixbpp)

# print psixbp

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
# id_num = [0, 1, 2, 3, 4]

# for num in id_num:
#     seq_pair = [gene_list[num], sequence_list[num]]
#     pair_list.append(seq_pair)

# genome_dict = dict(zip(genome_IdList, pair_list))
# print(genome_dict[genome_IdList[0]][0])
# print(len(genome_dict[genome_IdList[2]]))

# 6 BP GENOME SEGMENTS
gensixbp = []

full_gen = Seq("", generic_dna)
for sequence in sequence_list:
    full_gen += sequence
    full_gen_str = str(full_gen)

sixbpgen = regex.findall(r'([a-zA-Z]{6})', full_gen_str, overlapped=True)
gensixbp.extend(sixbpgen)


# print(genome_IdList)
#
# print(pair_list[0])
# genome = Entrez.efetch(db="biosample", id="5831916", rettype="gb", retmode="text")

# print(genome.read())
#
# rseq = seq.reverse_complement()
# # print seq.count("G")
# print type(full_seq_str)
# print rseq
