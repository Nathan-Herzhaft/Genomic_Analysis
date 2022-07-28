#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import Bio
from Bio import SeqIO




#%%

root = 'data/genes-'
augustus = pd.read_csv(root + 'augustus.csv')
genscan = pd.read_csv(root + 'genscan.csv')
refseq = pd.read_csv(root + 'refseq.csv')
xeno_refseq = pd.read_csv(root + 'xeno-refseq.csv')
ensembl = pd.read_csv(root + 'ensembl.csv')




#%%

def load_genome(n_sequences, display=False) :
    sequences = []
    count = 0
    for seq_record in SeqIO.parse('data/genome.fa','fasta') :
        if (count<n_sequences) :
            sequences.append(seq_record)
            if display :
                print("Id: " + seq_record.id + " \t" + "Length: " + str("{:,d}".format(len(seq_record))))
                print(repr(seq_record.seq) + "\n")
            count += 1
    return sequences


def load_chromosome(name, genome) :
    for sequence in genome :
        if sequence.name == name :
            return sequence
    print('Chromosome not found in genome')
    return None


def protein_sequence(chromosome) :
    coding_dna = chromosome.seq
    protein_sequence = coding_dna.translate(to_stop=True)
    return protein_sequence


genome = load_genome(6)
chromosome = load_chromosome('chr2L',genome)
len(protein_sequence(chromosome))
    





# %%
