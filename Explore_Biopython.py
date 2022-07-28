#%%
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from colorama import Back, Style, Fore
from Bio.SeqUtils import GC


# %%

def load_genome(file, n_sequences, display=False) :
    sequences = []
    count = 0
    for seq_record in SeqIO.parse(file,'fasta') :
        if (count<n_sequences) :
            sequences.append(seq_record)
            if display :
                print("Id: " + seq_record.id + " \t" + "Length: " + str("{:,d}".format(len(seq_record))))
                print(repr(seq_record.seq) + "\n")
            count += 1
    return sequences

sequences = load_genome('data/genome.fa',10,True)


# %%

color_switcher = {
    'A': Back.GREEN,
    'a': Back.GREEN,
    'C': Back.YELLOW,
    'c': Back.YELLOW,
    'G': Back.RED,
    'g': Back.RED,
    'T': Back.BLUE,
    't': Back.BLUE,
    ' ': Style.RESET_ALL
}

def nucleotide_seq(genome, len_sequence) :
    genes = []
    for i, nucleotide in enumerate(genome) :
        if i != 0 and i%len_sequence == 0 :
            genes.append(' ')
        genes.append(nucleotide)
    return genes

def display(genome, len_sequence) :
    genome = nucleotide_seq(genome, len_sequence)
    line_break = 0
    for i in range(len(genome)) :
        if genome[i] == ' ':
            line_break += 1
            if line_break != 0 and line_break%6 == 0 :
                text = '\n'
            else :
                text = color_switcher[genome[i]] + genome[i]
        else :
            text = color_switcher[genome[i]] + genome[i]
        print(text, end="")
    Style.RESET_ALL


print('Example of genome sequencing : \n')
display(sequences[0][0:300],10)




# %%
chr2L = load_genome(1)[0].seq

print('Length:\t' + str(len(chr2L)))
print("G Count:\t" + str(chr2L.count('G')))
print("GC%:\t\t" + str(100*(chr2L.count('G') + chr2L.count('C')) / len(chr2L)))




# %%

print("GC% Package:\t" + str(GC(chr2L)))
print("GgCcSs%:\t" + str(100*(chr2L.count("G") + chr2L.count("g") + 
                            chr2L.count("C") + chr2L.count("c") +
                            chr2L.count("S") + chr2L.count("s")) / len(chr2L)))


# %%

chr2Lshort = chr2L[0:21]

print("Orignial :" + chr2Lshort)
print("Complement: " + chr2Lshort.complement())
print("Reverse Complement: " + chr2Lshort.reverse_complement())



# %%

coding_dna = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
print("Coding DNA: \t" + coding_dna)
messenger_rna = coding_dna.transcribe()
print("Messenger RNA: \t" + messenger_rna)
print("Protein Sequence: \t\t" + messenger_rna.translate())
print("Direct Protein Sequence: \t" + coding_dna.translate() + "\n")

print("Vertebrate Mitochondrial Translation: " + coding_dna.translate(table="Vertebrate Mitochondrial"))





# %%

print("Standard translation with Stop: \t" + coding_dna.translate(to_stop=True))
print("Special Translation with Stop: \t\t" + coding_dna.translate(table=2, to_stop=True))




# %%


