"""
Simple method for classifying eukaryotic protein sequences into 1 of 5 categories.
The number of protein sequences for each of the 5 classes are shown in brackets.
1. Cytosolic (2463)
2. Extracellular/Secreted (1236)
3. Nuclear (2736)
4. Mitochondrial (1023)
5. None of the above (2002)
"""
from Bio import SeqIO
import numpy as np
from Bio.Data import IUPACData
from src import features
import matplotlib.pyplot as plt
from src.features import Features


raw_seq_count = {'cyto': 2463, 'secreted': 1236, 'mito': 1023, 'nuclear': 2736, 'other': 2002}


def read_and_validate_all_fasta_and_label(subset: str):
    """
    (For hard-coded relative path, assumes current working directory is `src`).
    :return:
    """
    valid_seqs, seq_lengths = list(), list()
    sequence_count = len(list(SeqIO.parse(f'../datasets/{subset}.fasta', 'fasta')))
    print(f'Expecting {raw_seq_count[subset]} cytosolic: {sequence_count}')

    for record in SeqIO.parse(f'../datasets/{subset}.fasta', 'fasta'):

        if set(record.seq) <= set(IUPACData.protein_letters):
            sequence = record.seq
            valid_seqs.append(sequence)
            seq_len = len(sequence)
            print(f'seq_len {seq_len}')
            seq_lengths.append(seq_len)
            # global_aa_comp = features.global_aa_comp(sequence)
            # print(global_aa_comp)
            local_aa_comp = features.local_aa_comp(sequence)

            # Global amino acid composition (i.e. percentages of all 20 amino acids present in whole sequence)
            # Local amino acid composition (i.e. over first 50 amino acids or last 50 amino acids)
            # Isoelectric point & molecular weight (e.g. see HERE or HERE)
            # Specific sequence patterns near the start or near the end of the sequence

        else:
            print(f'ID {record.id}, this {subset} sequence {record.seq} has invalid amino acid character. '
                  f'Therefore it is excluded from this dataset)')

    return valid_seqs, np.array(seq_lengths)





if __name__ == '__main__':

    import os
    print(os.getcwd())
    # read_all_fasta_and_label()
    # standard_aa = set(IUPACData.protein_letters)
    # print(len(standard_aa))

    cyto, cyto_lengths = read_and_validate_all_fasta_and_label('cyto')

    x = list(range(len(cyto)))
    y = sorted(cyto_lengths)
    plt.figure(figsize=(10, 6))  # Set the figure size
    plt.scatter(x, y, alpha=0.6)  # Plot the data
    plt.title('Distribution of Integer Values')  # Set the title of the plot
    plt.xlabel('Index')  # Set the x-axis label
    plt.ylabel('Integer Values')  # Set the y-axis label
    plt.grid(True)  # Show a grid for easier visualization
    plt.show()

    # mito = read_and_validate_all_fasta_and_label('mito')
    # nuclear = read_and_validate_all_fasta_and_label('nuclear')
    # other = read_and_validate_all_fasta_and_label('other')
    # secreted = read_and_validate_all_fasta_and_label('secreted')
