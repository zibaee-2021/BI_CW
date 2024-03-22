from src import main
import numpy as np
from Bio import SeqIO
from Bio.Data import IUPACData
import matplotlib.pyplot as plt
from src.features import Features

raw_seq_count = {'cyto': 2463, 'secreted': 1236, 'mito': 1023, 'nuclear': 2736, 'other': 2002}


def read_and_validate_all_fasta_and_label(subset_name: str):
    """
    (For hard-coded relative path, assumes current working directory is `src`).
    :return:
    """
    print(f'\n#################{subset_name}#################')
    valid_seqs, seq_lengths = list(), list()
    sequence_count = len(list(SeqIO.parse(f'../datasets/{subset_name}.fasta', 'fasta')))
    print(f'Expecting {raw_seq_count[subset_name]} {subset_name}: {sequence_count}')

    for record in SeqIO.parse(f'../datasets/{subset_name}.fasta', 'fasta'):

        if set(record.seq) <= set(IUPACData.protein_letters):
            sequence = record.seq
            valid_seqs.append(sequence)
            seq_len = len(sequence)
            # print(f'seq_len {seq_len}')
            seq_lengths.append(seq_len)
            # global_aa_comp = features.global_aa_comp(sequence)
            # print(global_aa_comp)
            # local_aa_comp = features.local_aa_comp(sequence)

            # Global amino acid composition (i.e. percentages of all 20 amino acids present in whole sequence)
            # Local amino acid composition (i.e. over first 50 amino acids or last 50 amino acids)
            # Isoelectric point & molecular weight (e.g. see HERE or HERE)
            # Specific sequence patterns near the start or near the end of the sequence

        else:
            print(f'ID {record.id}, this {subset_name} sequence {record.seq} has invalid amino acid character. '
                  f'Therefore it is excluded from this dataset)')

    print(f'After removing invalid sequences: {len(valid_seqs)} for {subset_name} from original: {sequence_count}')
    return valid_seqs, np.array(seq_lengths)


def build_features_for_all_proteins():
    # read in each class and exclude sequences with invalid characters (probably `X`)
    cyto_valid_seqs, cyto_lengths = read_and_validate_all_fasta_and_label('cyto')
    mito_valid_seqs, mito_lengths = read_and_validate_all_fasta_and_label('mito')
    # nuclear_valid_seqs, nuclear_lengths = read_and_validate_all_fasta_and_label('nuclear')
    # other_valid_seqs, other_lengths = read_and_validate_all_fasta_and_label('other')
    # secreted_valid_seqs, secreted_lengths = read_and_validate_all_fasta_and_label('secreted')

    # Use biopython's built-in property tools for feature engineering
    for cyto_valid_seq in cyto_valid_seqs:
        features = Features(cyto_valid_seq)
        global_aa_comp = features.global_aa_comp()
        mol_weight = features.mol_weight()
        pass


if __name__ == '__main__':
    build_features_for_all_proteins()

