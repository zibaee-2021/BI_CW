import numpy as np
from Bio import SeqIO
from Bio.Data import IUPACData
import matplotlib.pyplot as plt
from features import Features

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
            seq_lengths.append(seq_len)
        else:
            print(f'ID {record.id}, this {subset_name} sequence {record.seq} has invalid amino acid character. '
                  f'Therefore it is excluded from this dataset)')
    print(f'After removing invalid sequences, there are {len(valid_seqs)} {subset_name} proteins from a total of '
          f' {sequence_count} proteins.')
    return valid_seqs, seq_lengths


def build_features_for_all_proteins(cyto_valid_seqs):
    features_ary = np.zeros((len(cyto_valid_seqs), 69), dtype=np.float32)

    for i, cyto_valid_seq in enumerate(cyto_valid_seqs):

        feats = Features(cyto_valid_seq)

        features_ary[i, :20] = feats.global_aa_comp()

        Nterm_50_aa_comp, Cterm_50_aa_comp = feats.local_aa_comp(cyto_valid_seq)
        features_ary[i, 20:40] = Nterm_50_aa_comp
        features_ary[i, 40:60] = Cterm_50_aa_comp

        features_ary[i, 60] = feats.mol_weight()
        features_ary[i, 61] = feats.mean_flexibility()
        features_ary[i, 62] = feats.aromaticity()
        features_ary[i, 63] = feats.isoelectric_point()
        features_ary[i, 64], features_ary[i, 65] = feats.net_charge_and_mnc()
        features_ary[i, 66], features_ary[i, 67] = feats.total_charge_and_mtc()
        features_ary[i, 68] = feats.mean_hydrophobicity()

    return features_ary


def add_label(feats_ary, label):
    n = feats_ary.shape[0]
    labels_col = np.full((n, 1), label)
    labeled_feats = np.concatenate((labels_col, feats_ary), axis=1)
    return labeled_feats


if __name__ == '__main__':

    from time import time
    start = time()
    # load = True
    load = False

    if load:

        data = np.genfromtxt('../datasets/all_labelled_feats.csv', delimiter=',')
        print(f'labelled feats has shape {data.shape}')

    else:

        cyto, mito, secreted, nuclear, other = 'cyto', 'mito', 'secreted', 'nuclear', 'other'
        labels = {cyto: 0, secreted: 1, mito: 2, nuclear: 3, other: 4}

        cyto_valid_seqs, cyto_lengths = read_and_validate_all_fasta_and_label(cyto)
        cyto_feats = build_features_for_all_proteins(cyto_valid_seqs)
        cyto_feats = add_label(cyto_feats, labels[cyto])
        print(f'Completed building labelled cyto features array')

        sec_valid_seqs, sec_lengths = read_and_validate_all_fasta_and_label(secreted)
        sec_feats = build_features_for_all_proteins(sec_valid_seqs)
        sec_feats = add_label(sec_feats, labels[secreted])
        print(f'Completed building labelled {secreted} features array')

        mito_valid_seqs, mito_lengths = read_and_validate_all_fasta_and_label(mito)
        mito_feats = build_features_for_all_proteins(mito_valid_seqs)
        mito_feats = add_label(mito_feats, labels[mito])
        print(f'Completed building labelled {mito} features array')

        nuc_valid_seqs, nuc_lengths = read_and_validate_all_fasta_and_label(nuclear)
        nuc_feats = build_features_for_all_proteins(nuc_valid_seqs)
        nuc_feats = add_label(nuc_feats, labels[nuclear])
        print(f'Completed building labelled {nuclear} features array')

        other_valid_seqs, other_lengths = read_and_validate_all_fasta_and_label(other)
        other_feats = build_features_for_all_proteins(other_valid_seqs)
        other_feats = add_label(other_feats, labels[other])
        print(f'Completed building labelled {other} features array')

        all_labelled_feats = np.concatenate((cyto_feats, mito_feats, sec_feats, nuc_feats, other_feats), axis=0)

        np.savetxt('../datasets/all_labelled_feats.csv', all_labelled_feats, delimiter=",", fmt='%f')

    print(f'Completed in {round(time() - start)} seconds')
    