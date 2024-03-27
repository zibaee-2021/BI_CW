import local_signals
import model
import dataset
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from features import Features
from Bio.Data import IUPACData


def build_31_features_for_challenge(sequence):
    features_ary = np.zeros((1, 32), dtype=np.float32)

    feats = Features(sequence)
    features_ary[0, :20] = feats.global_aa_comp()

    Nterm_50_gravy, Cterm_50_gravy = feats.local_hydrophobicity(sequence)
    features_ary[0, 21] = Nterm_50_gravy
    features_ary[0, 22] = Cterm_50_gravy

    features_ary[0, 23] = feats.mol_weight()
    features_ary[0, 24] = feats.mean_flexibility()
    features_ary[0, 25] = feats.aromaticity()
    features_ary[0, 26] = feats.isoelectric_point()
    features_ary[0, 27], features_ary[0, 28] = feats.net_charge_and_mnc()
    features_ary[0, 29], _ = feats.total_charge_and_mtc()
    features_ary[0, 30] = feats.mean_hydrophobicity()
    features_ary[0, 31] = len(sequence)

    return features_ary


def read_and_validate_all_challenge_fasta_seq_and_id():
    # sequence, id = '', ''
    id_seq = {}
    sequence_count = len(list(SeqIO.parse(f'../datasets/fasta_files/challenge.fasta', 'fasta')))
    print(f'Expecting 20 sequences in challenge.fasta. Found {sequence_count}')

    for record in SeqIO.parse(f'../datasets/fasta_files/challenge.fasta', 'fasta'):
        if set(record.seq) <= set(IUPACData.protein_letters):
            sequence = record.seq
            id = record.id
            print(f'ID {id}, \n has sequence {sequence}')
            id_seq[id] = sequence
    return id_seq


if __name__ == '__main__':

    id_seq = read_and_validate_all_challenge_fasta_seq_and_id()
    for id, seq in id_seq.items():
        challenge_feats = build_31_features_for_challenge(seq)
        print(f'id {id}')
        print(f'Labelled features has shape {challenge_feats.shape}')

