from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score, balanced_accuracy_score, \
    roc_auc_score, precision_score, recall_score, f1_score
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from time import time
import pandas as pd
import matplotlib.pyplot as plt
import local_signals, dataset

cyto, mito, secreted, nuclear, other = 'cyto', 'mito', 'secreted', 'nuclear', 'other'
labels = {cyto: 0, secreted: 1, mito: 2, nuclear: 3, other: 4}

def split_train_test(ds):
    y = ds[:, 0]
    X = ds[:, 1:]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)
    return X_train, X_test, y_train, y_test


def is_mitochondrial(ds):
    noMito = ds[np.isin(ds[:, 0], [0, 1, 3, 4])]
    local_signals.has_MTS(noMito)


if __name__ == '__main__':
    cyto, mito, secreted, nuclear, other = 'cyto', 'mito', 'secreted', 'nuclear', 'other'
    labels = {cyto: 0, secreted: 1, mito: 2, nuclear: 3, other: 4}

    cyto_valid_seqs, _ = dataset.read_and_validate_all_fasta_and_label(cyto)
    secreted_valid_seqs, _ = dataset.read_and_validate_all_fasta_and_label(secreted)
    mito_valid_seqs, _ = dataset.read_and_validate_all_fasta_and_label(mito)
    nuclear_valid_seqs, _ = dataset.read_and_validate_all_fasta_and_label(nuclear)
    other_valid_seqs, _ = dataset.read_and_validate_all_fasta_and_label(other)

    #############################
    # % MITO ACCURACY
    #############################
    prots = len(mito_valid_seqs)
    found = 0
    for seq in mito_valid_seqs:
        if local_signals.has_MTS(str(seq)):
            found += 1
    print(f'\nMITO ACC: is {found} out of {prots} mitochondrial proteins. Accuracy {round(found/prots * 100, 2)}%. '
          f'So False negatives are percentage of unfound {round((prots-found)/prots * 100, 2)}%.')
    # Acc is 41 out of 1020 mitochondrial proteins. Accuracy 4.02%. So False negatives are percentage of unfound 1015.98%.

    # FALSE POS for MITO signal, (any positives in 'other' are FPs.)
    prots = len(other_valid_seqs)
    fps = 0
    for seq in other_valid_seqs:
        if local_signals.has_MTS(str(seq)):
            fps += 1
    print(f'MITO False positives {fps} out of {prots} prokayotic proteins. {round(fps/prots * 100, 2)}%')
    # False positives .. out of ... prokayotic proteins. ...%

    # POSSIBLY FALSE POS for MITO signal, (any positives in the other 3 eukaryotic proteins are possibly FP.)
    euk3_others = secreted_valid_seqs + cyto_valid_seqs + nuclear_valid_seqs
    prots = len(euk3_others)
    fps = 0
    for seq in euk3_others:
        if local_signals.has_MTS(str(seq)):
            fps += 1

    print(f'MITO False positives {fps} out of {prots} nuclear, secretory & cytoplasmic proteins. '
          f'{round(fps/prots * 100, 2)}%')
    # False positives .. out of .. nuclear, secretory & cytoplasmic proteins. ..%

    #############################
    # % NUCLEAR ACCURACY
    #############################
    prots = len(nuclear_valid_seqs)
    found = 0
    for seq in nuclear_valid_seqs:
        if local_signals.has_nuclear_signal(str(seq)):
            found += 1
    print(f'\nNUCLEAR ACC is {found} out of {prots} nuclear proteins. Accuracy {round(found/prots * 100, 2)}%. '
          f'So False negatives are percentage of unfound {round((prots-found)/prots * 100, 2)}%.')

    # FALSE POS for NUCLEAR signal in prokaryotes, (any positives in 'other' are FPs.)
    prots = len(other_valid_seqs)
    fps = 0
    for seq in other_valid_seqs:
        if local_signals.has_nuclear_signal(str(seq)):
            fps += 1
    print(f'NUCLEAR False positives {fps} out of {prots} prokaryotic proteins. {round(fps/prots * 100, 2)}%')
    # False positives .. out of ... prokaryotic proteins. .. %

    # FALSE POS for NUCLEAR signal, (any positives in the other 3 eukaryotic proteins are possibly FP)
    euk3_others = secreted_valid_seqs + cyto_valid_seqs + mito_valid_seqs
    prots = len(euk3_others)
    fps = 0
    for seq in euk3_others:
        if local_signals.has_nuclear_signal(str(seq)):
            fps += 1

    print(f'NUCLEAR False positives {fps} out of {prots} mito, secretory & cytoplasmic proteins. '
          f'{round(fps/prots * 100, 2)}%')
    # False positives ... out of ... mito, secretory & cytoplasmic proteins. ... %

    #############################
    # % SECRETORY ACCURACY
    #############################
    prots = len(secreted_valid_seqs)
    found = 0
    for seq in secreted_valid_seqs:
        if local_signals.has_secreted_signal(str(seq)):
            found += 1
    print(f'SECRETORY ACC is {found} out of {prots} secreted proteins. Accuracy {round(found/prots * 100, 2)}%. '
          f'So False negatives are percentage of unfound {round((prots-found)/prots * 100, 2)}%.')

    # FALSE POS for SECRETORY signal, (any positives in the other 3 eukaryotic proteins are possibly FP)
    euk3_others = nuclear_valid_seqs + cyto_valid_seqs + mito_valid_seqs
    prots = len(euk3_others)
    fps = 0
    for seq in euk3_others:
        if local_signals.has_secreted_signal(str(seq)):
            fps += 1

    print(f'SECRETORY False positives {fps} out of {prots} mito, nuclear & cytoplasmic proteins. '
          f'{round(fps/prots * 100, 2)}%')
    # False positives ... out of ... mito, nuclear & cytoplasmic proteins. ... %

