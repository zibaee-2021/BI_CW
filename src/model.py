from time import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score, \
    balanced_accuracy_score, roc_auc_score, precision_score, recall_score, f1_score


def split_train_test(ds):
    y = ds[:, 0]
    X = ds[:, 1:]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)
    return X_train, X_test, y_train, y_test


def train_RFC(X_train, y_train):

    clf = RandomForestClassifier()
    start = time()
    clf.fit(X_train, y_train)
    importances = clf.feature_importances_
    print(f'Trained in {round((time() - start))} secs')
    preds = clf.predict(X_train)
    acc = round(accuracy_score(y_train, preds) * 100, 1)
    # print(f'Train accuracy is {acc}')
    return clf, importances


def test_RFC(clf, X_test, y_test):
    preds = clf.predict(X_test)
    acc = round(accuracy_score(y_test, preds) * 100, 1)
    print(f'Test accuracy is {acc}')
    return acc


if __name__ == '__main__':

    data = np.genfromtxt('../datasets/fiveClass.csv', delimiter=',')
    print(f'Labelled features has shape {data.shape}')

    # Remove the 40 features of Nterm50 aa comp and Cterm50 aa comp. Remove the signal peptides too.
    data_less = np.concatenate([data[:, :21], data[:, 61:-3]], axis=1)
    # X_train, X_test, y_train, y_test = split_train_test(data)
    X_train, X_test, y_train, y_test = split_train_test(data_less)

    # print(f'Expecting train to be {int(0.8 * 9401)} rows by 73 columns: {X_train.shape}')
    print(f'Expecting TRAIN to be {int(0.8 * 9401)} rows by 32 columns: {X_train.shape}')
    clf, importances = train_RFC(X_train, y_train)

    std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
    feature_names = [f'{i}' for i in range(1, 33)]
    forest_importances = pd.Series(importances, index=feature_names)
    fig, ax = plt.subplots()
    forest_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title('Feature importances using MDI')
    ax.set_ylabel('Mean decrease in impurity')
    ax.set_xlabel('features', fontsize=10)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    plt.xticks(fontsize=10)
    plt.show()

    test_RFC(clf, X_test, y_test)
    # 'All 5 classes: Train accuracy is 100, Test accuracy is 62.8

    print('\nRun on 2 classes, prokaryotic vs eukarytoic')
    data = np.genfromtxt('../datasets/twoClass.csv', delimiter=',')
    print(f'labelled feats has shape {data.shape}')
    X_train, X_test, y_train, y_test = split_train_test(data)
    clf, importances = train_RFC(X_train, y_train)
    test_RFC(clf, X_test, y_test)
    # 'other' (prokayotic proteins) vs the rest: Train accuracy is 100.0, Test accuracy is 92.5

    print('\nRun on 4 classes, prokaryotic removed')
    data = np.genfromtxt('../datasets/noProk.csv', delimiter=',')
    print(f'labelled feats has shape {data.shape}')
    X_train, X_test, y_train, y_test = split_train_test(data)
    clf, importances = train_RFC(X_train, y_train)
    test_RFC(clf, X_test, y_test)
    # '4 eukaryotic only: Train accuracy is 100, Test accuracy is 61


