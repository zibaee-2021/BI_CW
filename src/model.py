from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score, balanced_accuracy_score, \
    roc_auc_score, precision_score, recall_score, f1_score
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from time import time


def split_train_test(ds):
    y = ds[:, 0]
    X = ds[:, 1:]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)
    return X_train, X_test, y_train, y_test



def train_RFC(X_train, y_train):

    clf = RandomForestClassifier()
    start = time()
    clf.fit(X_train, y_train)
    print(f'Trained in {round((time() - start))} secs')
    preds = clf.predict(X_train)
    acc = round(accuracy_score(y_train, preds) * 100, 1)
    print(f'Train accuracy is {acc}')
    return clf


def test_RFC(clf, X_test, y_test):
    preds = clf.predict(X_test)
    acc = round(accuracy_score(y_test, preds) * 100, 1)
    print(f'Test accuracy is {acc}')


if __name__ == '__main__':
    data = np.genfromtxt('../datasets/all_labelled_feats.csv', delimiter=',')
    print(f'labelled feats has shape {data.shape}')
    X_train, X_test, y_train, y_test = split_train_test(data)

    print(f'Expecting train to be {int(0.8 * 9401)} rows by 24 columns: {X_train.shape}')
    clf = train_RFC(X_train, y_train)
    test_RFC(clf, X_test, y_test)

