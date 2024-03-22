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
import matplotlib.pyplot as plt


def plot(x: list, y: list, x_label: str, y_label: str, title: str):
    x = list(range(len(x)))
    y = sorted(y)
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, alpha=0.6)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)
    plt.show()


if __name__ == '__main__':

    import os
    print(os.getcwd())
