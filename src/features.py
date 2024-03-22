from Bio import SeqIO
import numpy as np
from Bio.Data import IUPACData
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParam


class Features:

    def __init__(self, seq=None):
        if seq is not None:
            self._analysed_seq = ProteinAnalysis(seq)

    def global_aa_comp(self, seq=None):
        if seq is None:
            aa_compo = self._analysed_seq.get_amino_acids_percent()
        else:
            analysed_seq_local = ProteinAnalysis(seq)
            aa_compo = analysed_seq_local.get_amino_acids_percent()
        return {aa: round((100 * aa_comp), 1) for aa, aa_comp in aa_compo.items()}

    def local_aa_comp(self, sequence):
        Nterm_50_aa_comp, Cterm_50_aa_comp = 0, 0

        if len(sequence) >= 100:
            Nterm_50_aa_comp = self.global_aa_comp(sequence[:50])
            Cterm_50_aa_comp = self.global_aa_comp(sequence[-50:])

        return Nterm_50_aa_comp, Cterm_50_aa_comp

    def mol_weight(self):
        return round(self._analysed_seq.molecular_weight())

    def aromaticity(self):
        return self._analysed_seq.aromaticity()

    def flexibility(self):
        return self._analysed_seq.flexibility()

    def isoelectric_point(self):
        return self._analysed_seq.isoelectric_point()

