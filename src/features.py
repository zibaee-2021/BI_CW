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
            self.seq = seq

    def global_aa_comp(self, sequence=None):
        if sequence is None:
            aa_compo = self._analysed_seq.get_amino_acids_percent()
        else:
            aa_compo = ProteinAnalysis(sequence).get_amino_acids_percent()

        return {aa: round((100 * aa_comp), 1) for aa, aa_comp in aa_compo.items()}

    def local_aa_comp(self, sequence):
        Nterm_50_aa_comp, Cterm_50_aa_comp = 0, 0

        if len(sequence) >= 100:
            Nterm_50_aa_comp = self.global_aa_comp(sequence[:50])
            Cterm_50_aa_comp = self.global_aa_comp(sequence[-50:])

        return Nterm_50_aa_comp, Cterm_50_aa_comp

    def mol_weight(self, sequence=None):
        if sequence is None:
            return round(self._analysed_seq.molecular_weight())
        else:
            return round(ProteinAnalysis(sequence).molecular_weight())

    def aromaticity(self, sequence=None):
        if sequence is None:
            return self._analysed_seq.aromaticity()
        else:
            return ProteinAnalysis(sequence).aromaticity()

    def flexibility(self, sequence=None):
        if sequence is None:
            return self._analysed_seq.flexibility()
        else:
            return ProteinAnalysis(sequence).flexibility()

    def isoelectric_point(self, sequence=None):
        if sequence is None:
            return self._analysed_seq.isoelectric_point()
        else:
            return ProteinAnalysis(sequence).isoelectric_point()

    def net_charge_and_mnc(self, sequence=None):
        if sequence is None:
            net_charge = self._analysed_seq.charge_at_pH(7.4)
            return net_charge, net_charge / len(self.seq)
        else:
            anal_seq = ProteinAnalysis(sequence)
            net_charge = anal_seq.charge_at_pH(7.4)
            return net_charge, net_charge / len(sequence)

    def total_charge_and_mtc(self, sequence=None):
        if sequence is None:
            aa_counts = self._analysed_seq.count_amino_acids()
            tc = aa_counts['K'] + aa_counts['R'] + aa_counts['D'] + aa_counts['E']
            return tc, tc / len(self.seq)
        else:
            se = ProteinAnalysis(sequence)
            aa_counts = se.count_amino_acids()
            tc = aa_counts['K'] + aa_counts['R'] + aa_counts['D'] + aa_counts['E']
            return tc, tc / len(sequence)

    def mean_hydrophobicity(self, sequence=None):
        """Calculate the GRAVY (Grand Average of Hydropathy) according to Kyte and Doolitle, 1982."""
        if sequence is None:
            return round(100 * self._analysed_seq.gravy(scale='KyteDoolitle'))
        else:
            return round(100 * ProteinAnalysis(sequence).gravy(scale='KyteDoolitle'))


if __name__ == '__main__':
    fe = Features()
    hydro = fe.mean_hydrophobicity(sequence='DSNQDSNQDS')
    hydro2 = fe.mean_hydrophobicity(sequence='LLLVGILFWA')
    hydro3 = fe.mean_hydrophobicity(sequence='DSNQDILFWA')
    pass