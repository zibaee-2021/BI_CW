import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import local_signals


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
        glob_aa_comp = {aa: round(1000 * aa_comp) for aa, aa_comp in aa_compo.items()}
        global_aa_comp = np.array(list(glob_aa_comp.values()))
        return global_aa_comp

    def local_aa_comp(self, sequence):
        if len(sequence) > 50:
            Nterm_50_seq, Cterm_50_seq = sequence[:50], sequence[-50:]
            Nterm_50_aa_comp, Cterm_50_aa_comp = 0, 0
            if len(sequence) >= 100:
                Nterm_50_aa_comp = self.global_aa_comp(Nterm_50_seq)
                Cterm_50_aa_comp = self.global_aa_comp(Cterm_50_seq)
            return Nterm_50_aa_comp, Cterm_50_aa_comp
        else:
            return self.global_aa_comp(sequence), self.global_aa_comp(sequence)

    def local_hydrophobicity(self, sequence):
        if len(sequence) > 50:
            Nterm_50_seq, Cterm_50_seq = sequence[:50], sequence[-50:]
            Nterm_50_aa_comp, Cterm_50_aa_comp = 0, 0
            if len(sequence) >= 100:
                Nterm_50_aa_comp = self.mean_hydrophobicity(sequence=Nterm_50_seq)
                Cterm_50_aa_comp = self.mean_hydrophobicity(sequence=Cterm_50_seq)
            return Nterm_50_aa_comp, Cterm_50_aa_comp
        else:
            return self.mean_hydrophobicity(sequence), self.mean_hydrophobicity(sequence)

    def mol_weight(self, sequence=None):
        if sequence is None:
            mw = round(self._analysed_seq.molecular_weight())
        else:
            mw = round(ProteinAnalysis(sequence).molecular_weight())
        return round(mw / 1000)  # KDa

    def aromaticity(self, sequence=None):
        if sequence is None:
            arom = self._analysed_seq.aromaticity()
        else:
            arom = ProteinAnalysis(sequence).aromaticity()
        return int(round(1000 * arom))

    def mean_flexibility(self, sequence=None):
        if sequence is None:
            flex = self._analysed_seq.flexibility()
        else:
            flex = ProteinAnalysis(sequence).flexibility()
        mean_flex = np.mean(flex)
        return int(round(mean_flex * 1000))

    def isoelectric_point(self, sequence=None):
        if sequence is None:
            iep = self._analysed_seq.isoelectric_point()
        else:
            iep = ProteinAnalysis(sequence).isoelectric_point()
        return int(round(iep * 100))

    def net_charge_and_mnc(self, sequence=None):
        if sequence is None:
            nc = self._analysed_seq.charge_at_pH(7.4)
            nc = round(100 * nc)
            mnc = round((nc / len(self.seq)) * 100)
        else:
            anal_seq = ProteinAnalysis(sequence)
            nc = anal_seq.charge_at_pH(7.4)
            nc = round(100 * nc)
            mnc = round((nc / len(sequence)) * 100)
        return nc, mnc

    def total_charge_and_mtc(self, sequence=None):
        if sequence is None:
            aa_counts = self._analysed_seq.count_amino_acids()
            tc = aa_counts['K'] + aa_counts['R'] + aa_counts['D'] + aa_counts['E']
            mtc = round((tc / len(self.seq)) * 100)
        else:
            se = ProteinAnalysis(sequence)
            aa_counts = se.count_amino_acids()
            tc = aa_counts['K'] + aa_counts['R'] + aa_counts['D'] + aa_counts['E']
            mtc = round((tc / len(sequence)) * 100)
        return tc, mtc

    def mean_hydrophobicity(self, sequence=None):
        """Calculate the GRAVY (Grand Average of Hydropathy) according to Kyte and Doolitle, 1982."""
        if sequence is None:
            return round(100 * self._analysed_seq.gravy(scale='KyteDoolitle'))
        else:
            return round(100 * ProteinAnalysis(sequence).gravy(scale='KyteDoolitle'))

    def has_mito_signal_peptide(self, sequence):
        if local_signals.has_MTS(sequence):
            return 1
        else:
            return 0

    def has_nuclear_signal_peptide(self, sequence):
        if local_signals.has_nuclear_signal(sequence):
            return 1
        else:
            return 0

    def has_secreted_signal_peptide(self, sequence):
        if local_signals.has_secreted_signal(sequence):
            return 1
        else:
            return 0


if __name__ == '__main__':
    fe = Features()
    # hydro = fe.mean_hydrophobicity(sequence='DSNQDSNQDS')
    # hydro2 = fe.mean_hydrophobicity(sequence='LLLVGILFWA')
    # hydro3 = fe.mean_hydrophobicity(sequence='DSNQDILFWA')
    # nc1, mnc1 = fe.net_charge_and_mnc(sequence='AAAAA')
    # nc2, mnc2 = fe.net_charge_and_mnc(sequence='AAAEEE')
    # nc3, mnc3 = fe.net_charge_and_mnc(sequence='AAAKKK')
    # nc4, mnc4 = fe.net_charge_and_mnc(sequence='EEEKKK')

    tc1, mtc1 = fe.total_charge_and_mtc(sequence='AAAAA')
    tc2, mtc2 = fe.total_charge_and_mtc(sequence='AAAEEE')
    tc3, mtc3 = fe.total_charge_and_mtc(sequence='AAAKKK')
    tc4, mtc4 = fe.total_charge_and_mtc(sequence='EEEKKK')

    pass