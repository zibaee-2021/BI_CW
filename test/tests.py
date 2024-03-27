from unittest import TestCase
from BI_CW.src import main
import numpy as np
from numpy import testing as npt
from BI_CW.src.features import Features
from BI_CW.src import local_signals

class TestMain(TestCase):

    # def setUp(self):

    def test__has_return_to_ER_Cterm_signal(self):
        true_er_sig = ['XXXKDEL', 'XXXRDEL', 'XXXHDEL', 'XXXADEL', 'XXXDDEL']
        expected = True
        actual_0 = local_signals._has_return_to_ER_Cterm_signal(true_er_sig[0])
        self.assertEqual(expected, actual_0)
        actual_1 = local_signals._has_return_to_ER_Cterm_signal(true_er_sig[1])
        self.assertEqual(expected, actual_1)
        actual_2 = local_signals._has_return_to_ER_Cterm_signal(true_er_sig[2])
        self.assertEqual(expected, actual_2)
        actual_3 = local_signals._has_return_to_ER_Cterm_signal(true_er_sig[3])
        self.assertEqual(expected, actual_3)
        actual_4 = local_signals._has_return_to_ER_Cterm_signal(true_er_sig[4])
        self.assertEqual(expected, actual_4)

        expected = False
        actual = local_signals._has_return_to_ER_Cterm_signal('KDELX')
        self.assertEqual(expected, actual)
        expected = False
        actual = local_signals._has_return_to_ER_Cterm_signal('KXDEL')
        self.assertEqual(expected, actual)
        expected = False
        actual = local_signals._has_return_to_ER_Cterm_signal('EDEL')
        self.assertEqual(expected, actual)

    def test_global_aa_comp(self):
        seq = 'CACKDKK'
        features = Features(seq)
        expected = {'A': 143, 'C': 286, 'D': 143, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0,
                    'K': 429, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0,
                    'V': 0, 'W': 0, 'Y': 0}
        expected = np.array(list(expected.values()))
        actual = features.global_aa_comp(seq)
        npt.assert_equal(expected, actual)

        # npt.assert_array_equal(expected, actual)

    def test_local_aa_comp(self):
        seq = 'AAAAAAAAAAAAAAAAAAAA' \
              'AAAAAAAAAAAAAAAAAAAA' \
              'AAAAAAAAAAYYYYYYYYYY' \
              'YYYYYYYYYYYYYYYYYYYY' \
              'YYYYYYYYYYYYYYYYYYYY'
        expected_Nterm = {'A': 1000, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0,
                    'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0,
                    'V': 0, 'W': 0, 'Y': 0}
        expected_Cterm = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0,
                          'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0,
                          'V': 0, 'W': 0, 'Y': 1000}
        expected = np.array(list(expected_Nterm.values())), np.array(list(expected_Cterm.values()))
        features = Features()
        actual = features.local_aa_comp(seq)
        npt.assert_equal(expected, actual)

    def test_mol_weight(self):
        my_seq = 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQ'\
                 'LSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAKN'\
                 'WFVGLKKNGSCKRGPRTHYGQKAILFLPLPV'
        features = Features(my_seq)
        expected = 16.974
        actual = features.mol_weight()
        self.assertEqual(expected, actual)

    def test__has_ER_Nterm_signal(self):
        expected = True
        actual = local_signals._has_ER_Nterm_signal('KMSFVSLLLVGILFWATEAEQLTKCEVFQ')
        self.assertEqual(expected, actual)
        pass

    def test__has_MTS(self):
        expected = True
        actual = local_signals.has_MTS('MLSLRQSIRFFKPATRTLCSSRYLL')
        self.assertEqual(expected, actual)

        expected = False
        actual = local_signals.has_MTS('MKKKRQSIRFFKPATRTLCSSRYLL')
        self.assertEqual(expected, actual)
        expected = False
        actual = local_signals.has_MTS('QQQLRQSIRFFKPATRTLCSSRQQQ')
        self.assertEqual(expected, actual)
        expected = False
        actual = local_signals.has_MTS('MLSLEQSIRFFKPATRTLCSSRYLL')
        self.assertEqual(expected, actual)

    def test_has_secreted_signal(self):
        expected = True
        actual = local_signals.has_secreted_signal('MDSKGSSQKGSR'
                                                   'LLLLLVVSNLLL'
                                                   'CQGVVS')
        self.assertEqual(expected, actual)

    def test_has_nuclear_signal(self):
        expected = True
        actual = local_signals.has_nuclear_signal('KKKKK')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('KKKKKK')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('RRRRR')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('RRRRRR')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('KKKRRR')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('RRRKKK')
        self.assertEqual(expected, actual)

        expected = False
        actual = local_signals.has_nuclear_signal('RRRR')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('KKKK')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('AKAKAK')
        self.assertEqual(expected, actual)
        actual = local_signals.has_nuclear_signal('AKAKAK')
        self.assertEqual(expected, actual)

