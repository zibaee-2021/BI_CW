from unittest import TestCase
from src import main
import numpy as np
from numpy import testing as npt
from src.features import Features


class TestMain(TestCase):

    # def setUp(self):

    def test_global_aa_comp(self):
        seq = 'CACKDKK'
        features = Features(seq)
        expected = {'A': 14.3, 'C': 28.6, 'D': 14.3, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0,
                    'K': 42.9, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0,
                    'V': 0, 'W': 0, 'Y': 0}
        actual = features.global_aa_comp(seq)
        self.assertEqual(expected, actual)

        # npt.assert_array_equal(expected, actual)

    def test_local_aa_comp(self):
        seq = 'AAAAAAAAAAAAAAAAAAAA' \
              'AAAAAAAAAAAAAAAAAAAA' \
              'AAAAAAAAAAYYYYYYYYYY' \
              'YYYYYYYYYYYYYYYYYYYY' \
              'YYYYYYYYYYYYYYYYYYYY'
        expected_Nterm = {'A': 100.0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0,
                    'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0,
                    'V': 0, 'W': 0, 'Y': 0}
        expected_Cterm = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0,
                          'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0,
                          'V': 0, 'W': 0, 'Y': 100.0}
        expected = (expected_Nterm, expected_Cterm)
        features = Features()
        actual = features.local_aa_comp(seq)
        self.assertEqual(expected, actual)

    def test_mol_weight(self):
        my_seq = 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQ'\
                 'LSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAKN'\
                 'WFVGLKKNGSCKRGPRTHYGQKAILFLPLPV'
        features = Features(my_seq)
        expected = 16974
        actual = features.mol_weight()
        self.assertEqual(expected, actual)

    # def test_