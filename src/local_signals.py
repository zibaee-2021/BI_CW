"""
Some known localisation sequence motifs
"""
# Anywhere in sequence:
to_nuc =   'PPKKKRKV'
to_nuc_B = '--KKKRK-'
# I think it may also be PPRRRKRV and V is not considered to be an important hydrophobic aa
from_nuc   = 'MEELSQALASSF'
from_nuc_B = 'M--L---L---F'
# ----------------------------------------------------------------------------------------
# N-terminal only:
to_mito   = 'MLSLRQSIRFFKPATRTLCSSRYLL'
to_mito_B = '----R---R--K---R-----R---'

to_plastid   = 'MVAMAMASLQSSMSSLSLSSNSFLGQPLSPITLSPFLQG'  # present in all living things!
to_plastid_B = '-------S--SS-SS-S-SS-S------S--T-S-----'

to_ER =   'MMSFVSLLLVGILFWATEAEQLTKCEVFQ'
to_ER_B = '------LLLVGILFWA-------K-E---'
# ----------------------------------------------------------------------------------------
# C-terminal only:
to_peroxisomes = 'SKL'  # ABSENT IN PROKARYOTES

return_to_ER =   ['KDEL', 'RDEL', 'HDEL', 'ADEL', 'DDEL']

# "Whereas in mammals, the receptor binds to KDEL, RDEL and HDEL sequences,
# other organisms use alternative sequences,
# such as ADEL or DDEL signals (Lewis and Pelham, 1990; Pelham, 1992; Pidoux and Armstrong, 1992)."
# ----------------------------------------------------------------------------------------
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from BI_CW.src.features import Features


def _has_return_to_ER_Cterm_signal(seq):
    cterm_peptide = seq[-4:]
    assert len(cterm_peptide) == 4

    pattern = re.compile('[KRHAD]DEL')

    if pattern.match(cterm_peptide):
        print(f'{cterm_peptide} matches the pattern')
        return True
    else:
        print(f'{cterm_peptide} does not match the pattern')
        return False


# def _has_ER_Nterm_signal(seq):
#     """
#     Endoplasmic Reticulum (ER) Localization: Similar to secreted proteins, those that localize to the ER generally
#     have an N-terminal signal sequence that directs the protein to the ER.
#     N-term, 16-30 amino acid. 3 parts: Nterm +vely charged, Hydrophobic middle region, polar uncharged C-region.
#     to_ER =   'MMSFVSLLLVGILFWATEAEQLTKCEVFQ'
#     to_ER_B = '------LLLVGILFWA-------K-E---'
#     :return:
#     """
#     Nterm6 = seq[:6]
#     nc, _ = Features(Nterm6).net_charge_and_mnc()
#     if nc < 0:
#         return False
#     else:
#         mid = seq[6:16]
#         hydroph = Features(mid).mean_hydrophobicity()
#         if hydroph < 0:
#             return False
#         else:
#             end = seq[16:30]
#             nc, _ = Features(end).net_charge_and_mnc()
#             if nc != 0 and Features(end).mean_hydrophobicity() < 0:
#                 return False
#             else:
#                 return True


def _has_MTS(seq):
    """
    Alternating hydrophobic and positively-charged residue, VKARIK etc
    "A canonical mitochondrial localization signal (MLS) or mitochondrial targeting sequence (MTS) is a short peptide,
    about 15–70 amino acids long, bearing positively charged basic residues, that directs the transport of a protein
    to the mitochondria. These canonical sequences are located at the N-terminal of a given protein, consisting of an
    alternating pattern of hydrophobic and positively charged residues that form an amphipathic helix (Bolender,
    Sickmann, Wagner, Meisinger, & Pfanner, 2008; Brix, Dietmeier, & Pfanner, 1997)."
    Manipulation of mitochondrial genes and mtDNA heteroplasmy Bacmana et al. Chapter 19, Methods in Cell Biology 2020
    :param seq: Full sequence of protein.
    :return:
    """
    seq_Nterm80 = seq[:80]
    # template = '----R---R--K---R-----R---'
    positions_of_KR = [i for i, char in enumerate(seq_Nterm80) if char in ['K', 'R']]
    if len(positions_of_KR) < 5:
        return False
    else:
        pattern = r'([^\sKR]{2,3})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,3})'
        if re.search(pattern, seq_Nterm80) is None:
            print('No matching pattern')
            return False
        matches = re.finditer(pattern, seq_Nterm80)
        concatenated_segments = ''

        print(f'matches is None {matches is None}')
        for match in matches:
            print(f"Match: {match.group()}, Start: {match.start()}, End: {match.end()}")
            for i in range(1, 12, 2):  # Extract only the segments between K or R
                concatenated_segments += match.group(i)
            hydr = Features(concatenated_segments).mean_hydrophobicity()
            if hydr < 0:
                return False
            else:
                return True


def _has_other_signals(seq):
    """

    :param seq:
    :return:
    """
    def _has_peroxisomal_signal(seq):
        """
        Peroxisomal Localization: Peroxisomal Targeting Signals (PTS): These are sequences that target proteins to
        peroxisomes, typically at the C-terminus (PTS1) or N-terminus (PTS2).
        """
        pass


def _has_secreted_signal():
    """
    Secreted/Extracellular Localization:
    Signal Peptide: Proteins targeted for secretion typically begin with an N-terminal signal peptide that directs them to the endoplasmic reticulum (ER) for entry into the secretory pathway. This signal is usually 15-30 amino acids long and is characterized by a core of hydrophobic residues.
    Cleavage Sites: The signal peptide is cleaved off in the ER by signal peptidase.
    :return:
    """


def _has_nuclear_signal(seq):
    """
    Nuclear Localization Signals (NLS): These are short sequences rich in positively charged lysines or arginines
    that direct the protein into the nucleus through the nuclear pore complex.
    Proteins may also have Nuclear Export Signals (NES), which can direct them out of the nucleus.
    :return:
    """
    pass
