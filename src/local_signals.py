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
to_peroxisomes = 'SKL'

return_to_ER =   ['KDEL', 'RDEL', 'HDEL', 'ADEL', 'DDEL']

# "Whereas in mammals, the receptor binds to KDEL, RDEL and HDEL sequences,
# other organisms use alternative sequences,
# such as ADEL or DDEL signals (Lewis and Pelham, 1990; Pelham, 1992; Pidoux and Armstrong, 1992)."
# ----------------------------------------------------------------------------------------
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def has_nuclear_signal(seq):
    """
    Nuclear Localization Signals (NLS): These are short sequences rich in positively charged lysines or arginines
    that direct the protein into the nucleus through the nuclear pore complex.
    Proteins may also have Nuclear Export Signals (NES), which can direct them out of the nucleus.
    :return:
    """
    pattern = r'([KR]{5,6})'
    return bool(re.search(pattern, seq))


def has_MTS(seq):
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
    if len(seq) < 90: return False
    seq_Nterm80 = seq[:80]
    # template = '----R---R--K---R-----R---'
    positions_of_KR = [i for i, char in enumerate(seq_Nterm80) if char in ['K', 'R']]
    if len(positions_of_KR) < 5:
        return False
    else:
        pattern = r'([^\sKR]{2,3})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,5})(K|R)([^\sKR]{2,3})'
        if re.search(pattern, seq_Nterm80) is None:
            # print('No matching pattern')
            return False
        matches = re.finditer(pattern, seq_Nterm80)
        concatenated_segments = ''

        # print(f'matches is None {matches is None}')
        for match in matches:
            # print(f"Match: {match.group()}, Start: {match.start()}, End: {match.end()}")
            for i in range(1, 12, 2):  # Extract only the segments between K or R
                concatenated_segments += match.group(i)
            # hydr = Features(concatenated_segments).mean_hydrophobicity()
            hydr = ProteinAnalysis(concatenated_segments).gravy(scale='KyteDoolitle')
            if hydr < 0:
                return False
            else:
                return True


def has_secreted_signal(seq):
    """
    Secreted/Extracellular Localization:
    Signal Peptide: Proteins targeted for secretion typically begin with an N-terminal signal peptide that directs them to the endoplasmic reticulum (ER) for entry into the secretory pathway. This signal is usually 15-30 amino acids long and is characterized by a core of hydrophobic residues.
    Cleavage Sites: The signal peptide is cleaved off in the ER by signal peptidase.
    :return:
    """
    if seq is None: return False
    if len(seq) < 40: return False
    Nterm9 = (seq[:12]) # 12
    nc = ProteinAnalysis(Nterm9).charge_at_pH(7.4)
    if nc < 0:
        return False
    else:
        Mid15 = (seq[12:24]) # 12
        gravyMid = ProteinAnalysis(Mid15).gravy(scale='KyteDoolitle')  # positive means hydrophobic
        if gravyMid <= 1:
            return False
        else:
            C7 = (seq[24:33])  # 9
            for aa in C7:
                if aa in ['K', 'R', 'E', 'D']:
                    return False

            return True




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