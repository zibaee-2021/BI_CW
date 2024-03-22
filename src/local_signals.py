"""
Some known localisation sequence motifs
"""

def _has_MTS(seq):
    """
    Alternating hydrophobic and positively-charged residue, VKARIK etc
    "A canonical mitochondrial localization signal (MLS) or mitochondrial targeting sequence (MTS) is a short peptide,
    about 15–70 amino acids long, bearing positively charged basic residues, that directs the transport of a protein
    to the mitochondria. These canonical sequences are located at the N-terminal of a given protein, consisting of an
    alternating pattern of hydrophobic and positively charged residues that form an amphipathic helix (Bolender,
    Sickmann, Wagner, Meisinger, & Pfanner, 2008; Brix, Dietmeier, & Pfanner, 1997)."
    Manipulation of mitochondrial genes and mtDNA heteroplasmy Bacmana et al. Chapter 19, Methods in Cell Biology 2020
    :param seq:
    :return:
    """
    pass


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


def _has_ER_Nterm_signal(seq):
    """
    Endoplasmic Reticulum (ER) Localization: Similar to secreted proteins, those that localize to the ER generally
    have an N-terminal signal sequence that directs the protein to the ER.
    KDEL or HDEL Sequences: For retention in the ER, proteins often end with specific sequences such as
    KDEL (Lys-Asp-Glu-Leu) in eukaryotes or HDEL in yeast.
    :return:
    """
    pass


def _has_nuclear_signal(seq):
    """
    Nuclear Localization Signals (NLS): These are short sequences rich in positively charged lysines or arginines
    that direct the protein into the nucleus through the nuclear pore complex.
    Proteins may also have Nuclear Export Signals (NES), which can direct them out of the nucleus.
    :return:
    """
    pass
