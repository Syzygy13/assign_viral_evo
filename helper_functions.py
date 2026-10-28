import numpy as np

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    compared = 0
    different = 0
    allowed_bases = "ACTG"
    # Find number of differing bases is it is not a gap or ambiguous
    for ref, dist in zip(reference_sequence.upper(), distant_sequence.upper()):
        if ref not in allowed_bases or dist not in allowed_bases:
            continue
        compared += 1
        if ref != dist:
            different += 1
    if compared == 0:
        return 0.0

    # Error check to avoid ln of negative value
    percent_diff = different / compared
    if percent_diff >= 0.75:
        return float('inf')

    # d = (3/4) * ln(1 - (4/3)p)
    distance = -0.75 * np.log(1 - (4 / 3) * percent_diff)
    return float(distance * compared)

def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
