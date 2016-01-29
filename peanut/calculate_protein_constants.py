"""
Calculates extinction coefficients, molecular weight, and absorbance for proteins or 
protein segments of interest based on string representations of sequences

Dependencies : None
"""

def calculate_extinction_coeff(aa_sequence, verbose=True):
    """
    Calculates the two extinction coefficients (with and without incorporating Cys) for
    a given amino acid sequence

    Arguments:
    ----------
        aa_sequence : str
            1-letter code amino acid sequence of protein or protein segment of interest
        verbose : Bool
            if verbose==True, print number of Tyr, Trp, and Cys found in sequence
    Returns:
    --------
        secreted : float
            Extinction Coefficient in cm^-1 M^-1 assuming all cysteines form cystines
        cytosolic : float
            Extinction Coefficient in cm^-1 M^-1 ignoring cystines
    """
    count_Tyr = aa_sequence.count('Y')
    count_Trp = aa_sequence.count('W')
    count_Cys = aa_sequence.count('C')
    if verbose==True:
        print(str(count_Tyr)+" Tyr")
        print(str(count_Trp)+" Trp")
        print(str(count_Cys)+" Cys")

    cytosolic = 1490.0*count_Tyr + 5500.0*count_Trp
    secreted = 1490.0*count_Tyr + 5500.0*count_Trp + 125.0*count_Cys/2.0
    return secreted, cytosolic

def calculate_molecular_weight(aa_sequence):
    """
    Calculates molecular weight of a protein based on a given amino acid sequence
    Uses dictionary of molecular weight of amino acids and assumes protein weight
    is the sum of all component amino acids minus a water molecule per peptide bond

    Arguments:
    ----------
        aa_sequence : str
            1-letter code amino acid sequence of protein or protein segment of interest
    Returns:
    --------
        molecular_weight : float 
            Molecular Weight in g/mol
    """
    # wtf do you do with B and Z?
    aa_to_mw = {'A':89.0935, 'C':121.1590, 'D':133.1032, 'E':147.1299, 'F':165.1900, 'G':75.0669, 'H':155.1552, 'I':131.1736, 'K':146.1882, 'L':131.1736, 'M':149.2124, 'N':132.1184, 'P':115.1310, 'Q':146.1451, 'R':174.2017, 'S':105.0930, 'T':119.1197, 'V':117.1469, 'W':204.2262, 'Y':181.1894}
    molecular_weight = 0.0
    for k, char in enumerate(aa_sequence):
        if not aa_to_mw.has_key(char):
            if char == '*' and k == len(aa_sequence)-1:
                continue
            if char == '\n':
                continue
            else:
                raise IOError(str(char)+" is not a recognized amino acid")
        molecular_weight += aa_to_mw[char]
        if k != 0: # water mlc given off in peptide bond
            molecular_weight -= 18.01528
    return molecular_weight

def calculate_absorbance(aa_sequence, secreted=None, molecular_weight=None):
    """
    Calculates absorbance, based on secreted extinction coefficients and molecular weight
    Will calculate extinction coefficient and / or molecular weight if not given as inputs
    
    Arguments:
    ----------
        aa_sequence : str
            1-letter code amino acid sequence of protein or protein segment of interest
        Optional:
        ---------
            secreted : float
                Extinction Coefficient in cm^-1 M^-1 assuming all cysteines form cystines
            molecular_weight : float
                Molecular Weight in g/mol
    Returns:
    --------
        absorbance : float
    """
    if secreted==None:
        secreted, cytosolic = calculate_extinction_coeff(aa_sequence, verbose=False)
    if molecular_weight==None:
        molecular_weight = calculate_molecular_weight(aa_sequence)
    # which extinction coeff do you use?
    # or both?
    absorbance = secreted / molecular_weight
    return absorbance

