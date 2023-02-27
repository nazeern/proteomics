"""
Misc. analysis functions

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012 Casey W. Stark. All Rights Reserved.

  This file is part of `MIDA`.

"""

import numpy as np
import scipy

def renormalize(molecule, distribution):
    """
    Renormalize the fractional isotopomer distribution to match how the
    experimental data is handled.

    If the molecule weighs more than 2400 amu, we normalize to the sum of the
    M0 - M4 abundances. If the molecule weighs less than 2400 amu, we normalize
    to the sum of the M0 - M3 abundances.

    """
    if molecule.base_mass < 2400:
        cut = 4
    else:
        cut = 5

    # slice the distribution from M0 to the cut isotopomer mass (M3 or M4),
    # sum it, and divide the distribution by it.
    return distribution / (distribution[:cut]).sum()

def convert_p_to_abundances(p_values, natural_abundances):
    """
    Takes an array of p values and converts them to isotopic abundance values.
    Note that this only works for elements with only two isotopes. Otherwise,
    p values are not enough information.

    """
    # grab the "p" isotope natural abundance
    na0 = natural_abundances[0]
    na1 = natural_abundances[1]

    # validations before we continue
    if natural_abundances.shape[0] > 2:
        raise Exception("p values are ill-defined for elements with more than 2 elements. Please convert to abundances manually.")

    abundances = np.zeros((len(p_values), len(natural_abundances)))
    abundances[:, 0] = na0 - p_values
    abundances[:, 1] = na1 + p_values

    if (abundances < 0.0).any() or (abundances > 1.0).any():
        raise Exception("The supplied p values created bad abundance values! Please make sure that they are positive and small enough.")

    return abundances
