"""
The default isotope and amino acid data passed along to Molecule objects.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012, Casey W. Stark. All Rights Reserved.

  This file is part of `MIDA`.

"""

import numpy as np

from mida.data_types import ed_dtype, id_dtype, composition_dtype, \
    labile_dtype, aa_enrichment_dtype
from mida.molecule_objects import AminoAcid
from mida.data_containers import ChemicalDataContainer

###
# The element data
###

# format is the symbol, atomic number, and the number of (stable) isotopes.
element_data = np.array(
    [("H",  1, 2),
     ("C",  6, 2),
     ("N",  7, 2),
     ("O",  8, 3),
     ("S", 16, 4)],
    dtype=ed_dtype)

###
# The isotope data
###

# format is element index, A, mi (A - A0), mass, natural abundance
isotope_data = np.array(
    [(0,  1, 0,  1.007825, 0.999844),
     (0,  2, 1,  2.014101, 0.000156),
     (1, 12, 0, 12.000000, 0.9891  ),
     (1, 13, 1, 13.003355, 0.0109  ),
     (2, 14, 0, 14.003074, 0.99635 ),
     (2, 15, 1, 15.000108, 0.00365 ),
     (3, 16, 0, 15.994915, 0.99759 ),
     (3, 17, 1, 16.999132, 0.00037 ),
     (3, 18, 2, 17.999161, 0.00204 ),
     (4, 32, 0, 31.972071, 0.9493  ),
     (4, 33, 1, 32.971459, 0.0076  ),
     (4, 34, 2, 33.967867, 0.0429  ),
     (4, 36, 4, 35.967081, 0.0002  )],
    dtype=id_dtype)

###
# The amino acid data
###

amino_acid_data = {
    # Alanine
    "A": AminoAcid("A", np.array([7, 3, 1, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 4.00)], dtype=labile_dtype)),
    # Arginine
    "R": AminoAcid("R", np.array([14, 6, 4, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 3.43)], dtype=labile_dtype)),
    # Asparagine
    "N": AminoAcid("N", np.array([8, 4, 2, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 1.89)], dtype=labile_dtype)),
    # Aspartic acid
    "D": AminoAcid("D", np.array([7, 4, 1, 4, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 1.89)], dtype=labile_dtype)),
    # Cysteine (processed)
    "C": AminoAcid("C", np.array([10, 5, 2, 3, 1], dtype=composition_dtype),
                   labiles=np.array([(0, 1.62)], dtype=labile_dtype)),
    # Glutamic acid
    "E": AminoAcid("E", np.array([9, 5, 1, 4, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 3.95)], dtype=labile_dtype)),
    # Glutamine
    "Q": AminoAcid("Q", np.array([10, 5, 2, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 3.95)], dtype=labile_dtype)),
    # Glycine
    "G": AminoAcid("G", np.array([5, 2, 1, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 2.06)], dtype=labile_dtype)),
    # Histidine
    "H": AminoAcid("H", np.array([9, 6, 3, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 2.88)], dtype=labile_dtype)),
    # Isoleucine
    "I": AminoAcid("I", np.array([13, 6, 1, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 1.00)], dtype=labile_dtype)),
    # Leucine
    "L": AminoAcid("L", np.array([13, 6, 1, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.60)], dtype=labile_dtype)),
                   #aa_enrichments=np.array([(1, 4)], dtype=aa_enrichment_dtype)),
    # Lysine
    "K": AminoAcid("K", np.array([14, 6, 2, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.54)], dtype=labile_dtype)),
    # Methionine
    "M": AminoAcid("M", np.array([11, 5, 1, 2, 1], dtype=composition_dtype),
                   labiles=np.array([(0, 1.12)], dtype=labile_dtype)),
    # Phenylalanine
    "F": AminoAcid("F", np.array([11, 9, 1, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.32)], dtype=labile_dtype)),
    # Proline
    "P": AminoAcid("P", np.array([9, 5, 1, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 2.59)], dtype=labile_dtype)),
    # Serine
    "S": AminoAcid("S", np.array([7, 3, 1, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 2.61)], dtype=labile_dtype)),
    # Threonine
    "T": AminoAcid("T", np.array([9, 4, 1, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.20)], dtype=labile_dtype)),
    # Tryptophan
    "W": AminoAcid("W", np.array([12, 11, 2, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.08)], dtype=labile_dtype)),
    # Tyrosine
    "Y": AminoAcid("Y", np.array([11, 9, 1, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.42)], dtype=labile_dtype)),
    # Valine
    "V": AminoAcid("V", np.array([11, 5, 1, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.56)], dtype=labile_dtype)),
    # Methionine (processed)
    "m": AminoAcid("m", np.array([11, 5, 1, 3, 1], dtype=composition_dtype),
                   labiles=np.array([(0, 1.12)], dtype=labile_dtype)),
    # Glutamine (processed)
    "q": AminoAcid("q", np.array([7, 5, 1, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 3.95)], dtype=labile_dtype)),
    # Lysine (processed)
    "k": AminoAcid("k", np.array([16, 8, 2, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.54)], dtype=labile_dtype)),
    # Proline (processed)
    "p": AminoAcid("p", np.array([9, 5, 1, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 2.59)], dtype=labile_dtype)),
    # Asparagine (Deamidated)
    "n": AminoAcid("n", np.array([7, 4, 1, 4, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 1.89)], dtype=labile_dtype)),
    # Lysine (Guanidination)
    "k": AminoAcid("k", np.array([16, 7, 4, 2, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.54)], dtype=labile_dtype)),
    # Lysine (Carbamylated lysin)
    "k": AminoAcid("k", np.array([15, 7, 3, 3, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.54)], dtype=labile_dtype)),
    # Lysine (Ubiquitination)
    "k": AminoAcid("k", np.array([20, 10, 4, 4, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.54)], dtype=labile_dtype)),
    # Cysteine (S-Nitrosylation)
    "c": AminoAcid("c", np.array([6, 3, 2, 3, 1], dtype=composition_dtype),
                   labiles=np.array([(0, 1.62)], dtype=labile_dtype)),
   # Tyrosine (Photo-Decomposition)
    "y": AminoAcid("y", np.array([10, 9, 2, 4, 0], dtype=composition_dtype),
                   labiles=np.array([(0, 0.42)], dtype=labile_dtype)),    
}

###
# Create the chemical data container.
###

chemical_data = ChemicalDataContainer(element_data, isotope_data,
                                      amino_acid_data)

#enriched_aa_abundances = np.array([(0.01, 0.99)], dtype=np.float64)
#enriched_aa_fractions = np.array([(0.15)], dtype=np.float64)
