"""
Container for the chemical data.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012, Casey W. Stark. All Rights Reserved.

  This file is part of `MIDA`.

"""

import numpy as np

class ChemicalDataContainer:
    """
    A very simple container class for the sole purpose of importing and passing
    around one object instead of ~5 different arrays.

    """
    def __init__(self, element_data, isotope_data, amino_acid_data):
        self.element_data = element_data
        self.isotope_data = isotope_data
        self.amino_acid_data = amino_acid_data

        self.num_elements = self.element_data.shape[0]
        self.num_isotopes_array = element_data["num_isotopes"]
        self.element_symbols = self.element_data["symbol"]

        ###
        # Create isotope base masses.
        ###

        self.isotope_base_masses = np.zeros(self.num_elements)
        self.isotope_base_nominal_masses = np.zeros(self.num_elements)
        for i in xrange(self.num_elements):
            # get the isotopes of this element
            el_isotopes = self.isotope_data[ self.isotope_data["element_id"] == i ]
            # take the min value of these isotopes' atomic weights
            self.isotope_base_masses[i] = np.min(el_isotopes["mass"])
            self.isotope_base_nominal_masses[i] = np.min(el_isotopes["A"])

        ###
        # Create the isotope M_i tuple.
        ###

        mi_list = []
        for i in xrange(self.num_elements):
            # get the isotopes of this element
            el_isotopes = self.isotope_data[ self.isotope_data["element_id"] == i ]
            # add the mi's to the list
            mi_list.append(el_isotopes["mi"])

        self.isotope_mis = tuple(mi_list)

        ###
        # Create the natural abundances array.
        ###

        nabs_list = []
        for i in xrange(self.num_elements):
            # get the isotopes of this element
            el_isotopes = self.isotope_data[ self.isotope_data["element_id"] == i ]
            nabs_list.append(el_isotopes["natural_abundance"])

        self.natural_abundances = tuple(nabs_list)
