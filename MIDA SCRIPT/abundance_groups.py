"""
Abundance group functionality.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012 Casey W. Stark. All Rights Reserved.

  This file is part of `MIDA`.

"""

DEFAULT_CUTOFF = 15

import numpy as np
import scipy.misc

from mida.utils.numerics import binnings

class AbundanceGroup:
    """
    A group of atoms of the same element, used to compute the abundance vs.
    mass distribution.

    """
    def __init__(self, element_id, num_atoms, num_isotopes, isotope_mis):
        self.element_id = element_id
        self.num_atoms = num_atoms
        self.num_isotopes = num_isotopes
        self.isotope_mis = isotope_mis

        # make the combos
        self.combos = binnings(self.num_atoms, self.num_isotopes)
        self.combo_mis = (self.combos * self.isotope_mis).sum(axis=1, dtype=np.int32)

        # the multinomial coefficients (how many ways are there to make each
        # combo)
        coeffs = np.zeros(self.combos.shape)
        # fill in every row of coeffs with the sum of elements in the current
        # bin and previous ones.
        for i in xrange(self.num_isotopes):
            coeffs[:, i] = self.combos[:, 0:(i+1)].sum(axis=1)

        # multinomial coeffs -- row product of the combinations
        self.mn_coeffs = scipy.misc.comb(coeffs, self.combos).prod(axis=1)

    def __repr__(self):
        return "Abundance group: %i atoms of element %i, %i isotopes with masses %s" % (self.num_atoms, self.element_id, self.num_isotopes, self.isotope_mis)

    def __str__(self):
        return self.__repr__()

    def get_combo_abundances(self, abundances):
        """
        Compute the abundances of the combos based on the given isotopic
        abundances.

        Note that the (isotopic) abundances array must be in the shape
        (num_isotopes) or (num_enrichments, num_isotopes). The returned
        combination abundance array is in the shape (num_enrichments,
        num_combos).

        """
        # alias
        combos = self.combos

        # Simple case of one enrichment. The shape still has to be
        # (1, num_combos).
        if len(abundances.shape) == 1:
            return np.array([self.mn_coeffs * (abundances**combos).prod(axis=1)])

        # More than one enrichment.
        # Use numpy.newaxis to broadcast the abundances and combos in separate
        # dimensions. This results in an array with shape (num_enrichments,
        # num_combos, num_isotopes). We then multiply the isotope terms together
        # (axis 2), and multiply the combos by their multinomial coefficients.
        # The final product is an array of the fractional abundance of each
        # combo at every abundances value. The shape is (num_enrichments,
        # num_combos).
        return self.mn_coeffs * (abundances[:, np.newaxis]**combos[np.newaxis, :]).prod(axis=2)

    def get_distribution(self, abundances, mass_cutoff=DEFAULT_CUTOFF):
        """
        Compute the isotopomer distribution of the group.

        """
        # get the abundances of all combos at these isotopic abundances
        combo_abs = self.get_combo_abundances(abundances)

        # only go up to the highest mass combo or the cutoff
        num_mass_bins = min(np.max(self.combo_mis), mass_cutoff) + 1

        # make the distribution array and loop over the bins
        distribution = np.zeros((combo_abs.shape[0], num_mass_bins))
        for mass in xrange(num_mass_bins):
            abundances_at_mass = combo_abs[:, self.combo_mis == mass]
            distribution[:, mass] = abundances_at_mass.sum(axis=1)

        return distribution


class EnrichedAAGroup(AbundanceGroup):
    """
    An abundance group, but specifically for the case of an enriched amino acid.

    """
    def __init__(self, element_id, num_atoms, num_isotopes, isotope_mis):
        AbundanceGroup.__init__(self, element_id, num_atoms, num_isotopes,
                                isotope_mis)

    def get_en_aa_distribution(self, natural_abundances, enriched_abundances,
                               enriched_fraction, mass_cutoff=DEFAULT_CUTOFF):
        """
        Computes the group distribution using the natural and enriched isotopic
        abundances. Then combines the distributions into the total, weighting
        by the enriched fraction.

        """
        natural_dist = self.get_distribution(natural_abundances,
                                             mass_cutoff=mass_cutoff)
        enriched_dist = self.get_distribution(enriched_abundances,
                                              mass_cutoff=mass_cutoff)

        return ( (1.0 - enriched_fraction) * natural_dist
                 + enriched_fraction * enriched_dist )
