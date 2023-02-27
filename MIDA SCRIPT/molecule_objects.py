"""
Molecular data objects. General molecules and Peptides. AminoAcid objects are
just containers for the amino acid data.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012, Casey W. Stark. All Rights Reserved.

  This file is part of `MIDA`.

"""

DEFAULT_CUTOFF = 15

import numpy as np

from mida.abundance_groups import AbundanceGroup, EnrichedAAGroup
from mida.data_types import composition_dtype, labile_dtype, aa_enrichment_dtype


class Molecule:
    """
    A molecule with a specific chemical formula (elements, but not specific
    isotopes).

    """
    def __init__(self, composition, chemical_data, labiles=None,
                 aa_enrichments=None):
        """
        Create the molecule with the supplied elemental `composition`.
        Splits the atoms into abundances groups to compute the distribution
        given isotopic abundances for labile and enriched amino acid groups.

        """
        # attach arguments
        self.composition = composition
        self.labiles = labiles
        self.aa_enrichments = aa_enrichments

        self.chemical_data = chemical_data

        # check if there are any labile or amino acid enrichement groups
        self.active_labile_groups = False
        self.active_en_aa_groups = False
        if self.labiles is not None:
            self.active_labile_groups = True
        if self.aa_enrichments is not None:
            self.active_en_aa_groups = True

        ###
        # Make the abundances groups.
        ###

        # alias for shorthand and faster access in loops
        num_isotopes_array = self.chemical_data.num_isotopes_array
        isotope_mis = self.chemical_data.isotope_mis

        # Lists of the abundance groups
        self.na_groups = []
        self.labile_groups = []
        self.en_aa_groups = []

        # keep track of how many atoms of each element are left at natural
        # abundances
        self.na_composition = self.composition.copy()

        if self.active_labile_groups:
            # go through labile groups and make the abundance group objects,
            # correcting the na_composition along the way
            for group in self.labiles:
                # alias
                element_id = group["element_id"]
                # note that group["n"] can be a float. We need to round and int
                # cast before passing it on.
                num_atoms = int(round(group["n"]))

                self.labile_groups.append(AbundanceGroup(element_id, num_atoms,
                    num_isotopes_array[element_id], isotope_mis[element_id]))

                # subtract off from na_composition
                self.na_composition[element_id] -= num_atoms

        if self.active_en_aa_groups:
            # repeat the process for any amino acid enrichment groups
            for group in self.aa_enrichments:
                # alias
                element_id = group["element_id"]
                # note that group["n"] can be a float. We need to round and int
                # cast before passing it on.
                num_atoms = int(round(group["n"]))

                self.en_aa_groups.append(EnrichedAAGroup(element_id, num_atoms,
                    num_isotopes_array[element_id], isotope_mis[element_id]))

                self.na_composition[element_id] -= num_atoms

        # now make the groups at natural abundances
        for i, num_atoms in enumerate(self.composition):
            if num_atoms != 0:
                self.na_groups.append(AbundanceGroup(i,
                    self.na_composition[i], num_isotopes_array[i],
                    isotope_mis[i]))

    def get_distribution(self, labile_abundances=None,
                         en_aa_abundances=None, en_aa_fraction=None,
                         mass_cutoff=DEFAULT_CUTOFF):
        """
        Get distribution of all the groups and combine them into the total
        distribution for this molecule.

        """
        # List to store the distribution arrays in. We use a list here because
        # the distributions can be different shapes (much messier to handle for
        # modest performance boost).
        distributions = []

        natural_abs = self.chemical_data.natural_abundances

        # Compute natural abundances distributions. Note that the shape is
        # *always* (1, max_mass), but we broadcast to
        # (num_enrichments, max_mass).
        for group in self.na_groups:
            # compute the distribution
            dist = group.get_distribution(natural_abs[group.element_id],
                                          mass_cutoff=mass_cutoff)

            #if dist.shape[0] != num_enrichments:
            #    dist = np.ones((num_enrichments, dist.shape[1])) * dist

            # append the distribution
            distributions.append(dist)

        if self.active_labile_groups:
            # Note that we rely on the order of the groups and abundances being
            # the same!
            for group, group_it_abundances in zip(self.labile_groups,
                                                  labile_abundances):
                distributions.append(
                    group.get_distribution(group_it_abundances,
                                           mass_cutoff=mass_cutoff))

        if self.active_en_aa_groups:
            for group, group_it_abundances, group_en_fraction \
            in zip(self.en_aa_groups, en_aa_abundances, en_aa_fraction):
                # compute and append the distribution to aa_dists
                distributions.append(
                    group.get_en_aa_distribution(natural_abs[group.element_id],
                        group_it_abundances, group_en_fraction,
                        mass_cutoff=mass_cutoff))

        ###
        # Combine distributions of all groups
        ###

        # the `+ 1` is to account for the 0-th mass bin.
        max_mass_bins = mass_cutoff + 1

        # iterate over all distributions, combining them by mass. Note that we
        # catch the first iteration to create the `total_dist` var.
        for i, dist in enumerate(distributions):
            # first time through -- make the total_dist
            if i == 0:
                total_dist = dist

                # go on to the next distribution
                continue

            # find the sizes we will iterate over, paying attention to the max
            # number of mass bins.
            total_dist_size = total_dist.shape[1]
            dist_size = min(max_mass_bins, dist.shape[1])
            # the `- 1` is because the 0-th bin can only count once. When we
            # add the length of two distributions, we must get rid of one.
            new_dist_size = min(max_mass_bins, total_dist_size + dist_size - 1)

            num_enrichments = max(dist.shape[0], total_dist.shape[0])

            # Create a fresh array for the combined distribution.
            new_dist = np.zeros((num_enrichments, new_dist_size))

            for i in xrange(total_dist_size):
                for j in xrange(dist_size):
                    mass = i + j
                    if mass < new_dist_size:
                        new_dist[:, mass] += total_dist[:, i] * dist[:, j]
                    else:
                        break  # `j` will only increase, so go to next `i`

            # save the new abundances
            total_dist = new_dist

        return total_dist

    def __repr__(self):
        return self.formula

    def __str__(self):
        return self.__repr__()

    ### End of dunder methods

    @property
    def base_mass(self):
        """ Returns mass of the base isotopoologue. """
        if not hasattr(self, "_base_mass"):
            self._base_mass = (self.composition * self.chemical_data.isotope_base_masses).sum()
        return self._base_mass

    @property
    def nominal_mass(self):
        """
        Returns nominal mass (rounded and truncated) of the base isotopoologue.

        """
        if not hasattr(self, "_nominal_mass"):
            self._nominal_mass = int(round(self.base_mass))
        return self._nominal_mass

    # @todo: update for the new composition format
    @property
    def formula(self):
        """ Chemical formula for this molecular composition. """
        if not hasattr(self, "_formula"):
            formula_string = ""
            for i, num_atoms in enumerate(self.composition):
                symbol = self.chemical_data.element_symbols[i]
                formula_string += "%s%i" % (symbol, num_atoms)
            self._formula = formula_string
        return self._formula

class Peptide(Molecule):
    """ Representation of a peptide, a sequence of amino acids. """

    def __init__(self, sequence, chemical_data, h_index=0, o_index=3):
        # save the sequence, then process it
        self.sequence = sequence

        # We don't want to assume anything about the structure of the amino
        # acid data arrays. Instead, we use the first instance we find inside
        # the loop below.
        composition = None
        labiles = None
        aa_enrichments = None

        # shorthand for loop
        aa_data = chemical_data.amino_acid_data

        for i, aa_key in enumerate(sequence):
            aa = aa_data[aa_key]

            if composition is None:  # first iteration
                composition = aa.composition.copy()
            else:
                # we can add compositions directly
                composition = composition + aa.composition

            # make sure the aa object has labile groups before we worry about
            # it.
            if aa.labiles is not None:
                if labiles is None:  # first iteration
                    labiles = aa.labiles.copy()
                else:
                    # only add up the "n" column of these
                    labiles["n"] += aa.labiles["n"]

            # make sure the aa object has amino acid enrichment groups before we
            # worry about it.
            if aa.aa_enrichments is not None:
                if aa_enrichments is None:  # first iteration
                    aa_enrichments = aa.aa_enrichments.copy()
                else:
                    aa_enrichments["n"] += aa.aa_enrichments["n"]

            # dehydrogenation
            # @todo: make this smarter?
            if i != 0:
                composition[h_index] -= 2
                composition[o_index] -= 1

        # now init the molecule with the composition and groups generated here.
        Molecule.__init__(self, composition, chemical_data, labiles=labiles,
                          aa_enrichments=aa_enrichments)

    def __repr__(self):
        return "%s Peptide" % self.sequence

class AminoAcid(Molecule):
    """ Placeholder for amino acid data. """

    def __init__(self, one_code, composition, labiles=None,
                 aa_enrichments=None):
        self.one_code = one_code
        self.composition = composition
        self.labiles = labiles
        self.aa_enrichments = aa_enrichments

    def __repr__(self):
        return "%s Amino Acid" % self.one_code

    def __str__(self):
        return self.__repr__()

    def get_num_elements(self):
        return len(self.composition)
