"""
Provides definitions of numpy data types used in the package.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012, Casey W. Stark. All Rights Reserved.

  This file is part of `MIDA`.

"""

import numpy as np

# element data
# format is the symbol, atomic number, and the number of (stable) isotopes.
ed_dtype = np.dtype([("symbol", "a2"), ("Z", np.int32),
                     ("num_isotopes", np.uint32)])

# isotope data
# format is element index, A, mi (A - A0), mass, natural abundance
id_dtype = np.dtype([("element_id", np.int32), ("A", np.int32), ("mi", np.int32),
                     ("mass", np.float64), ("natural_abundance", np.float64)])

# elemental composition of a molecule
# the number of atoms of each element must be an integer
composition_dtype = np.int32

# labile groups
# format is the group label, which element, and the number of atoms
labile_dtype = np.dtype([("element_id", np.int32), ("n", np.float64)])

# amino acid enrichment groups
# format is the group label, which element, and the number of atoms
aa_enrichment_dtype = np.dtype([("element_id", np.int32), ("n", np.float64)])
