"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: alpha_sampling_wrapper.py
Project: helpers
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Wrapper around PySDM alpha sampling functions
(see generate_inputfiles.py script)
"""

import numpy as np
from PySDM import Formulae
from PySDM.backends import CPU

import sdm_1d_kid.common as cmn


def AlphaSamplingWrapper(probdistrib, alpha, alpha_seed, size_range):
    rng = np.random.default_rng(int(alpha_seed))
    alpha_sampling = cmn.get_alpha_sampling(probdistrib, alpha, size_range=size_range)

    def alpha_sampling_xi(radii, totxi):
        return alpha_sampling_xi.xi

    def alpha_sampling_radii(nsupers):
        backend = CPU(formulae=Formulae(seed=int(rng.integers(0, 2**32))))
        sampling_method = cmn.get_spectral_sampling_method(alpha_sampling)
        radii, alpha_sampling_xi.xi = sampling_method(nsupers, backend=backend)
        return radii

    return alpha_sampling_radii, alpha_sampling_xi
