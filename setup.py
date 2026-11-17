"""
Copyright (c) 2026 MPI-M, Clara Bayley

----- SDM-1D-KiD -----
File: setup.py
Project: SDM-1D-KiD
Created Date: Tuesday 17th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

from setuptools import find_packages, setup

setup(
    name="SDM-1D-KiD",
    version="X.Y.Z",
    packages=find_packages(),
    install_requires=[
        "pytest",
        "sphinx",
    ],
)
