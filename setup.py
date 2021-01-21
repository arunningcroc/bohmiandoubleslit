from setuptools import setup
from Cython.Build import cythonize
import numpy as np
setup(
    name="Fast potentials",
    ext_modules=cythonize("fastpotential.pyx", compiler_directives={'language_level': "3", 'boundscheck': False,
     'cdivision': True},annotate=True),
    include_dirs=[np.get_include()],
)
