import os
import sys
import subprocess
import setuptools

try:
    import numpy as np
except ImportError:
    subprocess.call([sys.executable, '-m', 'pip', 'install', 'numpy>=1.17'])
    import numpy as np
    
try:
    from Cython.Build import cythonize
except ImportError:
    subprocess.call([sys.executable, '-m', 'pip', 'install', 'cython>=0.29'])
    from Cython.Build import cythonize

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.core import setup

ext_modules = [
    Extension(
        "libs_qrem",
        sources=[
            "./libs_qrem/libs_qrem.pyx",
            "./cpp/eigen_utils.cpp",
            "./cpp/combinations.cpp",
            "./cpp/hamming.cpp",
            "./cpp/sgs_algorithm.cpp",
            "./cpp/qrem_filter_base.cpp",
            "./cpp/qrem_filter_delta.cpp",
            "./cpp/qrem_filter_nlp.cpp",
            "./cpp/qrem_filter_mooney_etal.cpp",
            "./cpp/qrem_filter_lnp.cpp",
        ],
        extra_compile_args=["-std=c++14"],
        language="c++"
    ),
]

setup(
    name="libs_qrem",
    version="0.1.4",
    description="efficient quantum readout error mitigation library",
    cmdclass={"build_ext": build_ext},
    ext_modules=cythonize(ext_modules, language_level=3),
    include_dirs=[np.get_include()],
    zip_safe=False,
    packages=["libs_qrem"],
    package_dir={
        "libs_qrem": "libs_qrem"
    },
)
