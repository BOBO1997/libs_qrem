import numpy as np
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.core import setup

ext_modules = [
    Extension(
        "mitigation",
        sources=[
            "./libs_qrem/mitigation.pyx",
            "./cpp/eigen_utils.cpp",
            "./cpp/combinations.cpp",
            "./cpp/hamming.cpp",
            "./cpp/sgs_algorithm.cpp",
            "./cpp/qrem_filter_base.cpp",
            "./cpp/qrem_filter.cpp",
            "./cpp/qrem_filter_nlp.cpp",
            "./cpp/qrem_filter_mooney_etal.cpp",
        ],
        extra_compile_args=["-std=c++14"],
        language="c++"
    ),
]

setup(
    name="libs_qrem",
    version="0.1.2",
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
