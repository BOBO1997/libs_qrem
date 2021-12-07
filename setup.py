from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.core import setup

ext_modules = [
    Extension(
        "hamming",
        sources=[
            "./libs_qrem/hamming.pyx",
            "./cpp/combinations.cpp",
            "./cpp/hamming.cpp"
        ],
        extra_compile_args=["-std=c++11"],
        language="c++"
    ),
    Extension(
        "sgs_algorithm",
        sources=[
            "./libs_qrem/sgs_algorithm.pyx",
            "./cpp/sgs_algorithm.cpp"
        ],
        extra_compile_args=["-std=c++11"],
        language="c++"
    ),
    Extension(
        "mitigation",
        sources=[
            "./libs_qrem/mitigation.pyx",
            "./cpp/mitigation.cpp",
            # "./libs_qrem/qrem_filter.pyx"
        ],
        extra_compile_args=["-std=c++11"],
        language="c++"
    ),
]

setup(
    name="libs_qrem",
    cmdclass={"build_ext": build_ext},
    ext_modules=cythonize(ext_modules, language_level=3),
    zip_safe=False,
    packages=["libs_qrem"],
    package_dir={
        "libs_qrem": "libs_qrem"
    },
)
