from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize

class build_ext(_build_ext):
    def finalize_options(self):
        super().finalize_options()
        import numpy  # ← ここで遅延 import（ビルド環境に既に入っている）
        self.include_dirs.append(numpy.get_include())

ext_modules = [
    Extension(
        "libs_qrem",
        sources=[
            "libs_qrem/libs_qrem.pyx",
            "cpp/eigen_utils.cpp",
            "cpp/combinations.cpp",
            "cpp/hamming.cpp",
            "cpp/sgs_algorithm.cpp",
            "cpp/harger_higham.cpp",
            "cpp/qrem_filter.cpp",
            "cpp/ignis_filter.cpp",
            "cpp/delta_filter.cpp",
            "cpp/least_norm_filter.cpp",
            "cpp/mooney_etal_filter.cpp",
            "cpp/nation_etal_filter.cpp",
        ],
        language="c++",
        extra_compile_args=["-std=c++17", "-pthread"],
        include_dirs=["./eigen"],  # numpy の include は build_ext で追加
    ),
]

setup(
    name="libs_qrem",
    version="0.1.6",
    description="efficient quantum readout error mitigation library",
    cmdclass={"build_ext": build_ext},
    ext_modules=cythonize(ext_modules, language_level=3),
    packages=["libs_qrem"],
    package_dir={"libs_qrem": "libs_qrem"},
    zip_safe=False,
)
