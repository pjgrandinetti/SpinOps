from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "spinOps.spinOps",
        ["spinOps/spinOps.pyx"],
        include_dirs=[numpy.get_include(),"c_code"],  # If you're using numpy
        extra_compile_args=["-O3"],  # Optional: extra arguments for the compiler
        extra_link_args=[],  # Optional: extra arguments for the linker
    )
]

setup(
    name="spinOps",
    version="0.1",
    ext_modules=cythonize(extensions),
)
