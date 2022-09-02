from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

examples_extension = Extension(
    name="spinOps",
    sources=["SpinOps.pyx"],
    libraries=["spinOps"],
    library_dirs=["lib"],
    include_dirs=["lib",numpy.get_include()]
)
setup(
    name="spinOps",
    ext_modules=cythonize([examples_extension])
)
