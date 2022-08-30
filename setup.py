from setuptools import Extension
from setuptools import setup
from Cython.Build import cythonize

examples_extension = Extension(
    name="myPackage",
    sources=["pySpinOp.pyx"],
    libraries=["spinOp"],
    library_dirs=["."],
)
setup(ext_modules=cythonize([examples_extension]))
