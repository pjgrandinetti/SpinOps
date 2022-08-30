from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="myPackage",
    sources=["pySpinOp.pyx"],
    libraries=["spinOp"],
)
setup(
    name="pySpinOp",
    ext_modules=cythonize([examples_extension])
)
