from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="pySpinOp",
    sources=["pySpinOp.pyx"],
    libraries=["spinOp"],
)
setup(
    name="pySpinOp",
    ext_modules=cythonize([examples_extension])
)
