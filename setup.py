from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "spinOps.spinOps",
        ["spinOps/spinOps.pyx", "spinOps/c_code/spin.c", "spinOps/c_code/spatial.c"], 
        include_dirs=[numpy.get_include(),"spinOps/c_code"],  # If you're using numpy
        extra_compile_args=["-O3"], 
        extra_link_args=['-fPIC'],
        library_dirs=["spinOps/c_code"],

    )
]

setup(
    name="spinOps",
    version="0.1.1",
    url='https://github.com/pjgrandinetti/spinOps',
    description='A Python package for spin operations',
    author='Philip J. Grandinetti',
    author_email='grandinetti.1@osu.edu',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    ext_modules=cythonize(extensions),
)


