from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "spinOps.spinOps",  # Module name
        ["spinOps/spinOps.pyx"],  # Path to the .pyx file
        include_dirs=["spinOps/c_code", numpy.get_include()],  # Include directory for C headers
        libraries=[],  # Add any required libraries here
    )
]

setup(
    name="spinOps",
    version="0.1.2",
    author="Philip Grandinetti",
    author_email="your-email@example.com",
    description="A Python package for operations on quantum spin systems.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/pjgrandinetti/SpinOps",
    packages=find_packages(),
    ext_modules=cythonize(extensions),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "Cython",
    ],
)
