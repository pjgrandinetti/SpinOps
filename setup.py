import glob
import os
import platform
import sys
from distutils import ccompiler  # ensure this import before build_ext
from distutils.sysconfig import get_python_inc
from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup
from setuptools.command.install import install as _install
import shutil, glob

# ---- C & Cython sources -------------------------------------------------
c_sources = glob.glob("spinOps/c_code/*.c")
pyx_source = "spinOps/spinOps.pyx"
sources = [pyx_source, *c_sources]

# choose compiler and language flags
# On Windows, force use of MinGW-w64 gcc; on all platforms compile as C99
if sys.platform == "win32":
    from distutils.sysconfig import get_config_var

    ptr_size = get_config_var("SIZEOF_VOID_P") or __import__("struct").calcsize("P")
    ptr_size = int(ptr_size)
    os.environ["CC"] = (
        "i686-w64-mingw32-gcc" if ptr_size == 4 else "x86_64-w64-mingw32-gcc"
    )
    os.environ["CXX"] = (
        "i686-w64-mingw32-g++" if ptr_size == 4 else "x86_64-w64-mingw32-g++"
    )
# Compile all sources as C99
ext_language = "c"
ext_args = ["-std=gnu99"]
# Always use MinGW on Windows, keep dynamic linking for MinGW runtimes
ext_link_args = ["-static", "-static-libgcc", "-static-libstdc++"] if sys.platform == "win32" else []

# configure define macros
define_macros = [
    ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")
]  # only numpy macro, drop SIZEOF_VOID_P override on Windows

# configure Cython extension with platform-specific flags
ext_modules = cythonize(
    Extension(
        "spinOps._spinOps",
        sources=sources,
        include_dirs=[np.get_include(), "spinOps/c_code", get_python_inc()],
        define_macros=define_macros,
        extra_compile_args=ext_args,
        extra_link_args=ext_link_args,
        language=ext_language,
    ),
    language_level="3",
    force=True,
)

from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler

# Custom build_ext: on Windows use MinGW-w64, else default
from setuptools.command.build_ext import build_ext as _build_ext


class CustomBuildExt(_build_ext):
    def build_extensions(self):
        if sys.platform == "win32":
            compiler = new_compiler(compiler="mingw32")
            customize_compiler(compiler)
            self.compiler = compiler
        super().build_extensions()


class CustomInstall(_install):
    def run(self):
        super().run()
        if sys.platform == "win32":
            gcc_bin = os.path.dirname(shutil.which("gcc"))
            dlls = glob.glob(os.path.join(gcc_bin, "libgcc_s_*.dll")) + glob.glob(os.path.join(gcc_bin, "libstdc++-6.dll"))
            pkg_dir = os.path.join(self.install_lib, "spinOps")
            for dll in dlls:
                shutil.copy(dll, pkg_dir)


# ---- Package metadata ---------------------------------------------------
setup(
    name="SpinOps",
    version="0.1.1",
    description="Angular-momentum / spin operators in C + Cython",
    author="Philip Grandinetti",
    license="MIT",
    python_requires=">=3.12",
    install_requires=["numpy>=1.21,<2", "matplotlib>=3.3.4"],
    packages=find_packages(include=["spinOps*"]),
    ext_modules=ext_modules,
    cmdclass={"build_ext": CustomBuildExt, "install": CustomInstall},
    include_package_data=True,
    zip_safe=False,
)
