import os
import sys

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize
macros = []

install_requires = ["scipy", "numpy", "natsort", "cython", "pysam"]


if os.getenv("TRAVIS"):
    install_requires.append("coveralls")


if sys.version_info[0] == 2:
    install_requires.append("functools32")

compile_options = ["-Ofast", "-Wall", "-std=c++11"] #, "-frename-registers", "-funroll-loops"] #

# bamtools_args = ["-I ./lib/bamtools-master/include/ -L ./lib/bamtools-master/src/"]
                   # -fprofile-generate
                   #"-fopenmp", "-D_GLIBCXX_PARALLEL"]


extensions = [Extension("SICER2.src.reads_to_bins",
                        ["SICER2/src/reads_to_bins.pyx"], language="c++", extra_compile_args=compile_options),
              Extension("SICER2.src.statistics",
                        ["SICER2/src/statistics.pyx"], language="c++", extra_compile_args=compile_options),
              Extension("SICER2.src.find_islands",
                        ["SICER2/src/find_islands.pyx"], language="c++", extra_compile_args=compile_options),
              Extension("SICER2.src.read_bam",
                        ["SICER2/src/read_bam.pyx"], language="c++", extra_compile_args=compile_options),
              Extension("SICER2.src.genome_info",
                        ["SICER2/src/genome_info.pyx"], language="c++", extra_compile_args=compile_options)]

setup(
    name="SICER2",
    packages=find_packages(),
    ext_modules = cythonize(extensions, annotate=True),
    scripts=["bin/SICER2"],
    package_data={'SICER2': ['effective_sizes/*.txt', 'chromsizes/*chromsizes'],
                  '': ['*.pyx', '*.pxd', '*.h', '*.c']},
    version="0.0.1",
    description="Ultraperformant ChIP-Seq broad peak/domain finder.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/SICER2",
    keywords=["ChIP-Seq"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment", "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    include_dirs=["."],
    long_description=
    ("Chip-Seq broad peak/domain finder based on SICER. See the url for more info."
     ))
