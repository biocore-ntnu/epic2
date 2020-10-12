import os
import sys
import pysam

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

__version__ = open("epic2/version.py").readline().split(" = ")[1].replace(
    '"', '').strip()
macros = []

install_requires = [
    "scipy", "numpy", "natsort", "cython", "pysam", "pandas", "pyranges"
]

if sys.version_info[0] == 2:
    install_requires.append("functools32")

if os.getenv("TRAVIS"):
    install_requires.append("coveralls")

if sys.version_info[0] == 2:
    install_requires.append("functools32")

from sys import platform

compile_options = [
    "-Ofast", "-Wall", "-flto"
]  #, "-frename-registers", "-funroll-loops"] # , "-lgzstream", "-lz"

link_options = [
    "-flto"
]

if platform == "linux" or platform == "linux2":
    compile_options.append("-std=c++11")
elif platform == "darwin":
    compile_options.extend("-stdlib=libc++ -std=c++11".split())

from subprocess import check_output

try:
    conda_path = check_output("which conda", shell=True).decode().strip()
except:
    conda_path = ""

conda_include = []
conda_lib = []
if conda_path:
    conda_base = conda_path.replace("bin/conda", "")
    conda_include.append(os.path.join(conda_base, "include"))
    conda_lib.append(os.path.join(conda_base, "lib"))

dir_path = os.path.dirname(os.path.realpath(__file__))
include_dirs = [dir_path + "/epic2/src", dir_path]
# print(conda_lib)

extensions = [
    Extension(
        "epic2.src.cpp_read_files",
        ["epic2/src/cpp_read_files.pyx", "epic2/src/gzstream.cpp"],
        language="c++",
        include_dirs=conda_include + include_dirs,
        library_dirs=conda_lib,
        extra_compile_args=compile_options,
        extra_link_args=link_options,
        libraries=["z"]),
    Extension(
        "epic2.src.reads_to_bins",
        ["epic2/src/reads_to_bins.pyx"],
        language="c++",
        include_dirs=conda_include,
        library_dirs=conda_lib,
        extra_compile_args=compile_options,
        extra_link_args=link_options),
    Extension(
        "epic2.src.SICER_stats", ["epic2/src/SICER_stats.pyx"],
        language="c++",
        extra_compile_args=compile_options,
        extra_link_args=link_options),
    Extension(
        "epic2.src.SICER_stats2", ["epic2/src/SICER_stats2.pyx"],
        language="c++",
        extra_compile_args=compile_options,
        extra_link_args=link_options),
    Extension(
        "epic2.src.statistics", ["epic2/src/statistics.pyx"],
        language="c++",
        extra_compile_args=compile_options,
        extra_link_args=link_options),
    # Extension("epic2.src.differential",
    #           ["epic2/src/differential.pyx"], language="c++",
    #           extra_compile_args=compile_options),
    Extension(
        "epic2.src.find_islands", ["epic2/src/find_islands.pyx"],
        language="c++",
        extra_compile_args=compile_options,
        extra_link_args=link_options),
    Extension(
        "epic2.src.read_bam",
        ["epic2/src/read_bam.pyx"],
        language="c++",
        include_dirs=conda_include + pysam.get_include(),
        library_dirs=conda_lib,
        extra_compile_args=compile_options,
        extra_link_args=link_options + pysam.get_libraries(),
        define_macros=pysam.get_defines(),
        libraries=["z"]),
    Extension(
        "epic2.src.genome_info", ["epic2/src/genome_info.pyx"],
        language="c++",
        extra_compile_args=compile_options,
        extra_link_args=link_options)
]

setup(
    name="epic2",
    packages=find_packages(),
    ext_modules=cythonize(extensions, annotate=True, language_level='2'),
    scripts=["bin/epic2", "bin/epic2-df", "bin/epic2-bw"],
    package_data={
        'epic2': [
            'effective_sizes/*.txt', 'chromsizes/*chromsizes',
            'examples/*.bed.gz'
        ],
        '': ['*.pyx', '*.pxd', '*.h', '*.c', '*.hpp', "*.bed.gz", 'epic2/src/gzstream.cpp']
    },
    version=__version__,
    description="Ultraperformant ChIP-Seq broad peak/domain finder.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/epic2",
    keywords=["ChIP-Seq"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta", "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    include_dirs=["."],
    long_description=
    ("Chip-Seq broad peak/domain finder based on SICER. See the url for more info."
     ))
