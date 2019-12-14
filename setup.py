# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
#
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import os, sys, re

try:
    from setuptools import setup, find_packages
    from setuptools.command.install import install as _install
    from codecs import open
    from os import path
except ImportError:
    from distutils.core import setup
    from distutils.command.install import install as _install

here = path.abspath(path.dirname(__file__))
name = "slamdunk"

#Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# now we have a `__version__` variable
exec(open(path.join(here, name, 'version.py')).read())

# Copy bin recursively
def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join("..", "..", path, filename))
    return paths

bin_files = package_files(name + os.sep + 'contrib')
plot_files = package_files(name + os.sep + 'plot')

def _runExternalBuilds(dir):

    import subprocess

    print("Building RNASeqReadSimulator.")
    syscall = "(cd " + os.path.join(dir, name, "contrib") + " ; ./build-rnaseqreadsimulator.sh)"
    subprocess.call([syscall], shell=True)

class install(_install):

    def initialize_options(self):
        _install.initialize_options(self)

    def finalize_options(self):
        _install.finalize_options(self)

    def run(self):
        _install.run(self)
        self.execute(_runExternalBuilds, (self.install_lib) ,msg="Installing external dependencies")

setup(
    name = name,

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='SLAMdunk suite for analyzing SLAM-seq data',
    long_description=long_description,

    # The project's main homepage.
    url='http://t-neumann.github.io/slamdunk',

    # Author details
    author='Tobias Neumann, Philipp Rescheneder',
    author_email='tobias.neumann.at@gmail.com, philipp.rescheneder@univie.ac.at',

    # Choose your license
    license='GNU Affero General Public License v3 or later (AGPLv3+)',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: C++',
        'Programming Language :: Java',
    ],

    # What does your project relate to?
    keywords='Next-Generation-Sequencing NGS QuantSeq SLAMSeq',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['doc', 'tests']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #py_modules=["slamdunk.main", "slamdunk.toolbox","slamdunk.simulate"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['joblib>=0.9.4','pybedtools>=0.6.4','intervaltree>=2.1.0','pandas>=0.13.1','biopython>=1.63','pysam>=0.8.3', 'Cython>=0.20.1'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
#     extras_require={
#         'dev': ['check-manifest'],
#         'test': ['coverage'],
#     },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'slamdunk.contrib': bin_files,
        'slamdunk.plot': plot_files,
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('bin', extra_files)],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
    'console_scripts': [
        'slamdunk=slamdunk.slamdunk:run',
        'alleyoop=slamdunk.alleyoop:run',
        'splash=slamdunk.splash:run',
    ],
    },
)
