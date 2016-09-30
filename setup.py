"""Slamdunk: SLAM-seq data analysis software 
Template taken from:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

import os, sys
try:
    from setuptools import setup, find_packages
    from setuptools.command.install import install as _install
    from setuptools.command.sdist import sdist as _sdist
    from setuptools.command.bdist_egg import bdist_egg as _bdist_egg
    from codecs import open
    from os import path
except ImportError:
    from distutils.core import setup
    from distutils.command.install import install as _install
    from distutils.command.sdist import sdist as _sdist
    from distutils.command.bdist_egg import bdist_egg as _bdist_egg

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()
    
def _runExternalBuilds(dir):
    print("BUILD:")
    print(dir)
    
class install(_install):
    def run(self):
        #from subprocess import call
        #call(["pip install -r requirements.txt --no-clean"], shell=True)
        self.execute(_runExternalBuilds, (self.install_lib,),msg="Installing stuff")
        _install.run(self)
        #self.execute(_run_build_tables, (self.install_lib,),
        #            msg="Build the lexing/parsing tables")
        
class bdist_egg(_bdist_egg):
    def run(self):
        #print("CALLING BDIST")
        #call(["pip install -r requirements.txt --no-clean"], shell=True)
        _bdist_egg.run(self)


class sdist(_sdist):
    def make_release_tree(self, basedir, files):
        print("Basedir: " + basedir)
        _sdist.make_release_tree(self, basedir, files)
        self.execute(_runExternalBuilds, (basedir,),msg="Dsitributing stuff")
#         self.execute(_run_build_tables, (basedir,),
#                      msg="Build the lexing/parsing tables")

setup(
    name='slamdunk',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.1.0',

    description='SLAMdunk application for analyzing SLAM-seq data',
    long_description=long_description,

    # The project's main homepage.
    url='https://slamdunk.readthedocs.org',

    # Author details
    author='Tobias Neumann, Philipp Rescheneder',
    author_email='tobias.neumann.at@gmail.com, philipp.rescheneder@univie.ac.at',

    # Choose your license
    license='GNU General Public License v3 (GPLv3)',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

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
        'Programming Language :: Java',
    ],

    # What does your project relate to?
    keywords='Next-Generation-Sequencing NGS QuantSeq SLAMSeq',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'doc', 'tests']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    py_modules=["slamdunk.main", "slamdunk.toolbox","slamdunk.simulate"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['joblib','pybedtools','intervaltree','pandas','numpy','biopython'],

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
#     package_data={
#         'sample': ['package_data.dat'],
#     },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
#     data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    scripts= ['bin/slamdunk', 'bin/alleyoop', 'bin/slamsim'],
    
    cmdclass={'install': install, 'sdist': sdist, 'bdist_egg': bdist_egg},
)