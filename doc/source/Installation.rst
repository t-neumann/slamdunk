Installation
============

This section covers the installation of *slamdunk*. Alternatively one could also run *slamdunk* from out of the box :ref:`docker-label` containers.

There are 3 different possibilities:

1. Installation with `conda <https://conda.io/docs/>`_ from the :ref:`conda-label` channel **(recommended)**

2. Installation with `pip <https://pypi.python.org/pypi/pip>`_ from the :ref:`pip-label`

3. Manual installation from :ref:`source-label`

.. _conda-label:

--------
Bioconda
--------

We recommend using virtual environments to manage your Python installation. Our favourite is `Anaconda <https://www.anaconda.com/>`_, a cross-platform tool to manage Python environments. You can installation instructions for Anaconda `here <http://conda.pydata.org/docs/install/quick.html>`_.
Then you can directly install `slamdunk <https://bioconda.github.io/recipes/slamdunk/README.html>`_ from the `Bioconda channel <https://bioconda.github.io/>`_.

""""""""""""
Requirements
""""""""""""

* `Miniconda <https://conda.io/miniconda.html>`_ or `Anaconda <https://www.anaconda.com/download/>`_

""""""""""""
Installation
""""""""""""

Once Anaconda is installed, you can create an environment containing the *slamdunk* package with the following command:

.. code:: bash

    conda create --name <name> -c bioconda slamdunk
    
..

.. _pip-label:

--------------------
Python Package Index
--------------------

`pip <https://pypi.python.org/pypi/pip>`_ is the recommended tool for installing Python packages. It allows a convenient and simple installation
of *slamdunk* from  the `Python Package Index PyPI <https://pypi.python.org/pypi>`_.

""""""""""""
Requirements
""""""""""""

* Python 2.7
* pip
* Java
* cmake
* `R <https://www.r-project.org/>`_ 3.2.2 for *alleyoop* and *splash*

""""""""""""
Installation
""""""""""""

.. code:: bash

    # Having root permissions

    pip install slamdunk
    
    # Latest development version
    
    pip install git+https://github.com/t-neumann/slamdunk.git
    
    # Local user only

    pip install --user slamdunk
    export PATH=$PATH:$HOME/.local/bin/
    
.. _source-label: