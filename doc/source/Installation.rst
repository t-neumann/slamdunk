Installation
============

This section covers the installation of *slamdunk*. Alternatively one could also run *slamdunk* from out of the box :ref:`docker-label` containers.

There are 2 different possibilities:

1. Installation from `PyPI <https://pypi.python.org/pypi>`_ using the :ref:`pip-label` **(recommended)**

2. Installation from :ref:`source-label`

------------
Requirements
------------

There are no major requirements for *slamdunk*. The python package will acquire all external dependencies by itself.

"""""""""
R runtime
"""""""""

Both *alleyoop* and *splash* utilize `R <https://www.r-project.org/>`_ for various calculations and diagnostic plots.
Therefore, a global installation of R is required to run those tools.

* Minimum required R version: **R 3.0.1**


* Required packages (from `CRAN <https://cran.r-project.org/>`_):

+----------------+
| R packages     |
+================+
| *getopt*       |
+----------------+
| *ggplot2*      |
+----------------+
| *gridExtra*    |
+----------------+
| *RColorBrewer* |
+----------------+
| *lattice*      |
+----------------+
| *matrixStats*  |
+----------------+
| *dplyr*        |
+----------------+
| *tidyr*        |
+----------------+

**Note:** All required R packages will be installed upon the first analysis run of *alleyoop* or *splash*.

.. _pip-label:

--------------------
Python Package Index
--------------------

`pip <https://pypi.python.org/pypi/pip>`_ is the recommended tool for installing Python packages. It allows a convenient and simple installation
of *slamdunk* from  the `Python Package Index PyPI <https://pypi.python.org/pypi>`_.


.. code:: bash

    # Having root permissions
    pip install --extra-index-url https://testpypi.python.org/pypi slamdunk
    
    # Local user only
    pip install --extra-index-url https://testpypi.python.org/pypi --user slamdunk
    export PATH=$PATH:$HOME/.local/bin/
    
    
**Note:** There is no official *slamdunk* release yet, so there will be a test version on `TestPyPI <https://testpypi.python.org/pypi>`_ not guaranteed to be in line with latest developments.
If you want to be absolutely sure you're pulling the latest version, install *slamdunk* from :ref:`source-label`.

.. _source-label:

------
Source
------

1. Clone from `Github <https://github.com/t-neumann/slamdunk>`_.

.. code:: bash

    git clone https://github.com/t-neumann/slamdunk.git

    cd slamdunk

2. Install required python modules:

.. code:: bash

    # Having root permissions
    pip install -r requirements.txt
    
    # Local user only
    pip install --user -r requirements.txt
    
3. Change to `contrib`

.. code:: bash

    cd slamdunk/contrib
    
4. Install NGM by following the `build-ngm.sh` instructions.

.. code:: bash

    ./build-ngm.sh

5. Install VarScan2 following the `build-varscan.sh` instructions.

.. code:: bash

    ./build-varscan.sh

6. Install Samtools following the `build-samtools.sh` instructions.

.. code:: bash

    ./build-samtools.sh

7. Install RNASeqReadSimulator following the `build-rnaseqreadsimulator.sh` instructions.

.. code:: bash

    ./build-rnaseqreadsimulator.sh
    
8. Run *slamdunk* (optionally put it in your *$PATH*  to run it from anywhere).

.. code:: bash

    cd slamdunk/bin 

    # Run it from directory
    ./slamdunk --help
   
    # Put it in your $PATH to run it from anywhere
    export PATH=$(pwd):$PATH
   
    slamdunk --help
