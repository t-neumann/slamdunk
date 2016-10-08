Installation
============

This section covers the installation of *slamdunk*. Alternatively one could also run *slamdunk* from out of the box :ref:`docker-label` containers.

------------
Requirements
------------

There are no major requirements for *slamdunk*. The python package will acquire all external dependencies by itself.

Both *alleyoop* and *splash*, however, utilize `R <https://www.r-project.org/>`_ for various calculations and diagnostic plots.

Therefore, a global installation of R is required to run those tools (minimum version **R 3.0.1**).

All required R libraries listed below will be installed upon the first run of any tool that requires them.

=================== 
Required library              
===================  
*getopt*
*ggplot2*
*gridExtra*
*RColorBrewer*
*lattice*
*matrixStats*
*dplyr*
*tidyr*
===================

-----------------
Pip (recommended)
-----------------

There is no official release yet, so there will be test versions on `TestPyPI <https://testpypi.python.org/pypi>`_ not guaranteed to be in line with latest developments.

If you want to be absolutely sure you're pulling the latest version, check out the instructions on how to install *slamdunk* from source.

.. code:: bash

    # Having root permissions
    pip install --extra-index-url https://testpypi.python.org/pypi slamdunk
    
    # Local user only
    pip install --extra-index-url https://testpypi.python.org/pypi --user slamdunk

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
