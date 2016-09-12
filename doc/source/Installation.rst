Installation
============

This section covers the installation of *slamdunk*. Alternatively one could also run *slamdunk* from out of the box :ref:`docker-label` containers.

1. Install required python modules `pysam`, `joblib`, `pybedtools`, `intervaltree`, `pandas` and `numpy`:

.. code:: bash

    # Having root permissions
    pip install pysam joblib pybedtools intervaltree pandas numpy
    
    # Local user only
    pip install --user pysam joblib pybedtools intervaltree pandas numpy
    
1. Install required R libraries `getopt`, `ggplot2`, `gridExtra`, `RColorBrewer`, `lattice` and `matrixStats`:

.. code:: bash

    R -e 'install.packages(c("getopt","ggplot2","gridExtra","RColorBrewer","lattice","matrixStats"),repos="https://cran.wu.ac.at/")'
    
2. Clone from `Github <https://github.com/t-neumann/slamdunk>`_.

.. code:: bash

    git clone https://github.com/t-neumann/slamdunk.git
    
3. Change to `bin`

.. code:: bash

    cd slamdunk/bin
    
4. Install NGM by following the `build-ngm.sh` instructions.

.. code:: bash

    ./build-ngm.sh

5. Install VarScan2 following the `build-varscan.sh` instructions.

.. code:: bash

    ./build-varscan.sh
    
5. Run *slamdunk* (optionally put it in your *$PATH*  to run it from anywhere).

.. code:: bash

    # Run it from directory
    ./slamdunk --help
   
    # Put it in your $PATH to run it from anywhere
    export PATH=$(pwd):$PATH
   
    slamdunk --help
