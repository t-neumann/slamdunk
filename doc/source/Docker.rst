.. _docker-label:

Docker
======

`Docker <https://github.com/t-neumann/slamdunk>`_ enables you to run *slamdunk* out-of-the box wihtout caring about prerequisites and installation procedures.
All you need is to install the Docker engine and run the *slamdunk* Docker image.

Docker images
-------------

You can directly use the *slamdunk* Docker images from `Docker hub <https://hub.docker.com/>`_. Two images are available:

.. * `slamdunk <https://hub.docker.com/r/tobneu/slamdunk/>`_: Build upon NextGenMap's `Docker image <https://hub.docker.com/r/philres/nextgenmap/>`_ (preferred)
* `slamdunk-full <https://hub.docker.com/r/tobneu/slamdunk-full/>`_: Build from scratch from a clean Ubuntu image

You can also build them from scratch using any of the two *Dockerfiles* in `slamdunk/docker`:

..    # Build NGM-based image
..    cd slamdunk/docker/nextgenmap
..    docker build -t <user>/<imagename> .

.. code:: bash
    
    # Build Ubuntu-based image
    cd slamdunk/docker/ubuntu
    docker build -t <user>/<imagename> .
    
    # View images
    docker images


Running the Docker container locally
------------------------------------

Once the container images is pulled/created, *slamdunk* can be run out of the box. The only difference is that the data containing directory has to be mounted on the container and 
all files are relative to that mounting point.
 
.. code:: bash

   # All data is contained in current working directory and available under /data in the container

    docker run -m 8g -v $(pwd):/data tobneu/slamdunk slamdunk map -o /data/map -r /data/GRCm38.fa -5 5 -n 1 /data/testset.fq.gz

Running a Docker container on AWS
---------------------------------

Many cloud providers provide Docker support these days. The most prominent cloud provider are the `Amazon Web Services <https://aws.amazon.com/>`_ (AWS) have a nice integration
into `Docker machine` which easily let's you install Docker engine on AWS machines.
AWS provides free *t2.micro* instances for free-tier accounts. These only come with 1 GiB of memory which is not enough for mammalian genomes (~ 8 GiB needed) where one would need
to go for at least *t2.large* instances.

.. code:: bash

    # Start a free t2.micro instance (small genomes only e.g. bacteria) on e.g. the us-west-1 region
    docker-machine create --driver amazonec2 --amazonec2-region us-west-1 aws-slamdunk
    
or
    
.. code:: bash

    # Start a paid t2.large instance (large genomes) on e.g. the us-west-1 region
    docker-machine create --driver amazonec2 --amazonec2-region us-west-1 --amazonec2-instance-type t2.large aws-slamdunk
    
Once the instance has started, one needs to connect to it before running any Docker images.

.. code:: bash

    eval $(docker-machine env aws-slamdunk)
    
Now everything is set up and one can run slamdunk on an AWS instance.

.. code:: bash

   docker run -m 8g -v $(pwd):/data tobneu/slamdunk slamdunk map -o /data/map -r /data/GRCm38.fa -5 5 -n 1 /data/testset.fq.gz
   
After the work is done, remember to shut down your AWS instances unless you do not care about paying for an idle instance.

.. code:: bash

   # Shut down
   docker-machine stop aws-slamdunk
   
   # Remove
   docker-machine rm -y aws-slamdunk


