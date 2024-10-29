Installing `GUSTOL08`
============================

Requirements
------------

This package has the following dependencies:

* `Python <http://www.python.org>`_ 3.8
* `Numpy <http://www.numpy.org>`_ 1.19 or later
* `Astropy <http://www.astropy.org>`__ 4.0 or later
* `scipy <https://www.scipy.org>`__ 1.7 or later
* `sphinx`__ 3.2 or later
* `myst-parser`__ 0.16 or later
* `readthedocs-sphinx-search`__ 0.1.1 or later
* `argparse`__ 1.4 or later

Installation
------------

Currently, the only path to install the `GUSTOL08` package is
from the git repositor (you must be logged-in)::

    git clone https://github.com/vtcloud/GUSTOL08.git
    cd GUSTOL08
    pip install .


The package needs a configuration file with the paths to several data
files. This configuration file, `GL2Pconfig.txt` is placed in the directory
`${HOME}/.gusto` with ${HOME} the environmental variable to the users home
directory. The file contains the following full path information:


###### Alternate installation or updating method:

You can also install the latest developer version in a single line with pip.
However, there might be a few adjustments to be made in order for the
program to work correctly,- see Installation::

    pip install git+https://github.com/vtcloud/GUSTOL08.git
