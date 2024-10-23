#gustoL09P

GUSTO L09P pipeline
===============================

## Installation

.. code-block:: python
``$ pip install -e gustoL09P --use-feature=in-tree-build``

The "-e" switch means that edits to the code in the local repository 
are used when executing the pipeline 
Also, create the following environmental variables so the programs can
find the relevant directories, ancillary data, and write the
processing results to the appropriate directories (the synthax below
is for Mac's .zshrc file):

# GUSTO exports for Mac
export GUSTO_ROOT='/Users/<user>/Projects/GUSTO'
export GUSTO_DATA='/Users/<user>/Projects/GUSTO/Data'



## Usage

Current use is reading in the Level 0.8 data and creating a Level 0.9 data product:

- execGL09P

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that 
this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`vtcloud_template_package` was created by vtcloud. It is licensed under the terms of the CC0 v1.0 Universal license.

## Credits

To all people who helped ...