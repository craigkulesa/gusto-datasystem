#gustoL09P

GUSTO L09P pipeline
===============================

## Installation

Change to the root directory of level0.9 (.../gusto-datasystem/level0.9) and run:

.. code-block:: python
``$ pip install -e .``

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

For limiting the number of files processed, use (change the scan numbers appropriately):

- execGL09P -r 9000 9600

The program uses a parameter setup file including information about where it can find
the data files and where to put them. The assumed data tree is like:

| <data root directory>
| ├── <level 0.7>
| ├── <level 0.8>
| ├── <level 0.9>         
| ├── <level 0.95>
| └── <level 1.0>
 

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that 
this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`gustoL09P` was created by vtcloud. It is licensed under the terms of the CC0 v1.0 Universal license.

## Credits

To all people who helped ...