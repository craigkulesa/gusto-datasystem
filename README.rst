#GUSTO_Pipeline

GUSTO pipeline
===============================

## Installation

Change to the root directory of level0.9 (.../gusto-datasystem/level0.9) and run:

.. code-block:: python
``$ pip install -e .``

The "-e" switch means that edits to the code in the local repository 
are used when executing the pipeline.   In python lingo, this is an "editable" package installation. 


## Usage

Current use is reading in the Level 0.8 data and creating a Level 0.9 data product:

.. code-block:: python
``runGUSTO``

For limiting the number of files processed, use (change the scan numbers appropriately):

.. code-block:: python
``runGUSTO -s 9000 9600``

More complex command:

.. code-block:: python
``runGUSTO -c /var/tmp/config.gusto -s 19000 19004 -b 2 -j 4 -l 0.7 1.0``

All command line switches can be checked with:

.. code-block:: python
``runGUSTO -h`` 

The program uses a config file including information about where it can find
the data files and where to put them.  The default is to look for the following:

$HOME/.config.gusto
gusto-datasystem/src/GUSTO_Pipeline/config.gusto

or wherever you point to with the -c argument, such as in the example above.

The assumed data tree is like:

| <data root directory>
| ├── <level0.7>
| ├── <level0.8>
| ├── <level0.9>         
| └── <level1>

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that 
this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`GUSTO_Pipeline` was created by the GUSTO data processing team. It is licensed under the terms of the CC0 v1.0 Universal license.

## Credits

To all people who helped ...
