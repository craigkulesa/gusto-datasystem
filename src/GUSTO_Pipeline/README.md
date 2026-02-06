# Pipeline Quick Start 
## tl;dr version for the impatient
1. Decide which pipeline steps you intend to run.  If you want to start with level 0.5, you will need a C compiler, CFITSIO, GNU make, and influxdb installed with the GUSTO flight telemetry database.  If you start with level 0.7, you still need influxdb.  If you start at level 0.8 or later, you only need python.
2. Install the pertinent GUSTO data into your local data folder.  The pipeline assumes a data tree structure that looks like this:

```   
root data directory
  ├── lags
  ├── level0.5
  ├── level0.7
  ├── level0.8
  ├── level0.9
  ├── level1
  ├── log
  └── udp
```

The pipeline will create any missing folders and will optionally erase the contents of a folder before performing operations in it.

3.  Change to the root directory of the pipeline (../gusto-datasystem) and run:

	`$ pip3 install -e .`

The "-e" switch makes the installation `editable` which means that edits to the code in the local repository are used when executing the pipeline.

4. There are a few Python dependencies: pyastronomy, astropy, influxdb (for level 0.7 only), configargparse, and tqdm.  You can install them using your system's package manager (apt, yum, pkgsrc, brew, etc.) , or set up a virtual environment to bypass the system python packages:

```
	$ python3 -m venv venv
	$ . ./venv/bin/activate
	$ pip3 install pyastronomy tqdm astropy influxdb configargparse
```

5.  Inspect the installed `config.gusto` file and either leave it in place or make a copy for your local edits.  Generally the defaults are OK but you will almost certainly need to alter the `path` setting to your data root directory.  The pipeline will look for the config.gusto file in the source tree or a .config.gusto file in the user's home directory.

## Usage 

To run the pipeline, you execute `runGUSTO`.  The combination of the default configuration file, and built-in defaults should allow it to run.  The default will run level 0.9 from level 0.8 data, using 1 CPU and for 10 scanIDs from 14920 to 14930 in both bands 1 and 2.

There are a number of configuration options.  All can be specified in the configuration file and on the command line.  The command line takes highest precedence.

Current arguments are, in no special order:

       -c (config) points to an alterative configuration file
       -e (erase) will erase the contents of any reused destination folders before starting
       -j specifies the number of CPUs to use (without it, python will best-guess)
       -b (band) will either process B1 or B2 data, depending on whether you specify 1 or 2, or both (-b 1 2)
       -s (scanID) sets the scanID range for the processing
       -p (path) sets the root data path (see above)
       -v (verbose) allows for more verbose text output (default is False)
       -d (debug) sets a debug option, which also allows for more output
       -l (level) sets the starting and ending processor level.  If only one parameter is used, the starting and ending levels are assumed to be the same.
       --polyorder sets the baseline fitting polynomial order (1 is default)
       --despurmethod sets the despurring method in level 0.9 ('polyRes' is default)
       --spurchannelfilter automatically applies a channel mask to the spectra (default is False)
       --calmethod selects the calibration method at level 0.9 (cal_scaledGainHOTs is default)

### Examples of use 

Minimal (most parameters come from config file or defaults):

	runGUSTO -c /tmp/config.gusto -b 1 -s 14000 15000 -l 0.7

will process Band 1 data from scanIDs 14000 to 15000, performing the Level 0.7 step before stopping.  It will read from an externally specified config file.

Maximal (command line overrides most things):

     runGUSTO -e -j 4 -b 1 2 -p /home/obs/data/GUSTO -s 14920 14930 -l 0.5 1.0

will erase output folders before writing, uses 4 CPUs for operation, processes data for both Bands 1 and 2, uses a root data folder of /home/obs/data/GUSTO, will process scanIDs from 14920 to 14930, and will perform all pipeline steps from level 0.5 to level 1.0.  It assumes the standard locations for the config.gusto file. 


### Runtime example

Below is the output from the `maximal` command shown in the last example:

```
09-10-2025 19:55:09 Executing GUSTO Pipeline.
Command Line Args:   -e -b 1 2 -s 14920 14930 -j 4 -l 0.5 1.0
Config File (/home/obs/src/gusto-datasystem/src/GUSTO_Pipeline/config.gusto):
  path:              /home/obs/data/GUSTO/
  polyorder:         3
  calmethod:         cal_scaledGainHOTs
  despurmethod:      polyRes
  spurchannelfilter: False
  
2025-10-09 19:55:09 - INFO - Started logging to /home/obs/data/GUSTO/log/pipeline_20251009195509.log
2025-10-09 19:55:09 - INFO - Executing pipeline levels: ['0.5', '0.7', '0.9', '1.0']
2025-10-09 19:55:09 - INFO - Scan range for data processing: [14920, 14930]
#############################################################################

2025-10-09 19:55:09 - INFO - Starting GUSTO Pipeline Level 0.5
2025-10-09 19:55:09 - WARNING - Erasing /home/obs/data/GUSTO/level0.5
2025-10-09 19:55:09 - INFO - Number of cores used for processing: 4
2025-10-09 19:55:09 - INFO - Processing Band 1
2025-10-09 19:55:10 - INFO - Processing 310 files, please wait...
2025-10-09 19:55:20 - INFO - Processing Band 2
2025-10-09 19:55:21 - INFO - Processing 351 files, please wait...
2025-10-09 19:55:36 - INFO - Pipeline step done. 661 sequences or files were processed.
2025-10-09 19:55:36 - INFO - Execution time: 0.01h  0.45m   27.19s

#############################################################################

2025-10-09 19:55:36 - INFO - Starting GUSTO Pipeline Level 0.7
2025-10-09 19:55:36 - WARNING - Erasing /home/obs/data/GUSTO/level0.7/
2025-10-09 19:55:36 - INFO - Number of cores used for processing: 2
2025-10-09 19:55:42 - INFO - Pipeline step done. 6 sequences or files were processed.
2025-10-09 19:55:42 - INFO - Execution time: 0.00h  0.11m   6.41s

#############################################################################

2025-10-09 19:55:42 - INFO - Starting GUSTO Pipeline Level 0.9
2025-10-09 19:55:42 - INFO - Number of cores used for processing: 4
2025-10-09 19:55:42 - WARNING - Erasing /home/obs/data/GUSTO/level0.9/
2025-10-09 19:55:42 - INFO - Processing band 1
2025-10-09 19:56:28 - INFO - Processing band 2
2025-10-09 19:57:55 - INFO - Pipeline step done. 10 sequences or files were processed.
2025-10-09 19:57:55 - INFO - Execution time: 0.04h  2.22m   133.03s

#############################################################################

2025-10-09 19:57:55 - INFO - Starting GUSTO Pipeline Level 1.0
2025-10-09 19:57:55 - INFO - Number of cores used for processing: 4
2025-10-09 19:57:55 - WARNING - Erasing /home/obs/data/GUSTO/level1/
2025-10-09 19:57:55 - INFO - Processing Band 1
2025-10-09 19:57:56 - INFO - Processing Band 2
2025-10-09 19:57:57 - INFO - Pipeline step done. 10 sequences or files were processed.
2025-10-09 19:57:57 - INFO - Execution time: 0.00h  0.03m   1.96s

#############################################################################

09-10-2025 19:57:57 Done running GUSTO Pipeline.
```

## Special notes for current version 
- New flag as of 1 Feb 2026:  mediansubtract.   This snippet of code is from Russ.  This applies a spectral baseline subtraction based on the median of all spectra in the sequence for that mixer.  The current default is False but it is provided for evaluation.  It is certain to improve baselines but we need to be careful to ensure that large scale diffuse emission is not also removed in the background subtraction.
- New quicklook regridder from Chris Martin: L20_regridder.  Codewise, it is merged into the tree but from a user interface standpoint, it remains a standalone module in the pipeline.  It will be integrated more properly in an updated branch, but as an interim workaround, to run it:

  ``` python3 -m GUSTO_Pipeline.L20_regridder.py ```

### Older notes 
- If you are using the influx database from DR1 Rev A or earlier, please download and install the new version (reduced size = 260 MB) with interpolated calibration load THOT values when they are missing: [see link to calibration folder](http://soral.as.arizona.edu/GUSTO/calibration/) and look at the instructions for installing it from the Level 0.5 pipeline readme.
- As of October 2025, the level 0.5-level 1 SDFITS format has changed slightly.  In the data table, there are two new columns for Tsys and rms which are ignored at levels 0.5-0.7 and filled in at level 0.9.  If you are processing data with the new pipeline, either generate new level 0.5 data yourself or use a new starting dataset for level 0.5 or 0.7 (October 2025 or newer).
- Level 0.8 and level 0.95 pipelines have been absorbed into Level 0.7 and Level 0.9 pipeline steps, respectively.  They are no more.  
- The sequences file changd in June 2025.  If you are using an older version, either delete any old sequences.txt file for the L07 pipeline to regenerate, or download a new sequences.txt from the calibration folder, above.
