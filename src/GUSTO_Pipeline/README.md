# Pipeline Quick Start 
## tl;dr version for the impatient
1. Decide which pipeline steps you intend to run.  If you want to start with level 0.5, you will need a C compiler, CFITSIO, GNU make, and influxdb installed with the GUSTO flight telemetry database.  If you start with level 0.7, you still need influxdb.  If you start at level 0.8 or later, you only need python.
2. Install the pertinent GUSTO data into your local data folder.  The pipeline assumes a data tree structure that looks like this:

| <data root directory>
| |--> lags
| |--> level0.5
| |--> level0.7
| |--> level0.9
| |--> level1
| |--> level2
| |--> logs

The pipeline will create any missing folders and will optionally erase the contents of a folder before performing operations in it.

3.  Change to the root directory of the pipeline (../gusto-datasystem) and run:

.. code-block:: python
``$ pip3 install -e .``

The "-e" switch makes the installation `editable` which means that edits to the code in the local repository are used when executing the pipeline.

4.  Inspect the installed `config.gusto` file and either leave it in place or make a copy for your local edits.  Generally the defaults are OK but you will almost certainly need to alter the `path` setting to your data root directory.  The pipeline will look for the config.gusto file in the source tree or a .config.gusto file in the user's home directory.

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

	runGUSTO -c /tmp/config.gusto -b 1 -s 14000 15000 -l 0.7 0.8

will process Band 1 data from scanIDs 14000 to 15000, performing the Level 0.7 and level 0.8 steps before stopping.  It will read from an externally specified config file.

Maximal (command line overrides most things):

     runGUSTO -e -j 4 -b 1 2 -p /home/obs/data/GUSTO -s 14920 14930 -l 0.5 1.0

will erase output folders before writing, uses 4 CPUs for operation, processes data for both Bands 1 and 2, uses a root data folder of /home/obs/data/GUSTO, will process scanIDs from 14920 to 14930, and will perform all pipeline steps from level 0.5 to level 1.0.  It assumes the standard locations for the config.gusto file. 


### Runtime example

Below is the output from the `maximal` command shown in the last example:

      2025-10-08 22:49:29,589 - pipelineLogger - WARNING - Started logging to /home/obs/data/GUSTO/log/pipeline_20251008224929.log
      08-10-2025 22:49:29 Executing GUSTO Pipeline.
      Namespace(config=None, erase=True, verbose=False, debug=False, path='/home/obs/data/GUSTO/', band=['1', '2'], level=['0.5', '1.0'], cpus='4', scanid=['14920', '14930'], polyorder='3', calmethod='cal_scaledGainHOTs', despurmethod='polyRes', spurchannelfilter=False)
      Command Line Args:   -e -b 1 2 -s 14920 14930 -j 4 -l 0.5 1.0
      Config File (/home/obs/src/gusto-datasystem/src/GUSTO_Pipeline/config.gusto):
      path:              /home/obs/data/GUSTO/
      polyorder:         3
      calmethod:         cal_scaledGainHOTs
      despurmethod:      polyRes
      spurchannelfilter: False


      Executing pipeline levels:  ['0.5', '0.7', '0.8', '0.9', '1.0']
      Scan range for data processing:  [14920, 14930]
      #############################################################################
      Executing Level 0.5: generating power spectra from lags and saving as SDFITS
      Erasing contents of /home/obs/data/GUSTO/level0.5
      Number of cores used for processing: 4

      Processing Band  1
      Processing  310  files, please wait...
      Processing Band  2
      Processing  351  files, please wait...
Level 0.5 to 0.7 done.  661  lag files were processed.

      #############################################################################
      Executing Level 0.7: generating telemetry headers and making sequence files
      Sequences file seemingly exists, skipping step...
      Erasing contents of /home/obs/data/GUSTO/level0.7/
      Processing Band 1 sequence 05702 from 14920 - 14922
      Processing Band 2 sequence 05702 from 14920 - 14922
      Processing Band 1 sequence 05704 from 14924 - 14926
      Processing Band 2 sequence 05704 from 14924 - 14926
      Processing Band 1 sequence 05707 from 14930 - 14932
      ERROR: no data!
      Sequence 05707 NOT OK, flag is 68
      Processing Band 2 sequence 05707 from 14930 - 14932
      ERROR: no data!
      Sequence 05707 NOT OK, flag is 68
      Processing Band 1 sequence 05703 from 14922 - 14924
      Processing Band 2 sequence 05703 from 14922 - 14924
      Processing Band 1 sequence 05705 from 14926 - 14928
      Processing Band 2 sequence 05705 from 14926 - 14928
      Processing Band 1 sequence 05706 from 14928 - 14930
      Processing Band 2 sequence 05706 from 14928 - 14930
      Number of cores used for processing: 2

      Level 0.7 to 0.8 done.  6  sequences were processed.

      #############################################################################
      Executing Level 0.8: channel flags for spurs and row flags for bad fringing
      Erasing contents of /home/obs/data/GUSTO/level0.8/
      Level 0.8 to 0.9 done.  10  SDFITS files were processed.

      #############################################################################

      Wed Oct  8 22:50:02 2025: Executing 0.9: calibrating spectra
      Number of cores used for processing: 4
      Erasing contents of /home/obs/data/GUSTO/level0.9/
      Processing band 1
      Processed: /home/obs/data/GUSTO/level0.8/NII_05702_14920_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/NII_05703_14922_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/NII_05704_14924_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/NII_05705_14926_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/NII_05706_14928_L08.fits
      getCalSpectra: Only REFHOT after OTF available
      getCalSpectra: Only REFHOT after OTF available
      getCalSpectra: Only REFHOT after OTF available
      getCalSpectra: Only REFHOT before OTF available
      getCalSpectra: Only REFHOT before OTF available
      getCalSpectra: Only REFHOT before OTF available

      Processing band 2
      Processed: /home/obs/data/GUSTO/level0.8/CII_05702_14920_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/CII_05703_14922_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/CII_05704_14924_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/CII_05705_14926_L08.fits
      Processed: /home/obs/data/GUSTO/level0.8/CII_05706_14928_L08.fits
      getCalSpectra: Only REFHOT after OTF available
      getCalSpectra: Only REFHOT after OTF available
      getCalSpectra: Only REFHOT before OTF available
      getCalSpectra: Only REFHOT before OTF available
      getCalSpectra: Only REFHOT before OTF available
      getCalSpectra: Only REFHOT before OTF available

      Level 0.9 pipeline done.

      Execution time: 0.04h  2.20m   131.85s
      Processed 10 files.

      #############################################################################
      Executing Level 1.0 pipeline: coordinate corrections
      Number of cores used for processing: 4

      Erasing contents of /home/obs/data/GUSTO/level1/

      Processing Band 1
      Processed: /home/obs/data/GUSTO/level0.9/NII_05702_14920_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/NII_05703_14922_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/NII_05704_14924_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/NII_05705_14926_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/NII_05706_14928_L09.fits

      Processing Band 2
      Processed: /home/obs/data/GUSTO/level0.9/CII_05702_14920_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/CII_05703_14922_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/CII_05704_14924_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/CII_05705_14926_L09.fits
      Processed: /home/obs/data/GUSTO/level0.9/CII_05706_14928_L09.fits
      Level 0.9 to 1.0 done.  10 SDFITS files were processed.

      #############################################################################

      08-10-2025 22:52:16 Done running GUSTO Pipeline.



## Special notes for current version 
- If you are using the influx database from DR1 Rev A or earlier, please download and install the new version (reduced size = 260 MB) with interpolated calibration load THOT values when they are missing: [see link to calibration folder](http://soral.as.arizona.edu/GUSTO/calibration/) and look at the instructions for installing it from the Level 0.5 pipeline readme.
- Finally, either delete any old sequences.txt file for the L07 pipeline to regenerate, or download a new sequences.txt from the calibration folder, above
