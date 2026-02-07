# Level 0.5 Quick Start 
## tl;dr version for the impatient
1. [Install influxdb 1.8.x](https://docs.influxdata.com/influxdb/v1/) using your OS's package manager, and run it.  You will also need a C compiler (gcc or clang),  CFITSIO, and GNU make installed at a minimum.
   
3. Download the influx database backup (influx-flight-consolidated.tar.gz) from the [calibration dataset](http://soral.as.arizona.edu/GUSTO/calibration/) and unarchive in a temporary folder (say, /path/to/backup).   
   
4. Restore the database using

       influxd restore -portable /path/to/backup

5. Now that the databsse is reconstituted, you no longer need the archive file influx-flight-consolidated.tar.gz or the influx database backup.  Feel free to delete or archive.
   
6. Build the level 0.5 pipeline code in the `src` directory.  Assuming your system does not need any Makefiile adjustments, run `make` (Linux, MacOS) or `gmake` (BSD) as appropriate.   Corrspec should build without errors and, depending on your compiler, should be clean of warnings too.

7. You need the [raw GUSTO correlator lag files](http://soral.as.arizona.edu/GUSTO/data/lags/) as the data input.  Decide where you want to put the input lags and output Level 0.5 data.  In the example below, I assert `/home/obs/data/GUSTO/` as the root folder for data, so the lags go under this directory as `lags` and the output level 0.5 data in `level0.5`.

8.  You can now run the GUSTO pipeline, which automatically invokes the level 0.5 step as needed.  It will generate the proper flags and file lists for use.  However, if you want to run it manually, see below:
  
### Standalone operation ###
To run the level 0.5 step without the pipeline, generate a file list in the level0.5 directory in the code tree and then run the shell script wrapper on it:
   
       ls -1 ../to/lags/*.dat > filelist.txt
       ./run_level05.sh -e -j 4 -f filelist.txt -o /home/obs/data/GUSTO/level0.5
   
where `-j` specifies the number of CPUs to process the data (1 is default), `-f` specifies the file list, `-e` erases the destination folder prior to processing (to keep from appending data to existing files), and `-o` sets the output directory base, creating it if it doesn't exist.
