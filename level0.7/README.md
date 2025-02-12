# Level 0.7 Quick Start 
## tl;dr version for the impatient
1. Ensure that you have python3 packages for influxdb, configargparse installed
2. If you have already installed and configured influxdb for the level 0.5 processor, skip to step 5.  Otherwise go to step 3.
3. If you haven't done so for level 0.5, [Install influxdb 1.8.x](https://docs.influxdata.com/influxdb/v1/) using your OS's package manager, and run it.
4. Download the influx database backup (influx-flight-consolidated.tar.gz) and the test lags (lags.tar.xz) from the [calibration dataset](http://soral.as.arizona.edu/GUSTO/calibration/) and unarchive in a temporary folder (say, /path/to/backup and /path/to/lags)
5. Restore the database using

       influxd restore -portable /path/to/backup

6. You no longer need influx-flight-consolidated.tar.gz or the influx backup.  Feel free to delete or archive.
7. Edit a L07-config file for the paths suitable for your system.  See the example in ./common/L07-config, copy and edit to suit
8. Ensure the calibration data products IF.txt and dataLog.txt are in your root input data folder pointed to by "datain".  Also untar the udp.tar.gz file in this folder.  

9. Run it on some level 0.5 data.  If you have everything set up in a config file, you can just run it without any arguments:

   time ./L07-pipeline.py

or you can override any of the configuration variables with command line arguments:
   time ./L07-pipeline.py -e -j 4 -b 2 -s 9500 9600

where the arguments are:
      -e will erase the contents of the destination folders before starting
      -j specifies the number of CPUs to use (without it, python pool() will guess)
      -b will process Band 1 or Band 2 data depending on whether you specify 1 or 2
      -s sets the scanID range for the processing
      -i sets the input data path (which should contain the calibration files, level0.5 data folder, and udp folder)
      -o sets the output data root folder (for the level 0.7 data)




   
