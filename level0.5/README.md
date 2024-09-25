# Level 0.5 Quick Start 
## tl;dr version for the impatient
1. Install influxdb using your OS's package manager, and run it. 
2. Download the influx database backup (influx-flight-consolidated.tar.gz) and the test lags (lags.tar.xz) from the [calibration dataset](http://soral.as.arizona.edu/GUSTO/calibration/) and unarchive in a temporary folder (say, /path/to/backup and /path/to/lags)
3. Restore the database using

       influxd restore -portable /path/to/backup

4. You no longer need influx-flight-consolidated.tar.gz or the influx backup.  Feel free to delete or archive.
5. Generate a file list and then run the processor on it:
   
       ls -1 ../to/lags/*.dat > filelist.txt
       ./run_level05.sh -c -j 4 -f filelist.txt
   
where `-j` specifies the number of CPUs to process the data (1 is default), `-f` specifies the file list, and `-c` cleans the destination folder prior to processing (to keep from appending data to existing files).


   
