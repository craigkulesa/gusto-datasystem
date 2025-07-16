# Level 0.8 Quick Start
## tl;dr version for the impatient
1. There are relatively few dependencies for the Level 0.8 pipeline step.  If you have a working python3 setup with astropy and a few other minor packages, it should run without concern.
2. You should have your data in a common path such as $HOME/data, and the level 0.7 data will be expected to be found in a "level0.7" directory, and the level 0.8 data will go into "level0.8".
3. Feel free to modify the config file (default is config.ini in the common/ folder) as needed.  Presently, it is used only to remember your data path.
4. Run the processor on some Level 0.7 data.  You can accept the defaults from the config file and simply the pipeline step from the level0.8/src directory without any arguments:

       time ./L08-pipeline.py
or you can override the defaults on the command line.   The example below will use 2 CPUs, process scanIDs from 1 to 30000 (e.g. the whole survey) from the level0.7 data, writing new files, and erases the output level 0.8 directory before starting.

       time ./L08-pipeline.py -j 2 -e -s 1 30000

       Namespace(config=None, erase=True, path='/home/obs/data/GUSTO/', update=False, cpus='2', force=False, scanid=['1', '30000'])
       Command Line Args:   -j 2 -e -s 1 30000
       Config File (../../common/config.ini):
       	      path:              /home/obs/data/GUSTO/

      Reusing directory /home/obs/data/GUSTO/level0.8/
      Erasing contents of /home/obs/data/GUSTO/level0.8/
      Using level 0.7 folder for input and creating new level 0.8 files
      Found 23728 files to process
      Done.
      
      real    47m7.321s
      user    89m3.860s
      sys     1m43.046s

where the arguments are (very similar to the Level 0.7 pipeline):

      -c or --config will specify the location of the config file
      -e or --erase will wipe the contents of the level0.8 directory before starting
      -p or --path with set the data path where the level0.7 and level0.8 directories are
      -j or --cpus will specify the number of cpu cores to use
      -f or --force will force channel and row flagging, even if already performed
      -s or --scanid has two arguments and will set the range of scanIDs to process
      -u or --update will update files in level0.8 folder "in place" (behavior of Abe's scripts)

Not specifying -u is the default; the pipeline will look in level 0.7 for the input files, and write new data to the level 0.8 directory.

## Notes for current version ##
- Note that the -e and -u flags conflict -- you don't want to wipe the level 0.8 directory while still expecting to process files in-place there.  At present, the code will not stop you from doing this. 


