# GUSTO DATA CUBES 
![CII of Galactic Plane](images/GalacticPlane_CII.png)

## Creating GUSTO L2 data cubes

The module GUSTOgridder.py will create level2  datacubes from level 1 GUSTO data spectra for both CII and NII.  The level 1 data are pipeline baseline corrected and mixer offset pointing corrected.

The regridder is based on the STO2 gridder which was originally based on the Green Bank Telescope of the NRAO gridder for single dish OTF data.
grid\_otf.py and optimized to work with multiple pixel mixer arrays.

The two additional extensions are added to the data cube:  WEIGHT and SCANID.  Weight is a cube of the weights used to grid the L1 scans onto
the final data cube.  The SCANID extenstion is an image made of the SCAN IDs of the L1 data.   


```
    options:
      -h, --help            show this help message and exit
      -v                    show program's version number and exit
      -s --source           Name of source directory in level1. galactic coords maps for source G???, RADEC otherwise
      -b --band             NII or CII
      -k --kernel           gridding kernel: gaussbessel, gauss, or nearest. Default is gaussbessel
      -x --mixer            mixer: NII 2, 3, 6 or CII 5, 8 : or 0 for all mixers in band
      -o --ofile            Output cube name. If none provided a name based on target and regridding parameters is created
      -P --pixBeam          pixels per beam, default is 3
      -Beam --BeamFWHM      beam FWHM in decimal arcmin, default use data header
      -dv --vel_spacing     Velocity spacing in output cube in km/s
      -l --VLSRrange --VLSRrange
                            minimum maximum velocity channel, default -200 200 km/s
      -f --wcsfile          Input fits cube to match WCS if not present the WCS will be made based on input L1 scans

```   
    


