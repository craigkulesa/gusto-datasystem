![gusto-hanging-small](https://github.com/user-attachments/assets/ce27f259-b34b-47d1-9655-598948e94dd7)
# GUSTO Ground Data System
GUSTO is a suborbital NASA Explorers Mission of Opportunity launched from McMurdo Station (Antarctica) on a long duration stratospheric balloon on 31 December 2023, flying for a record 57+ days before landing on Mac Robsertson Land near the Antarctic coast between Davis and Mawson stations (AUS).  This repository serves the development of the data system to process the downlinked science data and telemetry into deliverable science products. 
## GUSTO mission and science
GUSTO deployed a 0.9-meter telescope coupled to a 150 liter liquid helium hybrid cryostat which served to keep its suite of three 8-pixel heterodyne mixer arrays cooled to 4-6K over a 50-80 day mission.  Each heterodyne pixel delivers a high resolution spectrum, and the telescope served a mission of drift-scanning the Milky Way and LMC to amass a spectroscopic survey in the far-infrared light of the ionized carbon at 157.7 microns (1900 GHz) and ionized nitrogen at 205 microns (1461 GHz). 

## Data System Documentation
An overview slide sketching the highest-level data flow of the pipeline is shown here:
![pipeline-convergence-9_9_2024](https://github.com/user-attachments/assets/d1f62ea2-58e2-4aa1-a3dc-a75c7d061616)
### Data Level Definitions
#### Raw Level 0 lag files
Level 0 data is pulled from the autocorrelation spectrometers as raw lags coupled with a minimal spectrometer header, with instrument telemetry captured into a standalone Influx database.
#### Level 0.5-0.8 spectral bandpasses
At level 0.5, the lag files are converted to power spectra and incorporated into an SDFITS container by GUSTO scan ID.  At level 0.7, the instrument telemetry is concatenated into the FITS header.  At level 0.8, the data are despurred and channel flags updated.  
#### Level 0.9-1 calibrated spectra
At level 0.9, the bandpasses are normalized by calibration frames and a second round of despurring is performed.  At this point, the job of the calibration HOT and REF is over and only the OTF science frames remain.  At level 1, the data have per-pixel pointing offsets updated, and the frequency axis is converted to velocity in the LSR frame of reference.   The result of this final step generates the final level 1 products. 
#### Level 2 maps
Based on the spectral quality factors assigned before level 1, the level 2 processor assigns weights to each spectrum, and remaps the "random" sampling of the on-the-fly (drift scanning) observations into a regularly sampled grid.   The product of the level 2 processor are 3D FITS cubes with a tile size of the order of 1-2 square degrees. 
## Further reading
The [GUSTO Science Data Management Plan is available here](http://soral.as.arizona.edu/mediawiki/documents/GUSTO/GUSTO-UA-DOC-00043_GUSTO%20Science%20Data%20Management%20Plan_Rev_A.pdf).
