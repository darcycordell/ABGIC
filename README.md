# ABGIC
Geomagnetically Induced Currents Hazard Analysis Using MT Impedances

These scripts are designed to compute geoelectric fields over Alberta and British Columbia (Canada) using both a 1-D piecewise continuous conducitivty model (1-D impedance method) and an empirical magnetotelluric data array (3-D impedance method), given an incoming geomagnetic disturbance. Geoelectric fields are computed in the frequency-domain and then converted to the time domain. Finally, time-domain fields are used to compute path line integrals along the >240 kV transmission line network in Alberta.

The manuscript has been submitted to Space Weather and a DOI link will be provided here if the manuscript is accepted.

GIC SCRIPTS

Organization:

GIC_MAIN: This script is the main script to run. The first ~50 lines of code 
    contain user inputs with instructions on what options to use to re-create
    the figures in the paper. Once the user inputs are set, then the script
    can be run from the Command Window.

GIC_PLOTTING: After running GIC_MAIN, you will have all the necessary variables
    the MATLAB workspace to create the figures. Each figure is contained in a
    single code block although the whole script can also be run in the Command Window.

Folders-------------

000_DATA: This folder contains all the relevant data for the study;
    01_GEOGRAPHICAL_DATA: shape files for Canadian province outlines
    02_TRICHTCHENKO_ZONES: geojson files with outlines of piecewise continuous zones
    03_POWERLINE_DATA: Google Earth kml file with transmission line paths
    04_MAG_DATA: Magnetic field observatory data from NRCan, USGS, and CARISMA
        networks in IAGA format.
    05_MT_IMPEDANCE: Magnetotelluric impedance data saved as a MATLAB data file
        containing a structure "d" with all the variables related to the dataset
    06_FORWARD_MODELS: Folder containing the models, forward modelled data, and 
        other relevant files to run ModEM forward simulation.
        Two folders: one which contains the SABC and one which does not.


calc: Folder which contains all main calculations for E, MT data, line integrals, etc.

interpolation: Interpolation functions which interpolate in time (b) and frequency (Z)

load_functions: Functions which load impedances, magnetic data, powerlines, and geographic information

plotting: Functions to reproduce figures.

utils: General functions (e.g. fft, ifft)


