# Table of Contents
-------------------
[TOC]

# Repository Details
---------------

* LATEST RELEASE: v1.0
* LAST UPDATED: 2016-02-16
* TEAM: labprentice
* REPO: splash (public)

# Repository Structure
----------------------
## data/
This directory holds example data files (CSV and TXT). Note that the __py_verion__ directory contains a script for producing additional input data for SPLASH.

* __example_data.csv__

    * Example comma separated daily data for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)

* __daily_pn_2000_wfdei.txt__

    * Example daily precipitation data from the WATCH Forcing Data ERA Interim for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)

* __daily_sf_2000_cruts.txt__

    * Example daily fraction of bright sunshine hours based on CRU Time Series cloudiness for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)

* __daily_tair_2000_wfdei.txt__

    * Example daily air temperature from the WATCH Forcing Data ERA Interim for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)

## doc/
This directory holds the current documentation for the SPLASH code.

* __splash_doc.pdf__

    * The current PDF build of the documentation

* __splash_doc.tex__

    * The main LaTeX document file

* __splash.bib__

    * The BibLatex file for documentation references

* __img/__

    * Contains the EPS figures for the documentation

* __tex/__

    * Contains the modular LaTeX chapter files and appendix

## releases/v1.0/

### cpp_version/
This directory holds the C++ version of the SPLASH code.

* __DATA.cpp / DATA.h__

    * C++ class definition file. This class handles reading input data.

    * Utilizes *cstdlib*, *fstream*, *iostream*, *string*, and *vector* libraries

* __etr.h__

    * C++ header file for the __e__vapo__t__ranspiration and __r__adiation structure.

* __EVAP.cpp / EVAP.h__

    * C++ class definition file. This class calculates daily quantities of radiation, evaporation, and condensation.

    * Utilizes *cmath*, *stdio.h*, and *vector* libraries

* __global.cpp / global.h__

    * Script defining global constants.

* __main.cpp__

    * Main function for running the SPLASH code.

    * Utilizes *iostream*, *string* and *vector* libraries

* __Makefile__

    * Makefile for GNU C++ compiler.

* __smr.h__

    * C++ header file for the __s__oil __m__oisture and __r__unoff structure.

* __SOLAR.cpp / SOLAR.h__

    * C++ class definition file. This class calculates daily radiation fluxes.

* __SPLASH.cpp / SPLASH.h__

    * C++ class definition file. This class updates daily quantities of radiation, evapotranspiration, soil moisture and runoff based on the SPLASH methodology.

    * Utilizes *cmath*, *stdio.h* and *vector* libraries

* __srad.h__

    * C++ header file for solar radiation fluxes structure.

### f90_version/
This directory holds the FORTRAN90 version of the SPLASH code.

* __Makefile__

    * Makefile for gfortran compiler.

* __spinup_splash.f90__

    * Source code for

* __splash.mod.f90__

    * Source code for the SPLASH model to run for one year, following a spin-up of soil moisture, based on example meteorological data.

* __test_splash.f90__

    * Source code for consistency tests 1--3 (see Wiki).

### py_version/
This directory holds the Python version of the SPLASH code.

* __const.py__

    * Contains definitions for SPLASH global constants.

    * Utilizes *numpy* library

* __data.py__

    * DATA class definition for reading input data.

    * Utilizes *logging* and *numpy* libraries

* __evap.py__

    * EVAP class definition for calculating daily evapotranspiration quantities:

        * Equilibrium ET, mm/day

        * Potential ET, mm/day

        * Actual ET, mm/day

        * Condensation, mm/day

   * Utilizes *logging* and *numpy* libraries

* __main.py__

    * Main function for running the SPLASH code.

    * Utilizes *logging* library

* __solar.py__

    * SOLAR class definition for calculating daily radiation fluxes.

    * Utilizes *logging* and *numpy* libraries

* __splash.py__

    * SPLASH class definition for updating daily quantities of radiation, evapotranspiration, soil moisture and runoff.

    * Utilities *logging* and *numpy* libraries

* __splash_data.py__

    * SPLASH_DATA class and script to produce a CSV file with daily input data (i.e., sunshine fraction, air temperature, and precipitation)

    * Utilities *datetime*, *glob*, *numpy*, *os.path* and *scipy.io.netcdf* libraries

* __utilities.py__

    * Contains utility functions that are shared amongst classes.

    * Utilizes *logging* and *numpy* libraries

### r_version/
This directory holds the R version of the SPLASH code.

* __const.R__

    * This script contains the global constants defined in SPLASH.

* __data.R__

    * This script contains functions to handle the file IO for reading and writing data, i.e.:

        * read_csv(character fname, double y=-1)

        * read_txt(list my_data, character fname, character var, double y=-1)

* __evap.R__

    * This script contains functions to calculate daily condensation, and evapotranspiration.

* __main.R__

    * This script runs the SPLASH model for one year.

* __solar.R__

    * This script contains functions to calculate daily radiation fluxes.

* __splash.R__

    * This script contains functions for running SPLASH for point-based data, i.e.:

        * spin_up(list mdat, list dtot)

        * quick_run(double lat, double elv, double n, double y, double wn, double sf, double tc, double pn)

        * run_one_day(double lat, double elv, double n, double y, double wn, double sf, double tc, double pn)

* __test.R__

    * This script performs SPLASH consistency tests.

## working/
This directory contains the source code currently under development.

# SPLASH: Robust indices of radiation, evapotranspiration and plant-available moisture
----------------------------------------------------------------------------
## Theory
There is a growing need of global ecophysiological datasets for the study of vegetation dynamics under changing climate scenarios; however, simulation of natural processes is often necessary due to the lack of observations.
Bioclimatic indices, such as the climatic water deficit and the plant available water coefficient, are improvements over indices of mean annual temperature and precipitation. The algorithms to produce these indices are based on the STASH (STAtic SHell) model, developed as a simple process-based predictive model for the simulation of tree species distributions at the regional scale (Sykes and Prentice, 1995, 1996; Sykes et al., 1996).

In this work, we update, correct and improve the mechanistic processes of the STASH model, now entitled SPLASH (Simple Process-Led Algorithms for Simulating Habitats), to create simple, generic and robust algorithms of analytical formulae, developed under a set of practical assumptions and simplifications, to provide key modeling parameters (e.g., radiation, evapotranspiration, and soil moisture) at relevant time scales.
The model, currently designed to run for a specific location (i.e., latitude and elevation), operates on a minimum of three meteorological inputs including: precipitation (mm), air temperature (Celsius), and fraction of bright sunshine (unitless).

The methodology follows the steps outlined in Cramer & Prentice (1988) where daily soil moisture (*Wn*) is calculated based on the previous day's moisture content, incremented by daily precipitation (*Pn*) and condensation (*Cn*), and reduced by the daily actual evapotranspiration (*Ea*):

![equation](http://www.sciweavers.org/tex2img.php?eq=W_n%3DW_%7Bn-1%7D%2BP_n%2BC_n-E%5Ea_n&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

To solve the simple bucket model presented above, the following steps are taken at the daily timescale: calculate the radiation terms, estimate the condensation, estimate the evaporative supply, estimate the evaporative demand, calculate the actual evapotranspiration, and update the daily soil moisture.
At the end of each month, daily quantities may be aggregated into monthly totals and additional moisture indexes may be calculated.

## Key Outputs
* Daily
    * Extraterrestrial solar irradiation
    * Net surface radiation
    * Photosynthetic photon flux density
    * Condensation
    * Soil moisture
    * Runoff
    * Equilibrium evapotranspiration
    * Potential evapotranspiration
    * Actual evapotranspiration

## Model Inputs
For radiation, the basic geographic coordinates and time parameters needed are:

* year (*y*)
* day of year (*n*)
* latitude (φ), radians
* elevation (z), meters

For evapotranspiration, the basic meteorological variables needed are:

* daily air temperature (*Tair*), °C
* daily precipitation (*Pn*), mm
* daily fraction of bright sunshine hours (*Sf*), %

For spatial analyses, the 0.5° x 0.5° gridded [CRU TS 3.21](http://badc.nerc.ac.uk/view/badc.nerc.ac.uk__ATOM__ACTIVITY_0c08abfc-f2d5-11e2-a948-00163e251233) data sets may be used (Harris et al., 2014), for example: TMP (monthly mean daily air temperature); PRE (monthly precipitation totals); CLD (cloudiness fraction); and CRU TS 3.0 ELV (mean pixel elevation). Daily precipitation and air temperature are also available from the [WATCH](http://www.eu-watch.org/data_availability) dataset (Weedon et al., 2014).

# References
--------------------
* Cramer, W. and I. C. Prentice (1988) Simulation of regional soil moisture deficits on a European scale, _Norsk Geografisk Tidsskrift - Norwegian Journal of Geography_, 42:2-3, pp. 149-151.
* Harris, I., P. D. Jones, T. J. Osborn, and D. H. Lister (2014) Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset, _Int. J. Climatol._, 34, 623–642, doi:10.1002/joc.3711.
* Sykes, M. T. and I. C. Prentice (1995) Boreal forest futures: modelling the controls on tree species range limits and transient responses to climate change, _Water, Air, Soil Pollut._, 82, 415-428.
* Sykes, M. T. and I. C. Prentice (1996) Climate change, tree species distributions and forest dynamics: a case study in the mixed conifer/northern hardwoods zone in Northern Europe, _Clim. Change_, 34, 161-177.
* Sykes, M. T., I. C. Prentice, and W. Cramer (1996) A bioclimatic model for the potential distributions of north European tree species under present and future climates, _J. Biogeogr._, 23(2), 203-233.
* Weedon, G. P., G. Balsamo, N. Bellouin, S. Gomes, M. J. Best, and P. Viterbo (2014) The WFDEI meteorological forcing data set: WATCH Forcing Data methodology applied to ERA-Interim reanalysis data, _Water Resour. Res._, 50, doi:10.1002/2014WR015638.
