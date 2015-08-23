# README
---------------

* LAST UPDATED: 2015-08-23
* TEAM: labprentice
* REPO: stash (private)

## Contents
--------------------
### data/
This directory holds example data files (CSV and TXT) as well as a script for producing additional input data for SPLASH.

* __stash_getdata.py__
    * Script that produces a CSV with daily input data (i.e., sunshine fraction, air temperature, and precipitation)
* __example_data.csv__
    * Example comma separated daily data for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)
* __daily_pn_2000_wfdei.txt__
    * Example daily precipitation data from the WATCH Forcing Data ERA Interim for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)
* __daily_sf_2000_cruts.txt__
    * Example daily fraction of bright sunshine hours based on CRU Time Series cloudiness for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)
* __daily_tair_2000_wfdei.txt__
    * Example daily air temperature from the WATCH Forcing Data ERA Interim for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)

### doc/
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

### cpp_version/
This directory holds the C++ version of the SPLASH code.

* __DATA.cpp__
    * C++ class definition file. This class handles reading input data.

* __DATA.h__
    * C++ header file for DATA class.

* __etr.h__
    * C++ header file for the __e__vapo__t__ranspiration and __r__adiation structure.

* __EVAP.cpp__
    * C++ class definition file. This class calculates daily quantities of radiation, evaporation, and condensation.

* __EVAP.h__
    * C++ header file for EVAP class.

* __main.cpp__
    * Main function for running the SPLASH code.

* __Makefile__
    * Makefile for C++ SPLASH code.

* __smr.h__
    * C++ header file for the __s__oil __m__oisture and __r__unoff structure.

* __SPLASH.cpp__
    * C++ class definition file. This class updates daily quantities of radiation, evapotranspiration, soil moisture and runoff based on the SPLASH methodology.

* __SPLASH.h__
    * C++ header file for SPLASH class.

### f90_version/
This directory holds the FORTRAN90 version of the SPLASH code. 

* __Makefile__ 
    * Use to compile the splash.F script.

* __splash.F__ 
    * Runs SPLASH for one year, following a spin-up of soil moisture, based on example monthly meteorological data (hard-coded).

### py_version/
This directory holds the Python version of the SPLASH code. 

* __const.py__
    * Contains definitions for SPLASH global constants.

* __data.py__
    * DATA class definition for reading input data.

* __evap.py__
    * EVAP class definition for calculating daily radiation and evapotranspiration quantities:
        * Photosynthetic Photon Flux Density, mol/m^2/day
        * Equilibrium ET, mm/day
        * Potential ET, mm/day
        * Actual ET, mm/day
        * Condensation, mm/day

* __main.py__
    * Main function for running the SPLASH code.

* __splash.py__ 
    * SPLASH class definition for updating daily quantities of radiation, evapotranspiration, soil moisture and runoff.

* __stash_grid.py__ 
    * __NOT FULLY TESTED!__
    * Implements the EVAP_G class for grid-based processing 
    * Inputs include:
        * day of year
        * elevation (360x720 array), meters
        * sunshine fraction (360x720 array), decimal
        * mean daily air temperature (360x720 array), °C
        * evaporative supply rate (360x720 array), mm/h
        * year (optional)
    * CRU-based input data is used (user must have a copy of data files and specify their location)

### r_version/
This directory holds the R version of the SPLASH code. 

* __splash.R__ 
    * Function-based implementation of the SPLASH code, including:
        * reading input data from file
        * spin-up of soil moisture
        * running SPLASH at the daily time scale
    * Includes plotting examples of monthly and daily results


## Simple Process-Led Algorithms for Simulating Habitats (SPLASH): Modelling Radiation, Evapotranspiration and Plant-Available Moisture
----------------------------------------------------------------------------
### Theory
There is a growing need of global ecophysiological datasets for the study of vegetation dynamics under changing climate scenarios; however, simulation of natural processes is often necessary due to the lack of observations.
Bioclimatic indices, such as the climatic water deficit and the plant available water coefficient, are improvements over indices of mean annual temperature and precipitation. The algorithms to produce these indices are based on the STASH (STAtic SHell) model, developed as a simple process-based predictive model for the simulation of tree species distributions at the regional scale (Sykes and Prentice, 1995, 1996; Sykes et al., 1996). 

In this work, we update, correct and improve the mechanistic processes of the STASH model, now entitled SPLASH (Simple Process-Led Algorithms for Simulating Habitats), to create simple, generic and robust algorithms of analytical formulae, developed under a set of practical assumptions and simplifications, to provide key modeling parameters (e.g., radiation, evapotranspiration, and soil moisture) at relevant time scales. 
The model, currently designed to run for a specific location (i.e., latitude and elevation), operates on a minimum of three meteorological inputs including: precipitation (mm), air temperature (Celsius), and fraction of bright sunshine (unitless).

The methodology follows the steps outlined in Cramer & Prentice (1988) where daily soil moisture (*Wn*) is calculated based on the previous day's moisture content, incremented by daily precipitation (*Pn*) and condensation (*Cn*), and reduced by the daily actual evapotranspiration (*Ea*):

![equation](http://www.sciweavers.org/tex2img.php?eq=W_n%3DW_%7Bn-1%7D%2BP_n%2BC_n-E%5Ea_n&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

To solve the simple bucket model presented above, the following steps are taken at the daily timescale: calculate the radiation terms, estimate the condensation, estimate the evaporative supply, estimate the evaporative demand, calculate the actual evapotranspiration, and update the daily soil moisture. 
At the end of each month, daily quantities may be aggregated into monthly totals and additional moisture indexes may be calculated.

### Key Outputs
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
* Monthly
    * Photosynthetic photon flux density
    * Equilibrium evapotranspiration
    * Potential evapotranspiration
    * Cramer-Prentice bioclimatic moisture index
    * Climatic water deficit

### Model Inputs
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

## References
--------------------
* Cramer, W. and I. C. Prentice (1988) Simulation of regional soil moisture deficits on a European scale, _Norsk Geografisk Tidsskrift - Norwegian Journal of Geography_, 42:2-3, pp. 149-151.
* Harris, I., P. D. Jones, T. J. Osborn, and D. H. Lister (2014) Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset, _Int. J. Climatol._, 34, 623–642, doi:10.1002/joc.3711.
* Sykes, M. T. and I. C. Prentice (1995) Boreal forest futures: modelling the controls on tree species range limits and transient responses to climate change, _Water, Air, Soil Pollut._, 82, 415-428.
* Sykes, M. T. and I. C. Prentice (1996) Climate change, tree species distributions and forest dynamics: a case study in the mixed conifer/northern hardwoods zone in Northern Europe, _Clim. Change_, 34, 161-177.
* Sykes, M. T., I. C. Prentice, and W. Cramer (1996) A bioclimatic model for the potential distributions of north European tree species under present and future climates, _J. Biogeogr._, 23(2), 203-233.
* Weedon, G. P., G. Balsamo, N. Bellouin, S. Gomes, M. J. Best, and P. Viterbo (2014) The WFDEI meteorological forcing data set: WATCH Forcing Data methodology applied to ERA-Interim reanalysis data, _Water Resour. Res._, 50, doi:10.1002/2014WR015638.