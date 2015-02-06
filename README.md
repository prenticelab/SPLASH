# README
---------------

* LAST UPDATED: 2015-02-06
* TEAM: labprentice
* REPO: stash (private)

## Contents
--------------------
### data/
This directory holds a script for producing input data for the STASH 2.0 code.

* __stash_getdata.py__
    * Script that produces a CSV with daily input data (i.e., sunshine fraction, air temperature, and precipitation)
* __example_data.csv__
    * Example daily data for San Francisco, United States (37.7 N, 122.4 W, 142 m, 2000 CE)

### doc/
This directory holds the current documentation for the STASH 2.0 code.

* __stash_doc.pdf__
    * The current PDF build of the documentation
* __stash_doc.tex__
    * The main LaTeX document file
* __stash.bib__
    * The BibLatex file for documentation references
* __img/__
    * Contains the EPS figures for the documentation
* __tex/__
    * Contains the modular LaTeX chapter files (and appendix)

### cpp_version/
This directory holds the C++ version of the STASH 2.0 code.

* __EVAP.cpp__
    * C++ class definition file.

### f90_version/
This directory holds the FORTRAN90 version of the STASH 2.0 code. 

* __Makefile__ 
    * Use to compile the stash.F script.

* __stash.F__ 
    * TBA.

### py_version/
This directory holds the Python version of the STASH 2.0 code. 

* __stash.py__ 
    * Implements the EVAP class for point-based processing 
    * Inputs include:
        * latitude, degrees
        * day of year
        * elevation (optional), meters
        * year (optional)
        * sunshine fraction (optional), decimal
        * mean daily air temperature (optional), °C
        * evaporative supply rate (optional), mm/h
    * Input data must be imported separately by user (example data is available in the script).

* __stash_grid.py__ 
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
This directory holds the R version of the STASH 2.0 code. 

* __stash.R__ 
    * Implements the EVAP function for point-based processing 
    * Inputs include:
        * latitude, degrees
        * day of year
        * elevation (optional), meters
        * year (optional)
        * sunshine fraction (optional), decimal
        * mean daily air temperature (optional), °C
        * evaporative supply rate (optional), mm/h
    * Input data must be imported separately by user (example data is available in the script)
    * Includes plotting examples of monthly and daily results


## STASH 2.0: Radiation, Evapotranspiration and Moisture Availability
----------------------------------------------------------------------------
### Theory
With a growing need of global ecophysiological datasets for the study of vegetation dynamics under changing climate scenarios, the simple process-based bioclimatic model (STASH) presents an algorithm of analytical formulae, under a set of practical assumptions and simplifications, that provides key modeling parameters (e.g., radiation, evapotranspiration, and soil moisture) at the daily time scale. The model, currently designed to run for a specific location (i.e., latitude and elevation), operates on a minimum of three meteorological inputs including: precipitation (mm), air temperature (Celsius), and fraction of bright sunshine (unitless).

The methodology follows the pseudo-code presented by Cramer & Prentice (1988) where daily soil moisture (*Wn*) is calculated based on the previous day's moisture content, incremented by daily precipitation (*Pn*) and condensation (*Cn*), and reduced by the daily actual evapotranspiration (*Ea*):

![equation](http://www.sciweavers.org/tex2img.php?eq=W_n%3DW_%7Bn-1%7D%2BP_n%2BC_n-E%5Ea_n&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

To solve the simple bucket model presented above, the following steps are taken at the daily timescale: calculate the radiation terms, estimate the condensation, estimate the evaporative supply, estimate the evaporative demand, calculate the actual evapotranspiration, and update the daily soil moisture. At the end of each month, daily quantities may be aggregated into monthly totals and additional moisture indexes may be calculated.

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
* Weedon, G. P., G. Balsamo, N. Bellouin, S. Gomes, M. J. Best, and P. Viterbo (2014) The WFDEI meteorological forcing data set: WATCH Forcing Data methodology applied to ERA-Interim reanalysis data, _Water Resour. Res._, 50, doi:10.1002/2014WR015638.