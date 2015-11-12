#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# splash.py
#
# 2014-01-30 -- created
# 2015-11-11 -- last updated
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
# algorithms for simulating habitats (SPLASH): Modelling radiation evapo-
# transpiration and plant-available moisture, Geoscientific Model Development,
# 2015 (in progress)

###############################################################################
## IMPORT MODULES:
###############################################################################
import numpy

from const import kCw, kWm
from evap import EVAP


###############################################################################
## CLASSES
###############################################################################
class SPLASH:
    """
    Name:     SPLASH
    Features: This class updates daily quantities of radiation,
              evapotranspiration, soil moisture and runoff based on SPLASH.
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lat, elv):
        """
        Name:     SPLASH.__init__
        Input:    - float, latitude, degrees (lat)
                  - float, elevation, meters (elv)
        """
        # Error handle and assign required public variables:
        self.elv = elv
        #
        if lat > 90.0 or lat < -90.0:
            print "Latitude outside range of validity (-90 to 90)!"
            exit(1)
        else:
            self.lat = lat
        #
        # Initialize daily status variables:
        self.ho = 0.      # daily solar irradiation, J/m2
        self.hn = 0.      # daily net radiation, J/m2
        self.ppfd = 0.    # daily PPFD, mol/m2
        self.cond = 0.    # daily condensation water, mm
        self.wn = 0.      # daily soil moisture, mm
        self.precip = 0.  # daily precipitation, mm
        self.ro = 0.      # daily runoff, mm
        self.eet = 0.     # daily equilibrium ET, mm
        self.pet = 0.     # daily potential ET, mm
        self.aet = 0.     # daily actual ET, mm

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def spin_up(self, d):
        """
        Name:     SPLASH.spin
        Input:    - DATA class, (d)
                  - bool, (to_write)
        Output:   None.
        Features: Spins up the daily soil moisture.
        Depends:  quick_run
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Create a soil moisture array:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if isinstance(d.num_lines, list):
            n = d.num_lines[0]
            if (numpy.array(d.num_lines) == n).all():
                wn_vec = numpy.zeros((n,))
            else:
                print "Invalid number of lines read from DATA class!"
        else:
            n = d.num_lines
            wn_vec = numpy.zeros((n,))
        print "Created soil moisture array of length", len(wn_vec)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Run one year:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for i in xrange(n):
            # Get preceding soil moisture status:
            if i == 0:
                wn = wn_vec[-1]
            else:
                wn = wn_vec[i-1]
            #
            # Calculate soil moisture and runoff:
            sm, ro = self.quick_run(n=i+1,
                                    y=d.year,
                                    wn=wn,
                                    sf=d.sf_vec[i],
                                    tc=d.tair_vec[i],
                                    pn=d.pn_vec[i])
            wn_vec[i] = sm
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate change in starting soil moisture:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        start_sm = wn_vec[0]
        end_sm, ro = self.quick_run(n=1,
                                    y=d.year,
                                    wn=wn_vec[-1],
                                    sf=d.sf_vec[0],
                                    tc=d.tair_vec[0],
                                    pn=d.pn_vec[0])
        diff_sm = numpy.abs(end_sm - start_sm)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Equilibrate:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        spin_count = 1
        while diff_sm > 1.0:
            for i in xrange(n):
                # Get preceding soil moisture status:
                if i == 0:
                    wn = wn_vec[-1]
                else:
                    wn = wn_vec[i-1]
                #
                # Calculate soil moisture and runoff:
                sm, ro = self.quick_run(n=i+1,
                                        y=d.year,
                                        wn=wn,
                                        sf=d.sf_vec[i],
                                        tc=d.tair_vec[i],
                                        pn=d.pn_vec[i])
                wn_vec[i] = sm
            #
            start_sm = wn_vec[0]
            end_sm, ro = self.quick_run(n=1,
                                        y=d.year,
                                        wn=wn_vec[-1],
                                        sf=d.sf_vec[0],
                                        tc=d.tair_vec[0],
                                        pn=d.pn_vec[0])
            diff_sm = numpy.abs(end_sm - start_sm)
            spin_count += 1
        #
        print "Spun", spin_count, "years"
        self.wn_vec = wn_vec
        self.wn = wn_vec[-1]

    def quick_run(self, n, y, wn, sf, tc, pn):
        """
        Name:     SPLASH.quick_run
        Inputs:   - int, day of year (n)
                  - int, year (y)
                  - float, daily soil water content, mm (wn)
                  - float, daily fraction of bright sunshine (sf)
                  - float, daily air temperature, deg C (tc)
                  - float, daily precipitation, mm (pn)
        Output:   - float, soil moisture, mm (sm)
                  - float, runoff, mm (ro)
        Features: Returns daily soil moisture and runoff.
        Depends:  - kCw
                  - kWm
                  - EVAP
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate evaporative supply rate, mm/h
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sw = kCw*(wn/kWm)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate radiation and evaporation quantities
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        my_evap = EVAP(self.lat, n, self.elv, y, sf, tc, sw)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate today's soil moisture, mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sm = wn + pn + my_evap.cond - my_evap.aet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate runoff, mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if sm > kWm:
            # Bucket is too full
            #   allocate excess water to runoff
            #   set soil moisture to capacity (i.e., kWm)
            ro = sm - kWm
            sm = kWm
        elif sm < 0:
            # Bucket is too empty
            #   set soil moisture and runoff to zero
            sm = 0
            ro = 0
        else:
            ro = 0
        #
        return(sm, ro)

    def run_one_day(self, n, y, wn, sf, tc, pn):
        """
        Name:     SPLASH.run_one_day
        Inputs:   - int, day of year (n)
                  - int, year (y)
                  - float, soil water content, mm (wn)
                  - float, fraction of bright sunshine (sf)
                  - float, air temperature, deg C (tc)
                  - float, precipitation, mm (pn)
        Outputs:  None
        Features: Runs SPLASH model for one day.
        Depends:  - kCw
                  - kWm
                  - EVAP
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Set meteorological variables:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.precip = pn    # daily precipitation, mm
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate evaporative supply rate (sw), mm/h
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sw = kCw*(wn/kWm)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate radiation and evaporation quantities
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        my_evap = EVAP(self.lat, n, self.elv, y, sf, tc, sw)
        self.ho = my_evap.ra_d      # daily solar irradiation, J/m2
        self.hn = my_evap.rn_d      # daily net radiation, J/m2
        self.ppfd = my_evap.ppfd_d  # daily PPFD, mol/m2
        self.cond = my_evap.cond    # daily condensation water, mm
        self.eet = my_evap.eet_d    # daily equilibrium ET, mm
        self.pet = my_evap.pet_d    # daily potential ET, mm
        self.aet = my_evap.aet_d    # daily actual ET, mm
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate today's soil moisture (sm), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sm = wn + pn + my_evap.cond - my_evap.aet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate runoff (ro), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if sm > kWm:
            # Bucket is too full
            #   allocate excess water to runoff
            #   set soil moisture to capacity (i.e., kWm)
            ro = sm - kWm
            sm = kWm
        elif sm < 0:
            # Bucket is too empty
            #   reduce actual ET by discrepancy amount
            #   set soil moisture and runoff to zero
            self.aet += sm
            sm = 0
            ro = 0
        else:
            ro = 0
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Update soil moisture & runoff
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.wn = sm  # daily soil moisture, mm
        self.ro = ro  # daily runoff, mm

    def print_vals(self):
        """
        Name:     SPLASH.print_vals
        Inputs:   None
        Outputs:  None
        Features: Prints all daily values
        """
        print "Daily values:"
        print "  Ho: %0.6f  MJ/m^2" % ((1.0e-6)*self.ho)
        print "  HN: %0.6f MJ/m^2" % ((1.0e-6)*self.hn)
        print "  PAR: %0.6f mol/m^2" % (self.ppfd)
        print "  Cn: %0.6f mm" % (self.cond)
        print "  EET: %0.6f mm" % (self.eet)
        print "  PET: %0.6f mm" % (self.pet)
        print "  AET: %0.6f mm" % (self.aet)
        print "  Wn: %0.6f mm" % (self.wn)
        print "  RO: %0.6f mm" % (self.ro)
