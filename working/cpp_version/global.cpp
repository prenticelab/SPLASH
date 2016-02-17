#include "global.h"

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * global.cpp
 *
 * 2016-01-22 -- created
 * 2016-01-22 -- last updated
 *
 * ~~~~~~~~~
 * citation:
 * ~~~~~~~~~
 * T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
 * Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
 * algorithms for simulating habitats (SPLASH): Robust indices of radiation
 * evapotranspiration and plant-available moisture, Geoscientific Model
 * Development, 2016 (in progress)
 *
 * ~~~~~~~~~~~~
 * description:
 * ~~~~~~~~~~~~
 * This script holds global constants for the SPLASH model.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 *
 *
 * //////////////////////////////////////////////////////////////////////// */
namespace Global {
    extern const double A(107.0);         // (Monteith & Unsworth, 1990)
    extern const double alb_sw(0.17);     // (Federer, 1968)
    extern const double alb_vis(0.03);    // (Sellers, 1985)
    extern const double b(0.20);          // (Linacre, 1968)
    extern const double c(0.25);          // (Linacre, 1968)
    extern const double d(0.50);          // (Linacre, 1968)
    extern const double e(0.0167);        // (Berger, 1978)
    extern const double eps(23.44);       // (Berger, 1978)
    extern const double fFEC(2.04);       // (Meek et al., 1984)
    extern const double G(9.80665);       // (Allen, 1973)
    extern const double Gsc(1360.8);      // (Kopp & Lean, 2011)
    extern const double L(0.0065);        // (Allen, 1973)
    extern const double Ma(0.028963);     // (Tsilingiris, 2008)
    extern const double Mv(0.01802);      // (Tsilingiris, 2008)
    extern const double Po(101325);       // (Allen, 1973)
    extern const double R(8.31447);       // (Moldover et al., 1988)
    extern const double To(288.15);       // (Berberan-Santos et al., 1997)
    extern const double w(0.26);          // (Priestley & Taylor, 1972)
    extern const double omega(283.0);     // (Berger, 1978)
    extern const double Cw(1.05);         // (Federer, 1982)
    extern const double Wm(150.0);        // (Cramer & Prentice, 1988)
    extern const double PI(3.141592653589793);  // pi
    extern const double pir = (PI/180.0); // degrees to radians
}
