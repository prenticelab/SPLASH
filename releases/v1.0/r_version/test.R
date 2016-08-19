# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# test.R
#
# VERSION: 1.0-r2
# LAST UPDATED: 2016-08-19
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
#
# This file is part of the SPLASH model.
#
# SPLASH is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# SPLASH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
# algorithms for simulating habitats (SPLASH): Robust indices of radiation
# evapo-transpiration and plant-available moisture, Geoscientific Model
# Development, 2016 (in progress)
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script performs SPLASH consistency tests.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
#
#### IMPORT SOURCES ##########################################################
source("const.R")
source("data.R")
source("evap.R")
source("solar.R")
source("splash.R")


#### DEFINE FUNCTIONS ########################################################
# ************************************************************************
# Name:     plot_fig3
# Inputs:   -
# Returns:  None.
# Features: Plots Fig. 3 for manuscript
# ************************************************************************
plot_fig3 <- function(my_data, daily_totals, export_fig=FALSE,
                      out_file="fig3.tiff") {
    if (export_fig) {
        tiff(out_file, width=900, height=1000, units="px",
             compression="none", pointsize=16, res=72)
    }

    par(mfrow=c(8, 1))

    # [1]
    par(mar=c(1, 5, 1, 1))
    plot(my_data$sf, type="l", lwd=2, xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0.3, to=0.7, by=0.1))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=0.3, to=0.7, by=0.1))
    mtext(side=2, expression(italic(S[f])), line=3, cex=1.1)
    text(-12, 0.65, "(a)", pos=4, cex=1.7)

    # [2]
    par(mar=c(1, 5, 1, 1))
    plot(1e-6*daily_totals$hn, type="l", lwd=2, xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=3, to=18, by=3))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=3, to=18, by=3))
    mtext(side=2, expression(italic(H[N])~(MJ~m^{-2})), line=3, cex=1.1)
    text(-12, 17, "(b)", pos=4, cex=1.7)

    # [3]
    par(mar=c(1, 5, 1, 1))
    plot(daily_totals$cn, type="l", lwd=2, xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0.4, to=0.8, by=0.1))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=0.4, to=0.8, by=0.1))
    mtext(side=2, expression(italic(C[n])~(mm)), line=3, cex=1.1)
    text(-12, 0.75, "(c)", pos=4, cex=1.7)

    # [4]
    par(mar=c(1, 5, 1, 1))
    plot(my_data$pn, type="l", lwd=2, xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=-5, to=25, by=5))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=-5, to=25, by=5))
    mtext(side=2, expression(italic(P[n])~(mm)), line=3, cex=1.1)
    text(-12, 22, "(d)", pos=4, cex=1.7)

    # [5]
    par(mar=c(1, 5, 1, 1))
    plot(daily_totals$wn, type="l", lwd=2, xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0, to=150, by=30))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=0, to=150, by=30))
    mtext(side=2, expression(italic(W[n])~(mm)), line=3, cex=1.1)
    text(-12, 130, "(e)", pos=4, cex=1.7)

    # [6]
    par(mar=c(1, 5, 1, 1))
    plot(daily_totals$ro, type="l", lwd=2, xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=-5, to=20, by=5))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=-5, to=20, by=5))
    mtext(side=2, expression(italic(RO)~(mm)), line=3, cex=1.1)
    text(-12, 17, "(f)", pos=4, cex=1.7)

    # [7]
    par(mar=c(1, 5, 1, 1))
    plot(my_data$tair, type="l", lwd=2, xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0, to=25, by=5))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=0, to=25, by=5))
    mtext(side=2, expression(italic(T[air])~(degree*C)), line=3, cex=1.1)
    text(-12, 23, "(g)", pos=4, cex=1.7)

    # [8]
    par(mar=c(2, 5, 1, 1))
    plot(daily_totals$ep_n, type="l", lwd=2, xlab=NA, ylab=NA, axes=F,
         ylim=c(0, max(daily_totals$ep_n)))
    lines(daily_totals$ea_n, lty=2, lwd=2)
    axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
    axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=-60, to=720, by=60),
         cex.axis=1.6)
    axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=-1.5, to=6, by=1.5))
    axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
         at=seq(from=-1.5, to=6, by=1.5))
    mtext(side=2, expression(italic(E[n])~(mm)), line=3, cex=1.1)
    text(-12, 5, "(h)", pos=4, cex=1.7)

    if (export_fig) {
        dev.off()
    }
}


# ************************************************************************
# Name:     plot_fig4
# Inputs:   -
# Returns:  None.
# Features: Plots Fig. 4 for manuscript
# ************************************************************************
plot_fig4 <- function(monthly_totals, export_fig=FALSE, out_file="fig4.tiff") {
    dev.new()

    if (export_fig) {
        tiff(out_file, width=900, height=500, units="px",
             compression="none", pointsize=18, res=72)
    }

    par(mfrow=c(4,1))

    # [1]
    par(mar=c(1,5,1,1))
    plot(monthly_totals$ep_m, type="l", lwd=2, col="black",
         ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
    lines(monthly_totals$ea_m, lty=2, lwd=2)
    axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
    axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-50, to=175, by=25))
    axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-50, to=200, by=50),
         cex.axis=1.6)
    mtext(side=2, expression(italic(E[m])~(mm)), line=3, cex=1.1)
    legend("topright", legend=c(expression(italic(E[m]^{p})),
                                expression(italic(E[m]^{a}))),
           col=c("black", "black"), lty=c(1, 2), cex=1.6, inset=0.02,
           adj=c(0.5, 0.5), lwd=c(2, 2), horiz=TRUE, bty="n", seg.len=1)
    text(0.6, 150, "(a)", pos=4, cex=1.6)

    # [2]
    par(mar=c(1,5,1,1))
    plot(monthly_totals$cwd,type="l",lwd=2, col="black",
         ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
    axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-50, to=175, by=25))
    axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-50, to=150, by=50),
         cex.axis=1.6)
    mtext(side=2, expression(Delta*italic(E[m])~(mm)), line=3, cex=1.1)
    text(0.6, 150, "(b)", pos=4, cex=1.6)

    # [3]
    par(mar=c(1,5,1,1))
    plot(monthly_totals$eq_m, type="l", lwd=2, col="black",
         ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
    lines(monthly_totals$ea_m, lty=2, lwd=2)
    axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
    axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-50, to=175, by=25))
    axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-50, to=150, by=50),
         cex.axis=1.6)
    mtext(side=2, expression(italic(E[m])~(mm)), line=3, cex=1.1)
    legend("topright", legend=c(expression(italic(E[m]^{q})),
                                expression(italic(E[m]^{a}))),
           col=c("black", "black"), lty=c(1, 2), cex=1.6, inset=0.02,
           adj=c(0.5, 0.5), lwd=c(2, 2), horiz=TRUE, bty="n", seg.len=1)
    text(0.6, 150, "(c)", pos=4, cex=1.6)

    # [4]
    par(mar=c(2,5,1,1))
    plot(monthly_totals$cpa,type="l", lwd=2, col="black",
         ylim=c(0, 1.3), xlab=NA, ylab=NA, axes=F)
    axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
    axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=1, to=12, by=1),
         cex.axis=1.6)
    axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-0.3, to=1.2, by=0.3))
    axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-0.3, to=1.2, by=0.3),
         cex.axis=1.6)
    mtext(side=2, expression(italic(alpha[m])), line=3, cex=1.1)
    text(0.6, 1.1, "(d)", pos=4, cex=1.6)

    if (export_fig) {
        dev.off()
    }
}


# TEST 1: SOLAR ############################################################
solar <- calc_daily_solar(lat=37.7, n=172, elv=142, y=2000, sf=1, tc=23.0)
cat(sprintf("TEST 1---Solar values:\n"))
cat(sprintf("  kn: %d\n", solar$kN))
cat(sprintf("  nu: %0.6f degrees\n", solar$nu_deg))
cat(sprintf("  lambda: %0.6f degrees\n", solar$lambda_deg))
cat(sprintf("  rho: %0.6f\n", solar$rho))
cat(sprintf("  dr: %0.6f\n", solar$dr))
cat(sprintf("  delta: %0.6f degrees\n", solar$delta_deg))
cat(sprintf("  ru: %0.6f\n", solar$ru))
cat(sprintf("  rv: %0.6f\n", solar$rv))
cat(sprintf("  rw: %0.6f\n", solar$rw))
cat(sprintf("  hs: %0.6f degrees\n", solar$hs_deg))
cat(sprintf("  hn: %0.6f degrees\n", solar$hn_deg))
cat(sprintf("  tau_o: %0.6f\n", solar$tau_o))
cat(sprintf("  tau: %0.6f\n", solar$tau))
cat(sprintf("  Qn: %0.6f mol/m^2\n", solar$ppfd_mol.m2))
cat(sprintf("  Rnl: %0.6f w/m^2\n", solar$rnl_w.m2))
cat(sprintf("  Ho: %0.6f MJ/m^2\n", (1.0e-6) * solar$ra_j.m2))
cat(sprintf("  Hn: %0.6f MJ/m^2\n", (1.0e-6) * solar$rn_j.m2))
cat(sprintf("  Hnn: %0.6f MJ/m^2\n", (1.0e-6) * solar$rnn_j.m2))

# TEST 2: EVAP ############################################################
evap <- calc_daily_evap(lat=37.7, n=172, elv=142, y=2000, sf=1, tc=23.0, sw=0.9)
cat(sprintf("TEST 2---Evap values:\n"))
cat(sprintf("  s: %0.6f Pa/K\n", evap$s_pa.k))
cat(sprintf("  Lv: %0.6f MJ/kg\n", (1e-6) * evap$lv_j.kg))
cat(sprintf("  Patm: %0.6f bar\n", (1e-5) * evap$patm_pa))
cat(sprintf("  pw: %0.6f kg/m^3\n", evap$pw_kg.m3))
cat(sprintf("  gamma: %0.6f Pa/K\n", evap$gam_pa.k))
cat(sprintf("  Econ: %0.6f mm^3/J\n", (1e9) * evap$econ_m3.j))
cat(sprintf("  Cn: %0.6f mm\n", evap$cond_mm))
cat(sprintf("  rx: %0.6f\n", evap$rx))
cat(sprintf("  hi: %0.6f degrees\n", evap$hi_deg))
cat(sprintf("  EET: %0.6f mm\n", evap$eet_mm))
cat(sprintf("  PET: %0.6f mm\n", evap$pet_mm))
cat(sprintf("  AET: %0.6f mm\n", evap$aet_mm))

# TEST 3: SOIL MOISTURE (run one day) #########################################
t3 <- run_one_day(lat=37.7, elv=142, n=172, y=2000, wn=75, sf=1, tc=23, pn=5)
cat(sprintf("TEST 3---Soil moisture (run one day):\n"))
cat(sprintf("  Ho: %0.6f J/m2\n", t3$ho))
cat(sprintf("  Hn: %0.6f J/m2\n", t3$hn))
cat(sprintf("  PPFD: %0.6f mol/m2\n", t3$ppfd))
cat(sprintf("  EET: %0.6f mm/d\n", t3$eet))
cat(sprintf("  PET: %0.6f mm/d\n", t3$pet))
cat(sprintf("  AET: %0.6f mm/d\n", t3$aet))
cat(sprintf("  Cn: %0.6f mm/d\n", t3$cond))
cat(sprintf("  Wn: %0.6f mm\n", t3$wn))
cat(sprintf("  RO: %0.6f mm\n", t3$ro))

# # TEST 4:
daily_totals <- matrix(data=rep(0, 366), nrow=366, ncol=1)
daily_totals <- as.data.frame(daily_totals)
names(daily_totals) <- c("wn")
my_lat <- 37.7
my_elv <- 142
my_file <- "../../../data/example_data.csv"
my_data <- read_csv(my_file, 2000)
my_data$lat_deg <- my_lat
my_data$elv_m <- my_elv
daily_totals <- spin_up(my_data, daily_totals)
cat(sprintf("TEST 4---Spin-Up:\n"))
cat(sprintf("Day,Wn (mm)\n"))
for (i in seq(from=1, to=my_data$num_lines, by=1)) {
    cat(sprintf("%d,%0.6f\n", i, daily_totals$wn[i]))
}
