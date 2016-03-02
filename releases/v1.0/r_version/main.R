# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# main.R
#
# VERSION: 1.0
# LAST UPDATED: 2016-02-19
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
# This script runs the SPLASH model for one year.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - updated monthly and daily results & process [15.01.13]
# - updated plots of results [15.01.16]
# - added example data CSV file [15.01.16]
# - fixed Cramer-Prentice alpha definition [15.01.16]
# - updated monthly results plot [15.01.22]
# - added write out for daily/monthly results [15.01.27]
# - added example of yearly looping [15.01.27]
# - updated year extraction from filename in yearly loop example [15.01.30]
# - fixed plots for Figs. 3 & 4 of manuscript [15.11.23]
# - fixed directory paths [16.02.17]
#
#### IMPORT SOURCES ##########################################################
source("const.R")
source("data.R")
source("evap.R")
source("splash.R")


#### SPLASH ###################################################################
# Initialize daily results:
daily_totals <- matrix(data=rep(0, 3294), nrow=366, ncol=9)
daily_totals <- as.data.frame(daily_totals)
names(daily_totals) <- c("ho",   # daily solar irradiation, J/m2
                         "hn",   # daily net radiation, J/m2
                         "qn",   # daily PPFD, mol/m2
                         "cn",   # daily condensation, mm
                         "wn",   # daily soil moisture, mm
                         "ro",   # daily runoff, mm
                         "eq_n", # daily equilibrium ET, mm
                         "ep_n", # daily potential ET, mm
                         "ea_n") # daily actual ET, mm

# Initialize monthly results:
monthly_totals <- matrix(data=rep(0, 72), nrow=12, ncol=6)
monthly_totals <- as.data.frame(monthly_totals)
names(monthly_totals) <- c("eq_m",  # monthly equilibrium ET, mm
                           "ep_m",  # monthly potential ET, mm
                           "ea_m",  # monthly actual ET, mm
                           "cpa",   # Cramer-Prentice alpha, unitless
                           "cwd",   # climatic water deficit, mm
                           "q_m")   # monthly PPFD, mol/m2

# Location constants:
my_lat <- 37.7
my_elv <- 142

# Calculate days in the year
y <- 2000
ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)

# Example data (San Francisco, 2000 CE)
#   $file_name, file name
#   $sf, fractional sunshine hours, unitless
#   $tair, air temperature, deg C
#   $pn, daily precipitation, mm
#   $num_lines, number of lines of data
#   $year, year (Gregorian calendar)
my_file <- "../../../data/example_data.csv"
my_data <- read_csv(my_file, 2000)
my_data$lat_deg <- my_lat
my_data$elv_m <- my_elv

# Spin up the soil moisture content
daily_totals <- spin_up(my_data, daily_totals)

# Run one day:
daily_vals <- run_one_day(my_data$lat_deg,
                          my_data$elv_m,
                          172,
                          my_data$year,
                          145.0,
                          0.5,
                          17.3,
                          10.0)

cat(paste("Daily values:\n",
          "  Ho:", (1.0e-6)*daily_vals$ho, " MJ/m^2\n",
          "  HN:", (1.0e-6)*daily_vals$hn, " MJ/m^2\n",
          "  PAR:", daily_vals$ppfd, " mol/m^2\n",
          "  Cn:", daily_vals$cond, " mm\n",
          "  EET:", daily_vals$eet, " mm\n",
          "  PET:", daily_vals$pet, " mm\n",
          "  AET:", daily_vals$aet, " mm\n",
          "  Wn:", daily_vals$wn, " mm\n",
          "  RO:", daily_vals$ro, " mm"
          )
    )

# Run SPLASH for a full year:
all_months <- seq(from=1, to=12, by=1)
monthly_totals <- monthly_totals*0

# monthly
for (m in all_months) {
    # Calculate days of current month:
    nm <- julian_day(y, m + 1, 1) - julian_day(y, m, 1)

    # daily:
    for (i in seq(from=1, to=nm, by=1)) {
        # Calculate day of year:
        n <- julian_day(y, m, i) - julian_day(y, 1 , 1) + 1

        idx <- (n - 1)
        if (idx < 1) {
            idx <- ny
        }
        daily_vals <- run_one_day(my_data$lat_deg,
                                  my_data$elv_m,
                                  n,
                                  my_data$year,
                                  daily_totals$wn[idx],
                                  my_data$sf[n],
                                  my_data$tair[n],
                                  my_data$pn[n])

        # Update daily values:
        daily_totals$wn[n] <- daily_vals$wn
        daily_totals$ro[n] <- daily_vals$ro

        if ( (ny == 365) & (n == 365) ) {
            daily_totals$wn[n + 1] <- daily_totals$wn[n]
        }

        # Save daily results:
        daily_totals$ho[n] <- daily_vals$ho
        daily_totals$hn[n] <- daily_vals$hn
        daily_totals$qn[n] <- daily_vals$ppfd
        daily_totals$cn[n] <- daily_vals$cond
        daily_totals$eq_n[n] <- daily_vals$eet
        daily_totals$ep_n[n] <- daily_vals$pet
        daily_totals$ea_n[n] <- daily_vals$aet

        # Update monthly totals:
        monthly_totals$eq_m[m] <- monthly_totals$eq_m[m] + daily_vals$eet
        monthly_totals$ep_m[m] <- monthly_totals$ep_m[m] + daily_vals$pet
        monthly_totals$ea_m[m] <- monthly_totals$ea_m[m] + daily_vals$aet
        monthly_totals$q_m[m] <- monthly_totals$q_m[m] + daily_vals$ppfd
    } # end daily
    monthly_totals$cpa[m] <- monthly_totals$ea_m[m]/monthly_totals$eq_m[m]
    monthly_totals$cwd[m] <- monthly_totals$ep_m[m] - monthly_totals$ea_m[m]
} # end monthly

# Save results
daily_outfile <- "out/splash_results_daily.csv"
write.csv(daily_totals, file=daily_outfile)

monthly_outfile <- "out/splash_results_monthly.csv"
write.csv(monthly_totals, file=monthly_outfile)

# View results
to_plot <- FALSE
if (to_plot) {
    export_fig <- FALSE

    # Plot monthly ET results
    par(mfrow=c(1, 1))
    par(mar=c(4, 4.5, 1, 1))
    plot(monthly_totals$ep_m,type="l",lwd=2, col="black",
       ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
    lines(monthly_totals$ea_m, lty=2, lwd=2, col="cyan")
    lines(monthly_totals$eq_m, lty=1, lwd=2, col="green")
    lines(monthly_totals$cwd, lty=3, lwd=2, col="red")
    box(lwd=2)
    axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=1, to=12, by=1))
    axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=1, to=12, by=1))
    axis(side=2, las=1, tck=-0.02, labels=NA)
    axis(side=2, las=1, lwd=0, line=-0.4)
    mtext(side=1, "Month", line=2)
    mtext(side=2, as.expression(ET~(mm)), line=3)
    legend("topleft",legend=c("PET", "AET", "EET", "CWD"),
         col=c("black", "cyan", "green", "red"),
         lty=c(1, 2, 1, 3), lwd=c(2, 2, 2, 2), inset=0.01,
         y.intersp=1.2, horiz=FALSE, bty="n")

    # Plot monthly results, Figure 4
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

    # Plot daily ET results
    par(mfrow=c(1,1))
    par(mar=c(2.5,4.5,1,1))
    plot(daily_totals$ep_n, type="l", lwd=2, col="blue",
         ylim=c(0,1.25*max(daily_totals$ep_n)),
         xlab=NA, ylab=NA, axes=F)
    lines(daily_totals$eq_n, lty=1, lwd=2, col="green")
    lines(daily_totals$ea_n, lty=2, lwd=2, col="red")
    box(lwd=2)
    axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=360, by=60))
    axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=0, to=360, by=60),
         cex.axis=1.2)
    axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=8, by=1))
    axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=0, to=8, by=1),
         cex.axis=1.2)
    mtext(side=2, expression(list(Evapotranspiration, mm~d^{-1})), line=2,
          cex=1.2)
    legend("top",legend=c("Potential",
                          "Equilibrium   ",
                          "Actual"),
           col=c("blue", "green", "red"),
           lty=c(1,1,2), lwd=c(2,2,2), inset=0.01, x.intersp=1.1,
           y.intersp=2.0, horiz=TRUE, bty="n", cex=1.2)


    # Plot daily results, Figure 3
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
