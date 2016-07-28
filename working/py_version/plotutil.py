#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# plotutil.py
#
# VERSION: 1.1-dev
# LAST UPDATED: 2016-07-28
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
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
# led algorithms for simulating habitats (SPLASH): Robust indices of radiation,
# evapotranspiration and plant-available moisture, Geoscientific Model
# Development, 2016 (in progress)

###############################################################################
# IMPORT MODULES
###############################################################################
import matplotlib.pyplot as plt

from data import DATA
from splash import SPLASH


###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':

    example = 1
    my_data = DATA()
    if example == 1:
        # Example 1: read CSV file:
        my_file = '../../data/example_data.csv'
        my_data.read_csv(my_file)
    elif example == 2:
        # Example 2: read TXT files:
        my_sf_file = 'daily_sf_2000_cruts.txt'
        my_pn_file = 'daily_pn_2000_wfdei.txt'
        my_tair_file = 'daily_tair_2000_wfdei.txt'
        my_data.read_txt(my_sf_file, 'sf')
        my_data.read_txt(my_pn_file, 'pn')
        my_data.read_txt(my_tair_file, 'tair')

    # Consistency Test #4: Spin-Up
    my_lat = 37.7
    my_elv = 142.
    my_class = SPLASH(my_lat, my_elv, True)
    my_class.spin_up(my_data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ NEEDS UPDATED ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    #
    # Plot monthly ET results
    #
    my_months = range(1, 13)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
    ax1.plot(my_months, my_class.monthly_totals['ep_m'],
             'k-', linewidth=2, label='PET')
    ax1.plot(my_months, my_class.monthly_totals['ea_m'],
             'c--', linewidth=2, label='AET')
    ax1.plot(my_months, my_class.monthly_totals['eq_m'],
             'g-', linewidth=2, label='EET')
    ax1.plot(my_months, my_class.monthly_totals['cwd'],
             'r:', linewidth=2, label='CWD')
    ax1.set_xticks(my_months)
    ax1.set_xlabel('Months', fontsize=18)
    ax1.set_ylabel('ET (mm)', fontsize=18)
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=4, mode="expand", borderaxespad=0., fontsize=18)
    plt.xlim([1, 12])
    plt.show()
    #
    #
    #
    # Plot monthly results
    #
    my_months = range(1, 13)
    my_y_max = my_class.monthly_totals['ep_m'].max()
    fig = plt.figure()
    # [1]
    ax1 = fig.add_subplot(411)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.plot(my_months, my_class.monthly_totals['ep_m'], 'k-', linewidth=2,
             label='$E^p$')
    ax1.plot(my_months, my_class.monthly_totals['ea_m'], 'k--', linewidth=2,
             label='$E^a$')
    ax1.set_ylabel('$E_m$ (mm)', fontsize=14)
    ax1.set_xticks(range(1, 13, 1))
    ax1.set_yticks(range(0, 176, 50))
    plt.xlim([1, 12])
    plt.ylim([0, my_y_max])
    ax1.text(1.25, 150, '(a)', fontsize=12)
    g1 = ax1.legend(loc=1, fontsize=14)
    f1 = g1.get_frame()
    f1.set_linewidth(0)
    # [2]
    ax2 = fig.add_subplot(412)
    plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.plot(my_months, my_class.monthly_totals['cwd'], 'k-', linewidth=2)
    ax2.set_ylabel('$CWD$ (mm)', fontsize=14)
    ax2.set_xticks(range(1, 13, 1))
    ax2.set_yticks(range(0, 176, 50))
    plt.xlim([1, 12])
    plt.ylim([0, my_y_max])
    ax2.text(1.25, 150, '(b)', fontsize=12)
    # [3]
    ax3 = fig.add_subplot(413)
    plt.setp(ax3.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax3.plot(my_months, my_class.monthly_totals['eq_m'], 'k-', linewidth=2,
             label='$E^q$')
    ax3.plot(my_months, my_class.monthly_totals['ea_m'], 'k--', linewidth=2,
             label='$E^a$')
    ax3.set_ylabel('$E_m$ (mm)', fontsize=14)
    ax3.set_xticks(range(1, 13, 1))
    ax3.set_yticks(range(0, 176, 50))
    plt.xlim([1, 12])
    plt.ylim([0, my_y_max])
    ax3.text(1.25, 150, '(c)', fontsize=12)
    g3 = ax3.legend(loc=1, fontsize=14)
    f3 = g3.get_frame()
    f3.set_linewidth(0)
    # [4]
    ax4 = fig.add_subplot(414)
    plt.setp(ax4.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax4.get_xticklabels(), rotation=0, fontsize=12)
    ax4.plot(my_months, my_class.monthly_totals['cpa'], 'k-', linewidth=2)
    ax4.set_ylabel('$\\alpha$', fontsize=14)
    ax4.set_xticks(range(1, 13, 1))
    ax4.set_yticks([0.3*i for i in range(0, 5, 1)])
    plt.xlim([1, 12])
    plt.ylim([0, 1.3])
    ax4.text(1.25, 1.1, '(d)', fontsize=12)
    plt.show()
    #
    #
    #
    # Plot daily ET results
    #
    my_days = range(1, 367)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
    ax1.plot(my_days, my_class.daily_totals['ep_n'],
             'b-', linewidth=2, label='Potential')
    ax1.plot(my_days, my_class.daily_totals['eq_n'],
             'g-', linewidth=2, label='Equilibrium')
    ax1.plot(my_days, my_class.daily_totals['ea_n'],
             'r--', linewidth=2, label='Actual')
    ax1.set_ylabel('Evapotranspiration, mm d$^{-1}$', fontsize=18)
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0., fontsize=18)
    plt.show()
    #
    #
    #
    # Plot daily results
    #
    my_days = range(1, 367)
    my_xtks = range(0, 370, 60)
    fig = plt.figure()
    # [1]
    ax1 = fig.add_subplot(811)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.plot(my_days, my_data['sf'], 'k-', linewidth=2)
    ax1.set_ylabel('$S_f$', fontsize=14)
    ax1.set_xticks(my_xtks)
    ax1.set_yticks([0.1*i for i in range(4, 8, 1)])
    plt.xlim([-15, 380])
    ax1.text(-5, 0.7, '(a)', fontsize=12)
    # [2]
    ax2 = fig.add_subplot(812)
    plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.plot(my_days, (1e-6)*my_class.daily_totals['hn'], 'k-', linewidth=2)
    ax2.set_ylabel('$H_N$ (MJ m$^{-2}$)', fontsize=14)
    ax2.set_xticks(my_xtks)
    ax2.set_yticks(range(6, 19, 3))
    plt.xlim([-15, 380])
    ax2.text(-5, 16, '(b)', fontsize=12)
    # [3]
    ax3 = fig.add_subplot(813)
    plt.setp(ax3.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax3.plot(my_days, my_class.daily_totals['cn'], 'k-', linewidth=2)
    ax3.set_ylabel('$C_n$ (mm)', fontsize=14)
    ax3.set_xticks(my_xtks)
    ax3.set_yticks([0.1*i for i in range(5, 9, 1)])
    plt.xlim([-15, 380])
    ax3.text(-5, 0.75, '(c)', fontsize=12)
    # [4]
    ax4 = fig.add_subplot(814)
    plt.setp(ax4.get_yticklabels(), rotation=0, fontsize=12)
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.plot(my_days, my_data['ppt'], 'k-', linewidth=2)
    ax4.set_ylabel('$P_n$ (mm)', fontsize=14)
    ax4.set_xticks(my_xtks)
    ax4.set_yticks(range(0, 26, 5))
    plt.xlim([-15, 380])
    ax4.text(-5, 20, '(d)', fontsize=12)
    # [5]
    ax5 = fig.add_subplot(815)
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), rotation=0, fontsize=12)
    ax5.plot(my_days, my_class.daily_totals['wn'], 'k-', linewidth=2)
    ax5.set_ylabel('$W_n$ (mm)', fontsize=14)
    ax5.set_xticks(my_xtks)
    ax5.set_yticks(range(30, 151, 30))
    plt.xlim([-15, 380])
    ax5.text(-5, 120, '(e)', fontsize=12)
    # [6]
    ax6 = fig.add_subplot(816)
    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.setp(ax6.get_yticklabels(), rotation=0, fontsize=12)
    ax6.plot(my_days, my_class.daily_totals['ro'], 'k-', linewidth=2)
    ax6.set_ylabel('$RO$ (mm)', fontsize=14)
    ax6.set_xticks(my_xtks)
    ax6.set_yticks(range(0, 21, 5))
    plt.xlim([-15, 380])
    ax6.text(-5, 18, '(f)', fontsize=12)
    # [7]
    ax7 = fig.add_subplot(817)
    plt.setp(ax7.get_xticklabels(), visible=False)
    plt.setp(ax7.get_yticklabels(), rotation=0, fontsize=12)
    ax7.plot(my_days, my_data['tc'], 'k-', linewidth=2)
    ax7.set_ylabel('$T_{air}$ ($^{\circ}C$)', fontsize=14)
    ax7.set_xticks(my_xtks)
    ax7.set_yticks(range(10, 26, 5))
    plt.xlim([-15, 380])
    ax7.text(-5, 23, '(g)', fontsize=12)
    # [8]
    ax8 = fig.add_subplot(818)
    plt.setp(ax8.get_xticklabels(), rotation=0, fontsize=12)
    plt.setp(ax8.get_yticklabels(), rotation=0, fontsize=12)
    ax8.plot(my_days, my_class.daily_totals['ep_n'], 'k-', linewidth=2)
    ax8.plot(my_days, my_class.daily_totals['ea_n'], 'k--', linewidth=2)
    ax8.set_ylabel('$E_n$ (mm)', fontsize=14)
    ax8.set_xticks(my_xtks)
    ax8.set_yticks([1.5*i for i in range(5)])
    plt.xlim([-15, 380])
    ax8.text(-5, 5, '(h)', fontsize=12)
    #
    plt.show()
