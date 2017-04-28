#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# bioindex.py
#
# VERSION: 1.1-dev
# LAST UPDATED: 2017-04-28
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
#
###############################################################################
# IMPORT MODULES:
###############################################################################
import datetime
import os

import numpy

from utilities import find_files
from utilities import print_progress


###############################################################################
# FUNCTIONS
###############################################################################
def add_one_month(dt0):
    """
    Name:     add_one_month
    Inputs:   datetime.date (dt0)
    Outputs:  datetime.date (dt3)
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32)
    dt3 = dt2.replace(day=1)
    return dt3


def calc_monthly_indices(data):
    """
    Computes monthly bioclimatic indices from monthly values
    """
    assert isinstance(data, numpy.ndarray)

    # Instantiate return array
    mo_indices = numpy.array([])

    data_lines = len(data)
    for i in range(data_lines):
        mo_date = data['time'][i]
        mo_aet = data['aet'][i]
        mo_eet = data['eet'][i]
        mo_pet = data['pet'][i]
        mo_pre = data['pre'][i]

        # Calcuate moisture index:
        mi = None
        if mo_pet != 0:
            mi = mo_pre / mo_pet

        # Calculate climatic water deficit:
        cwd = mo_pet - mo_aet

        # Calculate Priestley-Taylor coefficient:
        alpha = None
        if mo_eet != 0:
            alpha = mo_aet / mo_eet

        if i == 0:
            mo_indices = numpy.array(
                (mo_date, mi, cwd, alpha),
                dtype={'names': ('time', 'MI', 'CWD', 'alpha'),
                       'formats': ('O', 'f4', 'f4', 'f4')},
                ndmin=1)
        else:
            temp = numpy.array(
                (mo_date, mi, cwd, alpha),
                dtype={'names': ('time', 'MI', 'CWD', 'alpha'),
                       'formats': ('O', 'f4', 'f4', 'f4')},
                ndmin=1)
            mo_indices = numpy.append(mo_indices, temp, axis=0)

    return mo_indices


def daily_to_monthly(data):
    """
    Calculates monthly sums from daily values
    """
    assert isinstance(data, numpy.ndarray)

    # Instantiate return array
    mo_data = numpy.array([])

    start_date = data['time'][0]
    end_date = data['time'][-1]
    if start_date != end_date:
        cur_date = start_date.replace(day=1)
        mo_counter = 0
        while cur_date < end_date:
            # Increment one month and find daily indexes for given month:
            nxt_date = add_one_month(cur_date)
            mo_idx = numpy.where(
                (data['time'] >= cur_date) & (data['time'] < nxt_date))

            # Compute monthly summations:
            mo_aet = data[mo_idx]['aet'].sum()  # mm/mo
            mo_eet = data[mo_idx]['eet'].sum()  # mm/mo
            mo_pet = data[mo_idx]['pet'].sum()  # mm/mo
            mo_pre = data[mo_idx]['pre'].sum()  # mm/mo

            # Repack into structured array:
            if mo_counter == 0:
                mo_data = numpy.array(
                    (cur_date, mo_aet, mo_eet, mo_pet, mo_pre),
                    dtype={'names': ('time', 'aet', 'eet', 'pet', 'pre'),
                           'formats': ('O', 'f4', 'f4', 'f4', 'f4')},
                    ndmin=1)
            else:
                temp = numpy.array(
                    (cur_date, mo_aet, mo_eet, mo_pet, mo_pre),
                    dtype={'names': ('time', 'aet', 'eet', 'pet', 'pre'),
                           'formats': ('O', 'f4', 'f4', 'f4', 'f4')},
                    ndmin=1)
                mo_data = numpy.append(mo_data, temp, axis=0)

            mo_counter += 1
            cur_date = nxt_date
    return mo_data


def get_date(year, doy):
    """
    Takes a year and day of year and returns a datetime object for that date.
    """
    my_date = datetime.datetime(year, 1, 1)
    my_date += datetime.timedelta(days=doy)
    my_date -= datetime.timedelta(days=1)
    return my_date.date()


def get_output_file(input_file):
    """
    Returns an output file name for a given input file.
    """
    f_name, _ = os.path.splitext(input_file)
    output_file = "%s_bioindex.csv" % (f_name)
    return output_file


def load_file(file_path):
    """
    Name:     load_file
    Inputs:   str, SPLASH output file path (file_path)
    Outputs:  numpy.ndarray, structured data array
    Features: Loads SPLASH output data into a structured array
    Depends:  get_date
    """
    # Instantiate return array
    data = numpy.array([])

    # Header row from the SPLASH class:
    # 0. Year, 1. Day, 2. Lat, 3. Lon, 4. Elv, 5. Precip_mm, 6. Tair_degC,
    # 7. Sf, 8. Cond_mm, 9. EET_mm, 10. PET_mm, 11. AET_mm, 12. SoilMoist_mm,
    # 13. Runoff_mm, 14. NetRadDay_MJ_m2, 15. NetRadNight_MJ_m2
    file_data = numpy.loadtxt(
        fname=file_path,
        dtype={'names': ('year', 'day', 'pre', 'eet', 'pet', 'aet'),
               'formats': ('i4', 'i4', 'f4', 'f4', 'f4', 'f4')},
        delimiter=',',
        usecols=(0, 1, 5, 9, 10, 11),
        skiprows=1
    )

    # Convert year + doy to datetime and repack into structured array:
    data_lines = len(file_data)
    for i in range(data_lines):
        my_year = int(file_data['year'][i])
        my_day = int(file_data['day'][i])
        my_date = get_date(my_year, my_day)
        my_aet = file_data['aet'][i]   # mm/d
        my_eet = file_data['eet'][i]   # mm/d
        my_pet = file_data['pet'][i]   # mm/d
        my_pre = file_data['pre'][i]   # mm/d
        if i == 0:
            data = numpy.array(
                (my_date, my_aet, my_eet, my_pet, my_pre),
                dtype={'names': ('time', 'aet', 'eet', 'pet', 'pre'),
                       'formats': ('O', 'f4', 'f4', 'f4', 'f4')},
                ndmin=1)
        else:
            temp = numpy.array(
                (my_date, my_aet, my_eet, my_pet, my_pre),
                dtype={'names': ('time', 'aet', 'eet', 'pet', 'pre'),
                       'formats': ('O', 'f4', 'f4', 'f4', 'f4')},
                ndmin=1)
            data = numpy.append(data, temp, axis=0)
    return data


def run(input_file):
    """
    Performs all the operations to read, convert and save SPLASH daily output
    to monthly bioclimatic indices.
    """
    if os.path.isfile(input_file):
        out_file = get_output_file(input_file)
        my_data = load_file(input_file)
        monthly_data = daily_to_monthly(my_data)
        monthly_indices = calc_monthly_indices(monthly_data)
        save_to_file(monthly_indices, out_file)


def save_to_file(data, outfile_path):
    """
    Saves monthly bioclimatic indices to file.
    """
    header_str = "Month,MI,CWD,alpha\n"
    data_lines = len(data)
    if os.path.isdir(os.path.dirname(outfile_path)):
        f = open(outfile_path, 'w')
        f.write(header_str)
        for i in range(data_lines):
            mo_time = data['time'][i]
            mo_mi = data['MI'][i]
            mo_cwd = data['CWD'][i]
            mo_alpha = data['alpha'][i]
            out_str = "%s,%f,%f,%f\n" % (mo_time, mo_mi, mo_cwd, mo_alpha)
            f.write(out_str)
        f.close()


###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    home_dir = os.path.expanduser("~")
    data_dir = os.path.join(home_dir, "Data", "splash", "out_for_gepisat")
    my_files = find_files(data_dir, "*.out")
    num_files = len(my_files)

    print_progress(0, num_files)
    for i in range(num_files):
        my_file = my_files[i]
        run(my_file)
        print_progress(i+1, num_files, prefix=os.path.basename(my_file))
