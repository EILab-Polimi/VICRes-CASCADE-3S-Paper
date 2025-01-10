# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 10:47:19 2023

initialize_sed_input calculate the daily sed input in each reach, given the annual sed.input.
Annual sed.input is partitioned proportionally to the daily discharge.

This is the python version of initialize_sed_input.m released by Marco Tangi (github: mtangi/DCASCADE_3S)
@author: bruno
"""

import numpy as np

def initialize_sed_input(ReachData, Fi_yield, Q, dates_Q, psi, roundpar, sed_yield, directAd):
    
    # variable initialization
    Qbi_input = [np.zeros((len(ReachData), len(psi))) for _ in range(len(Q))]
    dates_Q = dates_Q[:, :min(len(Q), len(dates_Q[0]))]  # check if dates_Q and Q are of the same length
    Q = Q[:min(len(Q), len(dates_Q[0])), :]

    years_ts = np.unique(dates_Q[0, :])  # years in the simulation

    # sediment yield in t km-2 y-1
    rho_s = 2.650  # sediment density (t/m3)
    sed_yield_year = sed_yield / rho_s * directAd  # sediment yield (m3/yr)

    # calculate daily sed. yield from catchment in Qbi_input
    sed_yield_day = np.zeros((len(ReachData), len(Q), len(psi)))
    j = 0
    i=1
    n=0
    # for each year...
    for i in range(len(years_ts)):
        length_year = np.sum(dates_Q[0, :] == years_ts[i])

        # for each reach...
        for n in range(len(ReachData)):
            daily_Q_perc_n = Q[dates_Q[0, :] == years_ts[i], n] / np.sum(Q[dates_Q[0, :] == years_ts[i], n])  # discharge for each day of the year [Kg/m3]
            sed_yield_year_n_class = sed_yield_year[n] * Fi_yield[n, :]  # annual sediment yield for reach n for each sed.class (m3/yr)

            sed_yield_day_tot_n = np.round(daily_Q_perc_n[:, np.newaxis] * sed_yield_year[n] * Fi_yield[n, :], roundpar)  # daily sediment yield for reach n for each sed.class (m3/day)
            
            # If i just round the sediment input, I risk losing some sediment volumes,
            # espacially for less represented sediment classes.
            # Thus, i need to add the sediment volume lost in the round
            # operation back to the inputs. I do so by adding the difference
            # (diff) between the rounded and non-rounded yield back to the
            # rounded yield matrix.
            diff = np.round(sed_yield_year_n_class - np.sum(sed_yield_day_tot_n,axis=0))
            
            # for each sed.class...
            for d in range(len(diff)):
                # Select the days in the year with the highest yield
                index = np.argsort(-sed_yield_day_tot_n[:, d])  # Sorting in descending order

                # Add or subtract the difference in these days, spread so that the added or sub. value is 1m3 per day
                selected_days = index[:int(abs(diff[d]))]  
                sed_yield_day_tot_n[selected_days, d] += (diff[d] / abs(diff[d])) 
            if np.any(sed_yield_day_tot_n < 0):
                print('Warning: Negative value in reach', n, ':', np.min(sed_yield_day_tot_n))

            sed_yield_day[n, j:j + length_year, :] = sed_yield_day_tot_n

        j = j + length_year
        
    for i in range(len(Qbi_input)):
        Qbi_input[i] = sed_yield_day[:, i, :] 


    return Qbi_input
