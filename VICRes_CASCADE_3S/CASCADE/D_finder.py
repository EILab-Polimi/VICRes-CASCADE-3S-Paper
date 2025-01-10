# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:50:21 2023

D_FINDER find the value of granulometry for the specified D_values, for
the sediment distribution Fi_r, in m

This is the python version of D_Finder.m released by Marco Tangi (github: mtangi/DCASCADE_3S)
@author: bruno
"""

import numpy as np

def D_finder(Fi_r, psi, D_values=None):

    # If D_values not given, find only the D50
    if D_values is None:
        D_values = [50]
    
    if np.isscalar(D_values):
       D_values = [D_values]

    dmi = 2**(-psi) / 1000  # Sediment classes [m]
    Fi_r=Fi_r[:,np.newaxis]

    D_changes = np.zeros((len(D_values), Fi_r.shape[1]))
    Perc_finer = np.zeros((len(dmi), Fi_r.shape[1]))
    Perc_finer[0, :] = 100

    for i in range(1, Perc_finer.shape[0]):
        Perc_finer[i, :] = Perc_finer[i - 1, :] - (Fi_r[i - 1, :] * 100)

    for k in range(Perc_finer.shape[1]):
        for j in range(len(D_values)):
            a = min(np.where(Perc_finer[:, k] > D_values[j])[0][-1], len(psi) - 2)
            D_changes[j, k] = (D_values[j] - Perc_finer[a + 1, k]) / (
                Perc_finer[a, k] - Perc_finer[a + 1, k]
            ) * (-psi[a] + psi[a + 1]) - psi[a + 1]
            D_changes[j, k] = 2 ** D_changes[j, k] / 1000
            D_changes[j, k] = D_changes[j, k] * (D_changes[j, k] > 0) + dmi[-1] * (D_changes[j, k] < 0)

    return D_changes
