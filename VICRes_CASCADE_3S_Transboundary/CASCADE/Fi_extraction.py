# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 10:33:32 2023

FI_EXTRACTION returns the sediment distribution for the considered sediment
classes based on the D50 and the inverse value of spread s

The function uses the Rosin sediment distribution (Shih and Komar, 1990)

INPUT : 
 D50   = vector containing the D50 value for all reaches of the network
 psi   = vector containing the mean grain size of each sediment class in phi, it must be coherent with the sediment distribution in Fi_Sn (same number of sediment class)  
 s     = inverse measure of the spread of the GSD curve

OUTPUT :
 Fi_Sn = Frequency of sediment for each sediment class, defined for each source node.

This is the python version of FI_extraction.m released by Marco Tangi (github: mtangi/DCASCADE_3S)
@author: bruno
"""

import numpy as np

def Fi_extraction(D50, s, psi):
    # GSD calculation
    # convert D50 in mm
    D50_finer = D50 * 1000

    # sediment classes diameter (mm)
    dmi = 2 ** (-psi)

    # find k parameter
    k = D50_finer / ((-np.log(1 - 50 / 100)) ** (1 / s))

    # find Fi_r
    F = 1 - np.exp(-(np.flip(dmi) / k) ** s)
    F = np.reshape(F, (1, 5))
    F[:, -1] = 1

    Fi_r = np.flip(np.concatenate((F[:, 0].reshape(-1, 1), np.diff(F, axis=1)), axis=1), axis=1)

    return Fi_r
