# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:35:25 2023

VELOCITY_EH returns the velocity of the sediment (in m/s) for each sediment
class for each reach using the Engelund Hansen equations (1967)

OUTPUTS:
 
 v_sed: [cxn] matrix reporting the velocity for each sediment class c for
        each reach n [m/s]
References
Engelund, F., and E. Hansen (1967), A Monograph on Sediment Transport in Alluvial Streams, Tekniskforlag, Copenhagen.

This is the python version of velocity_EH.m released by Marco Tangi (github: mtangi/DCASCADE_3S)
@author: bruno
"""

import numpy as np

def velocity_EH(Slope_reach, Wac_reach, v_reach, h_reach, minvel, phi, psi):
    # active layer definition
    L_a = 0.1 * h_reach  # characteristic vertical length scale for transport.
    
    # sediment velocity with fractional transport capacity
    # by measuring the trasport capacity for the single sed.classes
    # indipendenty, we obtain different values of sed. velocity
    rho_s = 2650  # sediment density [kg/m^3]
    rho_w = 1000  # water density [kg/m^3]
    g = 9.81
    
    dmi = (2 **(-psi)) / 1000
    dmi=np.round(dmi,19)
    C = (2 * g * Slope_reach * h_reach) / (v_reach) ** 2
    Slope_reach = Slope_reach[:, np.newaxis]
    h_reach = h_reach[:, np.newaxis]
    # friction factor
   
    # dimensionless shear stress
    tauEH = (Slope_reach * h_reach) / ((rho_s / rho_w - 1) * dmi)
    tauEH = tauEH.transpose(1, 0)
    
    # dimensionless transport capacity
    qEH = 0.05 / C * (tauEH) ** (5 / 2)
    qEH=qEH.transpose(1,0)
    # dimensionful transport capacity
    qEH_dim = qEH * np.sqrt((rho_s / rho_w - 1) * g * (dmi) ** 3)  # m2/s
    qEH_dim=qEH_dim.transpose(1,0)
    QS_kg = qEH_dim * Wac_reach * rho_s  # kg/s
    QS_EH = QS_kg / rho_s  # m3/s

    # calculate velocity
    v_sed = np.maximum(QS_EH / (Wac_reach * L_a * (1 - phi)), minvel)

    v_sed[:, L_a == 0] = minvel

    return v_sed
