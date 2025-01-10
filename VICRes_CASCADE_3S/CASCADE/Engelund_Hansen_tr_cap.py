# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:28:09 2023

ENGELUND_HANSEN_TR_CAP returns the value of the transport capacity (in Kg/s)
for each sediment class in the reach measured using the Engelund and Hansen equations

references
Engelund, F., and E. Hansen (1967), A Monograph on Sediment Transport in Alluvial Streams, Tekniskforlag, Copenhagen.

This is the python version of Engelung_Hansen_tr_cap.py released by Marco Tangi (github: mtangi/DCASCADE_3S)
@author: bruno
"""

import numpy as np

def Engelund_Hansen_tr_cap(Fi_r_reach, Slope_reach, Wac_reach, v_reach, h_reach, psi):
    
    # Transport capacity from Engelund-Hansen equations using bed material fraction approach (BMF, Molinas and Wu, 2000)
    dmi = 2**(-psi) / 1000  # sediment classes diameter (m)

    rho_s = 2650  # sediment density [kg/m^3]
    rho_w = 1000  # water density [kg/m^3]
    g = 9.81

    # Initialize Qtr_cap
    Qtr_cap = np.zeros_like(Fi_r_reach)

    for d in range(len(dmi)):
        # friction factor
        C = (2 * g * Slope_reach * h_reach) / v_reach**2
        # dimensionless shear stress
        tauEH = (Slope_reach * h_reach) / ((rho_s / rho_w - 1) * dmi[d])
        # dimensionless transport capacity
        qEH = 0.05 / C * tauEH**(5 / 2)
        # dimensionful transport capacity m3/s
        qEH_dim = qEH * np.sqrt((rho_s / rho_w - 1) * g * dmi[d]**3)  # m3/s (formula from the original cascade paper)

        QS_kg = qEH_dim * Wac_reach * rho_s  # kg/s

        QS_EH = QS_kg / rho_s  # m3/s

        Qtr_cap[d] = QS_EH * Fi_r_reach[d]

    return Qtr_cap
