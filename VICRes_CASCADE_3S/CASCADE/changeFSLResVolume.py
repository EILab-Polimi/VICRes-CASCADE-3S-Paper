# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 14:00:24 2023

@author: bruno
"""

import numpy as np

def changeFSLResVolume(Qbi_tr_t, Qbi_dep_t, NH,ResFlag,ResID,numres, phi, vol):
    
    
    """
    changeFSLResVolume changes the reservoir volume at full supply level according to the sediment
    storage in the reservoir; this function is developed to work with the simplified input for the flushing model.

    Parameters:
        DamDatabase_active (list): List of dictionaries containing reservoir data.
        Qbi_tr_t (ndarray): Sediment transport at time t for all reaches (size: M x N).
        Qbi_dep_t (ndarray): Sediment deposition at time t for all reaches (size: M x N).
        flooded_reaches_full (list): List of lists containing indexes of flooded reaches for each reservoir.
        phi (float): Sediment porosity.

    Returns:
        FSL_ResVolumenew (ndarray): Updated reservoir volume at full supply level for all reservoirs (size: 1 x N).
        sed_storage (ndarray): Reservoir storage at time t for all reservoirs (size: 1 x N).
    """
    sed_storage=np.zeros((numres))
    sed_Volume=np.zeros((numres))
    for n in NH:
        if ResFlag[n] == 1 and vol[ResID[n]-1]!=0.001:
            sed_storage_reach = np.sum(Qbi_dep_t[n,:])+np.sum(Qbi_tr_t[n,:])
            sed_storage[ResID[n]-1]=sed_storage_reach#Sediment storage in m3
            sed_Volume[ResID[n]-1]=sed_storage[ResID[n]-1] / (1 - phi) / 1E3 #[1000m3]
    return sed_Volume, sed_storage

    