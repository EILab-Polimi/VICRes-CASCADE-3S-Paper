# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 14:16:24 2023

SED_TRANSFER_FAST takes the matrix of the mobilized layers(V_mob) and the vector of
the sed velocity for each class(v_sed_id) in a reach (n) and returns the 3D matrices containing the
sed volumes in V_mob in their new position at the end of the timestep.

This simple version of the function represents the volume as a point
sediment parcel delivered from the ToN of the reach n. Thus the volume
has a single destination reach and it never get split.

This is the python version of sed_transfer_fast.py released by Marco Tangi (github: mtangi/DCASCADE_3S)
@author: bruno
"""

import numpy as np
from track_sed_position import track_sed_position

def sed_transfer_fast(n, v_sed_day, Lngt, Network,NH):
    # Initialize parameters
    outlet = NH[-1]
    
    DownNode=list(Network['Downstream']['Node'][0])
    DownNode=DownNode[0]
    DownNode=np.squeeze(DownNode)
    
    DownDist=list(Network['Downstream']['Distance'][0])
    DownDist=DownDist[0]
    DownDist=np.squeeze(DownDist)
    
    # Find start and end reach of the sediment volume after the timestep
    # Reach_dest is the ID of the reach where the sediment volume stops after the timestep
    # P_dest is the position from the from_node of the ID reach where the sediment volume stops after the timestep

    if n == outlet:
        reach_dest = np.full((v_sed_day.shape[0],), n)
        p_dest = Lngt[n] + v_sed_day[:, n]
        p_dest = p_dest[:,np.newaxis]
    else:
        # To find p_end, track the position of a sediment parcel starting
        # from the To_Node of the reach n (i.e. the From_node of the downstream reach).
        p_dest, reach_dest = track_sed_position(DownNode[n][0][0], v_sed_day, Lngt, Network, NH,0)
        p_dest = p_dest + Lngt[n]
     
    # Downdist contains the distance from the starting reach to all reaches downstream
    downdist = DownDist[n][0]
    a=Lngt[reach_dest]
    b=downdist[reach_dest]
    a=a[:, np.newaxis]
    b=b[:, np.newaxis]

    # Find position of the sediment volume
    #setout is equal to 1 if the volume in the sed.class left the network
    #via the outlet
    setout = p_dest - a - b > 0
    setout=setout.flatten()

    #in each row, setout is equal to 1 in the reach where the sed. volume
    #of each class is delivered 
    return reach_dest, setout
