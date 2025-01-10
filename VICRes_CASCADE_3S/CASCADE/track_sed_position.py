# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 14:25:32 2023

TRACK_SED_POSITION_TRCAP finds the position of a sediment parcel starting
from reach n after the timestep has passed, defined as the reach ID and
the position from the From_node of the starting reach. 

To satisfy the transport capacity in the ToNode section, the starting
position of the volume is positioned in a location that guarantees that
all of it passes through the ToNode and leaves the reach

This is the python version of track_sed_position.m released by Marco Tangi (github: mtangi/DCASCADE_3S)
@author: bruno
"""

import numpy as np

def track_sed_position(n, v_sed_day, Lngt, Network, NH, start_pos):
    
    # Define starting position
    # The starting position is the position on reach n from which the parcel starts,
    # defined as a fraction of the reach length.
    # If start_pos = 0, the parcel starts from the From_node.
    # If start_pos = 1, the parcel starts from the To_Node.
    n=n-1
    n_class = v_sed_day.shape[0]
    start_pos=0

    # Find path downstream
    # Start_pos (between 0 and 1) defines the position in reach n where the sediment parcel starts,
    # if 1 start from the From_node, if 0 starts from the To_Node.

    timestep = 1

    outlet = NH[-1]

    path2out=list(Network['Downstream']['Path'][0])
    path2out=path2out[0]
    path2out=np.squeeze(path2out)
    path2out=np.squeeze(path2out[n])
    path2out=path2out[outlet]
    path2out=np.squeeze(path2out)
    path2out=path2out-1
    
    downdist_path=list(Network['Downstream']['Distance'][0])
    downdist_path=downdist_path[0]
    downdist_path=np.squeeze(downdist_path)
    downdist_path=np.squeeze(downdist_path[n])
    downdist_path=downdist_path[path2out]
    

    # Path2out contains the path from reach n to the outlet, defined as the IDs of the reach downstream ordered.
    # Downdist_path contains the distance from the reaches in path2out to the reach n From_node
    

    # Find position and destination reach ID
    # Isolate the velocity of the downstream reaches
    v_sed_path = v_sed_day[:, path2out]

    # Change the length of the starting reach according to the starting position, different for each tr.cap
    Lngt_pathout = np.tile(Lngt[path2out], (n_class, 1))
    Lngt_pathout[:, 0] = Lngt_pathout[:, 0] * (1 - start_pos)
    #Lngt_pathout=Lngt_pathout.flatten()

    # Calculate the time (in days) it takes to completely cross a reach
    transit_time = Lngt_pathout / v_sed_path

    # The cumulative of transit_time defines how long it takes to reach each downstream To_Node considering the whole path to the reach
    if transit_time.ndim == 1:
    # If transit_time is 1-dimensional, no need to perform cumulative sum, just return it as cum_tr_time
        cum_tr_time = transit_time
    else:
    # If transit_time is multi-dimensional, perform the cumulative sum along axis 1
        cum_tr_time = np.cumsum(transit_time, axis=1)

   

    # Given cum_tr_time, find the reach where the parcel is after the timestep
    find_point = cum_tr_time - timestep
    
    if find_point.ndim==1:
        find_point = find_point[:,np.newaxis]
        indx_pos = np.nanargmin(find_point, axis=1)
        find_point[find_point<0] = np.nan;
    else:
        find_point[find_point[:, 0] < 0, 0] = 100000
        find_point[find_point<0] = np.nan;
        indx_pos = np.nanargmin(find_point, axis=1)
        find_point[find_point == 100000] = np.nan
        
    
    
    indx_pos = indx_pos.reshape(-1,)
    indx_pos[np.isnan(find_point[:, -1])] = (path2out.size-1)  # If the whole row in find_point is nan, the parcel left the network
    
    # I can find the time remaining for the parcel after it enters the reach indx_pos, needed to find the position of the parcel
    find_time = timestep - cum_tr_time
    if find_time.ndim==1:
        find_time = find_time[:,np.newaxis]
    find_time[find_time[:, 0] < 0, 0] = 100000
    find_time[find_time < 0] = np.nan
   

    # Find the index of the minimum value in the original array (find_time)
    indx_t = np.nanargmin(find_time, axis=1)  # indx_t is the reach before indx_pos

    time_left = find_time[np.arange(find_time.shape[0]), indx_t]
    time_left[time_left == 100000] = np.nan
    time_left[np.isnan(time_left)] = timestep  # If time_left is nan, it means that the parcel remained in the starting reach

    if n == outlet:
        # If the parcel is already at the outlet, use the velocity of the outlet to determine the final position
        # (that will be outside the network)
        v_sed_path=v_sed_path[:,np.newaxis]
        sed_pos = time_left * v_sed_path[np.arange(v_sed_path.shape[0]), 0] + downdist_path
        sed_pos[np.isnan(find_point[:, -1])] = downdist_path + Lngt[outlet] + v_sed_path[np.isnan(find_point[:, -1]), -1] * time_left[np.isnan(find_point[:, -1])]

    else:
        # Find the position of the parcel, given the length of the reaches already passed and the remaining time and velocity in the last reach (indx_pos).
        sed_pos = time_left * v_sed_path[np.arange(v_sed_path.shape[0]), indx_pos] + downdist_path[indx_pos]

        # If the whole row in find_point is nan (the parcel left the network),
        # use the velocity of the outlet to determine the final position (that will be outside the network)
        sed_pos[np.isnan(find_point[:, -1])] = downdist_path[-1] + Lngt[outlet] + v_sed_path[np.isnan(find_point[:, -1]), -1] * time_left[np.isnan(find_point[:, -1])]

    # Outind tells for which sed. size the point fell outside the network (1 - outside, 0 - inside)
    outind = np.isnan(find_point[:, -1])

    sed_pos = sed_pos.reshape(-1, 1)

    # Sed_pos = sed_pos + Lngt[n] * (1 - start_pos)
    if path2out.ndim == 0:
        path2out = np.full((5,), path2out)
        end_reach_ID=path2out[indx_pos]
    else:
        end_reach_ID = path2out[indx_pos] # I find the ID of the destination reach from indx_pos, given the order defined by path2out
      

    return sed_pos, end_reach_ID
