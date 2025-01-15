# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 09:48:53 2023

This is a faster version of DCASCADE_3S. It has been used in this work to perform 
evolutionary multi-objective optimization, significantly speeding up the computational 
time required and allowing the optimization to run for many evaluation functions.
The functions are compiled exploiting the numba library, optimizing the computational time.

@author: bruno
"""
import numba
import sys
import numpy as np
import pandas as pd
import csv
import scipy.io

from numba import jit, prange
from numba import njit
from Fi_extraction import Fi_extraction
from initialize_sed_input import initialize_sed_input
from sockets import Driver, InterfaceSocket, Status

@jit((numba.float64[:,:], numba.float64[:,:],numba.uint16[:],numba.uint8[:],numba.uint16[:],numba.int64,numba.float64,numba.float32[:]),nopython=True, cache=True)
def changeFSLResVolume(Qbi_tr_t, Qbi_dep_t, NH,ResFlag,ResID,numres, phi, vol):
    sed_storage=np.zeros((numres))
    sed_Volume=np.zeros((numres))
    for n in NH:
        if ResFlag[n] == 1 and vol[ResID[n]-1]!=0.001:
            sed_storage_reach = np.sum(Qbi_dep_t[n,:])+np.sum(Qbi_tr_t[n,:])
            sed_storage[ResID[n]-1]=sed_storage_reach#Sediment storage in m3
            sed_Volume[ResID[n]-1]=sed_storage[ResID[n]-1] / (1 - phi) / 1E3 #[1000m3]
    return sed_Volume, sed_storage


@jit((numba.int64,numba.uint16[:],numba.float32[:],numba.int64,numba.float64[:,:]),nopython=True, cache=True)
def exchangedata(numres,ResID,b,t,sed_Volume):
    #print(sed_storage[t, :])
    for i in range(numres+1):
        if i in ResID:
            b[i-1] = sed_Volume[t, i] #1000 m3
    
    return b

@jit((numba.float32[:,:],numba.int64,numba.int64,numba.int64,numba.float64[:],numba.float64[:],numba.float64[:],numba.int64[:],numba.int64[:],numba.int64[:],numba.float64[:,:]),nopython=True, cache=True)
def invertQ(Q,n_reaches,numres,t,n_Man, Wac, Slope, startflush,startsluicing,reserlength,ResVolume):
    Q1=Q[0,0:n_reaches]
    Q2=np.copy(Q1)
    Qfinal=Q2[n_reaches-1]
    for l in range (75,n_reaches):
        Q1[l]=Q2[l-1]
    Q1[74]=Qfinal
    vol=Q[0,n_reaches:n_reaches+numres]
    #vol[np.isclose(vol, -1, atol=0.001)] = -1
    #vol[np.isclose(vol, -2, atol=0.001)] = -2
    #vol[np.isclose(vol, 0.001, atol=0.0001)] = 0.001
    #vol[np.isclose(vol, 10, atol=0.001)] = 10
    #print(vol)
    #vol = vol.astype(np.float64)
    reserflushed = np.where(np.abs(vol + 1) < 0.0001)[0] #index of the reservoir
    resersluiced = np.where(vol < -3)[0] #index of the reservoir
    startflush[reserflushed] = 1 #Determine the reservoir that is starting sediment management
    startsluicing[resersluiced] = 1 
    #DEBUG
    reserlength[resersluiced] = -vol[resersluiced]   
    finishedflush = np.where(np.abs(vol + 2) < 0.0001)[0] #index of reservoir that finishes flushing
    startflush[finishedflush] = 0 #flushing is finished for that reservoir
    startsluicing[finishedflush] = 0 #Sluicing is finished
    #From VIC the volume is in 1000 m3
    ResVolume[t,:]=vol*1000 #m3
    #print('VOLUME UGUALE A',ResVolume[t,:])
    # variables initialization
    h = (Q1[:] * n_Man / (Wac * np.sqrt(Slope))) ** (3 / 5)
    v = 1 / n_Man * h ** (2 / 3) * np.sqrt(Slope)

    return Q1, vol, startflush, startsluicing, reserlength, ResVolume,h,v

@jit(numba.types.Tuple((numba.float64[:,:,:],numba.float64[:,:,:]))(numba.float64[:,:,:],numba.float64[:,:,:],numba.float64[:,:],numba.float64[:],numba.float64[:],numba.float64[:],numba.int64,numba.int64,numba.float64[:],numba.float64,numba.float64[:],numba.uint8[:],numba.uint16[:],numba.float32[:],numba.int64[:],numba.int64[:],numba.float32[:],numba.float64[:,:], numba.int64, numba.float64[:],numba.boolean[:],numba.boolean[:],numba.boolean[:],numba.int64[:],numba.int64[:],numba.float64[:,:,:],numba.float64[:,:,:]),nopython=True,parallel=True, cache=True)
def parallelloop(V_dep_numba,Qbi_incoming_numba,Qbi_input_t,EnergySlope,h,v,rho_s,rho_w,dmi,g,Wac,ResFlag,ResID,vol,startflush,PSI,Q1, ResVolume,t, psi_3S,fine_sand_indices,coarse_sand_indices,medium_indices,startsluicing,reserlength,Qbi_dep,Qbi_mob):
    
    for n in prange(463):
        # extract the deposit layer of the reach from the relative cell in the previous timestep

        V_dep_old =V_dep_numba[:, n, :]

        # 1) extract the deposit layer from the storage matrix and load the incoming cascades
        Qbi_incoming = Qbi_incoming_numba[:, n, :]
        Qbi_incoming[n, :] = Qbi_input_t[n, :]


        # 2) find cascades to be included into the active layer according to the limit V_lim_tot,
        # and use the cumulative GSD to compute tr_cap
        Fi_r_reach = np.sum((V_dep_old + Qbi_incoming), axis=0) / np.sum((V_dep_old + Qbi_incoming))
        Fi_r_reach[np.isnan(Fi_r_reach)] = 0

         # Transport capacity from Engelund-Hansen equations using bed material fraction approach (BMF, Molinas and Wu, 2000)


        # Initialize Qtr_cap


        #for d in range(len(dmi)):
            # friction factor
        #    C = (2 * g * EnergySlope[n] * h[n]) / v[n]**2
        #    if C < 0.0001:
        #        C=0.0001
        C = (2 * g * EnergySlope[n] * h[n]) / v[n]**2
        if C < 0.0001:
            C=0.0001
            # dimensionless shear stress
        #    tauEH = (EnergySlope[n] * h[n]) / ((rho_s / rho_w - 1) * dmi[d])
        tauEH = (EnergySlope[n] * h[n]) / ((rho_s / rho_w - 1) * dmi)    
            # dimensionless transport capacity
        #    qEH = 0.05 / C * tauEH**(5 / 2)
        qEH = 0.05 / C * tauEH**(5 / 2)    
            # dimensionful transport capacity m3/s
        #    qEH_dim = qEH * np.sqrt((rho_s / rho_w - 1) * g * dmi[d]**3)  # m3/s (formula from the original cascade paper)
        qEH_dim = qEH * np.sqrt((rho_s / rho_w - 1) * g * dmi**3)
        QS_kg = qEH_dim * Wac[n] * rho_s  # kg/s

        QS_EH = QS_kg / rho_s  # m3/s

        #    Qtr_cap[d] = QS_EH * Fi_r_reach[d]
        Qtr_cap = QS_EH * Fi_r_reach

        tr_cap = Qtr_cap * 24 * 60 * 60
        tr_cap[np.isnan(tr_cap)] = 0
        # If flushing is active, we compute the transport capacity using the IRTCES formula 
        if ResFlag[n] == 1 and abs(vol[ResID[n]-1] + 1) > 0.0001 and vol[ResID[n]-1]>0.0001 and startflush[ResID[n]-1]==1:
            tr_cap = PSI*vol[ResID[n]-1] #t/s
            tr_cap = tr_cap * 1000 / rho_s  #m3/s
            tr_cap = tr_cap * Fi_r_reach* 24 * 60 * 60 #m3/d
        else: 
            tr_cap = tr_cap


        if ResFlag[n] == 1 and abs(vol[ResID[n]-1] - 0.001) > 0.0001 and startflush[ResID[n]-1]==0:
            
            # V_mob contains the volume mobilized in the reach, that is about to be transferred downstream, without the presence of the reservoir
            #Q in m3/s, volume in m3, residence time in d^-1
            Qday=Q1[n]*3600*24
            if Qday == 0:
                Qday = 0.0001
            ResidenceTime = ResVolume[t,ResID[n]-1]/Qday
            if ResidenceTime == 0:
                ResidenceTime = 0.0001
            #print(vol[ResID[n]-1])
            #print(ResID[n]-1)
            #print(ResVolume[t,ResID[n]-1])
            
            Te=np.zeros(len(psi_3S))
            
            # Fine Sand or finer: Gill equation for fine grained sediments (< 0.125 mm)
            Te[fine_sand_indices] = (ResidenceTime**(3))/(1.0265*ResidenceTime**(3)+0.02621*ResidenceTime**(2)-0.133*10**(-3)*ResidenceTime+0.1*10**(-5))
            # Coarse Sand or Gravel: Jothiprakash and Grag (2008), coarse sediment
            
            Te[coarse_sand_indices] = (8000 - 36 * (ResidenceTime**(-0.78))) / (78.85 + (ResidenceTime**(-0.78)))
            Te[coarse_sand_indices] /= 100

            # Medium: Jothiprakash and Grag (2008), medium sediment
            
            Te[medium_indices] = ResidenceTime / (0.00013 + 0.01 * ResidenceTime + 0.0000166 * (ResidenceTime**0.5))
            Te[medium_indices] /= 100

            # Clamp values to [0, 1]
            Te[Te < 0] = 0
            Te[Te > 1] = 1

            # Handle NaN values
            Te[np.isnan(Te)] = 0
            # If we have sluicing, we use the Churchill curve to estimate TE
            if startsluicing[ResID[n]-1]==1:
                si=((ResVolume[t,ResID[n]-1]/Q1[n])**2)/reserlength[ResID[n]-1]
                Te[:] = (112 - 800 * (0.3048 * si) ** (-0.2)) / 100
                Te[Te < 0] = 0
                Te[Te > 1] = 1
                Te[np.isnan(Te)] = 0
            V_mob = Qbi_incoming * (1-Te)
            
            V_dep = np.maximum(0, (V_dep_old + Qbi_incoming - V_mob))
            
            Qbi_dep[:, n, :] = V_dep
            # Qbi_mob contains the volume mobilized in the reach, that is about to be transferred downstream
            Qbi_mob[:, n, :] = V_mob
        # 3) Deposit the cascades in the active layer until the volume mobilized for each class is equal to the tr_cap                                                
        else:
            prc = np.minimum(tr_cap / np.sum((V_dep_old + Qbi_incoming), axis=0), 1)
            prc[np.isnan(prc)] = 0
            V_mob = (V_dep_old + Qbi_incoming) * prc
            
            V_dep = np.maximum(0, (V_dep_old + Qbi_incoming - V_mob))
            
            Qbi_dep[:, n, :] = V_dep
            # Qbi_mob contains the volume mobilized in the reach, that is about to be transferred downstream
            Qbi_mob[:, n, :] = V_mob
    return Qbi_dep, Qbi_mob
    
@jit(numba.types.Tuple((numba.float64[:,:,:],numba.float64[:,:,:],numba.float64[:,:,:],numba.float64[:,:]))(numba.uint16[:],numba.float64[:,:,:],numba.float64[:,:,:],numba.float64[:,:],numba.float64[:],numba.float64[:],numba.float64[:],numba.int64,numba.float64[:],numba.uint8[:],numba.uint16[:],numba.float32[:],numba.int64[:],numba.int64[:],numba.float32[:],numba.float64[:,:], numba.int64, numba.float64[:],numba.int64[:],numba.int64[:],numba.float64,numba.float64,numba.int64[:,:], numba.float64[:],numba.int64[:],numba.float64[:,:],numba.uint16, numba.float64[:,:,:], numba.float64[:,:,:], numba.float64[:,:,:], numba.float64[:,:]),nopython=True, cache=True)
def looping(NH,V_dep_numba,Qbi_incoming_numba,Qbi_input_t,EnergySlope,h,v,rho_s,Wac,ResFlag,ResID,vol,startflush,PSI,Q1,ResVolume,t,psi_3S,startsluicing,reserlength,phi,minvel,all_path2oute, Lngt, DownNode_float,DownDist_float,outlet,Qbi_tr,Qbi_dep,Qbi_mob,Q_out):
    zero = np.zeros((1,), dtype=np.bool_)
    fine_sand_indices = 2**(-psi_3S) < 0.125
    coarse_sand_indices = 2**(-psi_3S) > 1
    medium_indices = ~fine_sand_indices & ~coarse_sand_indices                                                                                                                         
    v[v < 0.0001] = 0.0001
    dmi = 2**(-psi_3S) / 1000  # sediment classes diameter (m)
    start_pos=0
    
    timestep = 1
    rho_s = 2650  # sediment density [kg/m^3]
    rho_w = 1000  # water density [kg/m^3]
    g = 9.81
    
    Q_out.fill(0)  
    Qbi_dep.fill(0)  
    Qbi_mob.fill(0)  
    Qbi_tr.fill(0)
    Qbi_dep, Qbi_mob=parallelloop(V_dep_numba,Qbi_incoming_numba,Qbi_input_t,EnergySlope,h,v,rho_s,rho_w,dmi,g,Wac,ResFlag,ResID,vol,startflush,PSI,Q1, ResVolume,t, psi_3S,fine_sand_indices,coarse_sand_indices,medium_indices,startsluicing,reserlength,Qbi_dep,Qbi_mob)

    L_a = 0.1 * h  # characteristic vertical length scale for transport.
    
    C = (2 * g * EnergySlope * h) / (v) ** 2
    mask = C < 0.0001

    C[mask] = 0.00001
    EnergySlope = EnergySlope[:, np.newaxis]
    h = h[:, np.newaxis]
    # friction factor
   
    # dimensionless shear stress
    tauEH = (EnergySlope * h) / ((rho_s / rho_w - 1) * dmi)
    tauEH = tauEH.transpose(1, 0)
    
    # dimensionless transport capacity
    qEH = 0.05 / C * (tauEH) ** (5 / 2)
    qEH=qEH.transpose(1,0)
    # dimensionful transport capacity
    qEH_dim = qEH * np.sqrt((rho_s / rho_w - 1) * g * (dmi) ** 3)  # m2/s
    qEH_dim=qEH_dim.transpose(1,0)
    QS_kg = qEH_dim * Wac * rho_s  # kg/s
    QS_EH = QS_kg / rho_s  # m3/s

    # calculate velocity
    v_sed = np.maximum(QS_EH / (Wac * L_a * (1 - phi)), minvel)

    v_sed[:, L_a == 0] = minvel
    for n in range(463):
        if ResFlag[n]==1 and abs(vol[ResID[n]-1] - 0.001) > 0.0001:
            if startflush[ResID[n]-1]==1 and abs(vol[ResID[n]-1] + 1) > 0.0001: #flushing is active and drawdown is terminated
                v_sed[:,n]=v_sed[:,n]
            else:
                v_sed[:,n]=(1-0.99)*v_sed[:,n]
        else:
            v_sed[:,n]=v_sed[:,n]
    v_sed_day=v_sed*60*60*24
    n_class = v_sed_day.shape[0]
    reach_dest=np.empty(v_sed_day.shape[0], dtype=np.int64)
    
    for n in range(463):
        V_mob = Qbi_mob[:, n, :]
        if np.sum(V_mob) > 0:
            maske = all_path2oute[n]
            mask = maske != 0
            path2oute = all_path2oute[n,mask]
            if n == outlet:
                
                reach_dest[:]=n
                
                p_dest = Lngt[n] + v_sed_day[:, n]
                
            else:
                # To find p_end, track the position of a sediment parcel starting
                # from the To_Node of the reach n (i.e. the From_node of the downstream reach).
                
                e = DownNode_float[n]
                e=e-1
                
                
                #path2oute=np.squeeze(path2out[e])
                
                #path2oute=np.squeeze(path2oute[outlet])
                
                path2oute=path2oute-1
                
                
                downdist_pathe=DownDist_float[e,path2oute]
                downdist_pathe=downdist_pathe+zero
                

                # Path2out contains the path from reach n to the outlet, defined as the IDs of the reach downstream ordered.
                # Downdist_path contains the distance from the reaches in path2out to the reach n From_node
                

                # Find position and destination reach ID
                # Isolate the velocity of the downstream reaches
                v_sed_path = v_sed_day[:, path2oute]
                
                path2oute= path2oute+zero
                # Change the length of the starting reach according to the starting position, different for each tr.cap
                #Lngt_pathout = np.tile(Lngt[path2oute], (n_class, 1))
                #if e == outlet:
                #    Lngt_pathouti = np.empty((n_class,))
                    #for i in range(n_class):
                #    Lngt_pathouti[:] = Lngt[path2oute]
                #    transit_time = Lngt_pathouti / v_sed_path
                #else:
                Lngt_pathout = np.atleast_2d(np.empty((n_class,)+ np.shape(Lngt[path2oute])))
                if Lngt_pathout.shape[0]==1:
                    Lngt_pathout=Lngt_pathout.T
                Lngt_pathout[:] = Lngt[path2oute]
                Lngt_pathout[:, 0] = Lngt_pathout[:, 0] * (1 - start_pos)
                    #Lngt_pathout = np.vstack([Lngt[path2oute]] * n_class)
                    #if len(Lngt_pathout.shape) == 1:
                    #    Lngt_pathout = np.expand_dims(Lngt_pathout, axis=1)
                    #Lngt_pathout[:, 0] = Lngt_pathout[:, 0]  #* (1 - start_pos)
                    #Lngt_pathout=Lngt_pathout.flatten()
                #print('lngtpathout',Lngt_pathout)
                    # Calculate the time (in days) it takes to completely cross a reach
                transit_time = np.atleast_2d(Lngt_pathout / v_sed_path)
                v_sed_path = v_sed_day[:, path2oute]
                #print('transittime',transit_time)
                # The cumulative of transit_time defines how long it takes to reach each downstream To_Node considering the whole path to the reach
                if transit_time.shape[0] == 1:
                # If transit_time is 1-dimensional, no need to perform cumulative sum, just return it as cum_tr_time
                    cum_tr_time = transit_time.T
                else:
                # If transit_time is multi-dimensional, perform the cumulative sum along axis 1
                    #cum_tr_time = np.cumsum(transit_time, axis=1)
                    cum_tr_time = np.empty_like(transit_time)
                    for i in range(transit_time.shape[1]):
                        cum_tr_time[:, i] = np.sum(transit_time[:, :i+1], axis=1)

               

                # Given cum_tr_time, find the reach where the parcel is after the timestep
                find_point = cum_tr_time - timestep


                find_point[find_point[:, 0] < 0, 0] = 100000
                find_point = np.where(find_point < 0, 100000000, find_point)
                #find_point[find_point<0] = np.nan;
                indx_pos = np.argmin(find_point, axis=1)
                #find_point[find_point == 100000] = np.nan
                find_point = np.where(find_point == 100000, np.nan, find_point)
                find_point = np.where(find_point == 100000000, np.nan, find_point)
                
                
                indx_pos = indx_pos.reshape(-1,)
                indx_pos[np.isnan(find_point[:, -1])] = (path2oute.size-1)  # If the whole row in find_point is nan, the parcel left the network
                # I can find the time remaining for the parcel after it enters the reach indx_pos, needed to find the position of the parcel
                find_time = timestep - cum_tr_time
                

                find_time[find_time[:, 0] < 0, 0] = 100000
                find_time= np.where(find_time < 0, 100000000, find_time)
                #find_time[find_time < 0] = np.nan
                   
                indx_t = np.argmin(find_time, axis=1)  # indx_t is the reach before indx_pos
                #find_time= np.where(find_time == 100000, np.nan, find_time)
                time_left = np.zeros_like(indx_t, dtype=find_time.dtype)
                for i in range(find_time.shape[0]):
                    time_left[i] = find_time[i, indx_t[i]]
                #time_left = find_time[np.arange(find_time.shape[0]), indx_t]
                time_left[time_left == 100000] = np.nan
                time_left[time_left == 100000000] = np.nan
                time_left[np.isnan(time_left)] = timestep  # If time_left is nan, it means that the parcel remained in the starting reach
                
                sed_pos = np.zeros((indx_pos.shape[0]), dtype=find_time.dtype)
                if e == outlet:
                #     # If the parcel is already at the outlet, use the velocity of the outlet to determine the final position
                #     # (that will be outside the network)
                    v_sed_path=np.atleast_2d(v_sed_path)
                    sed_pos = time_left * v_sed_path[np.arange(v_sed_path.shape[0]), 0] + downdist_pathe
                    sed_pos[np.isnan(find_point[:, -1])] = downdist_pathe + Lngt[outlet] + v_sed_path[np.isnan(find_point[:, -1]), -1] * time_left[np.isnan(find_point[:, -1])]

                else:
                    
                    # Find the position of the parcel, given the length of the reaches already passed and the remaining time and velocity in the last reach (indx_pos).
                    for i in range(v_sed_path.shape[0]):
                        sed_pos[i] = time_left[i] * v_sed_path[i, indx_pos[i]] + downdist_pathe[indx_pos[i]]
                    #sed_pos = time_left * v_sed_path[np.arange(v_sed_path.shape[0]), indx_pos] + downdist_pathe[indx_pos]
                    # If the whole row in find_point is nan (the parcel left the network),
                    # use the velocity of the outlet to determine the final position (that will be outside the network)
                    sed_pos[np.isnan(find_point[:, -1])] = downdist_pathe[-1] + Lngt[outlet] + v_sed_path[np.isnan(find_point[:, -1]), -1] * time_left[np.isnan(find_point[:, -1])]

                # Outind tells for which sed. size the point fell outside the network (1 - outside, 0 - inside)
                #outind = np.isnan(find_point[:, -1])
                p_dest = sed_pos
                if path2oute.shape[0] == 1:
                    prova=np.empty((n_class,),dtype=np.int64)
                    prova[:]=path2oute
                    end_reach_ID=prova[indx_pos]
                else:
                    end_reach_ID = path2oute[indx_pos] # I find the ID of the destination reach from indx_pos, given the order defined by path2out
                
                

                
                reach_dest = end_reach_ID
                

                p_dest = p_dest + Lngt[n]
             
            # Downdist contains the distance from the starting reach to all reaches downstream
            
            downdisti = DownDist_float[n]
            a=Lngt[reach_dest]
            b=downdisti[reach_dest]
            
            # Find position of the sediment volume
            setout = p_dest - a - b > 0
            setout=setout.flatten()


            
            for c in range(len(psi_3S)):   
                if setout[c] == 0:
                    Qbi_tr[:, reach_dest[c], c] += V_mob[:, c]
                else:
                    Q_out[:, c] =Q_out[:, c] + V_mob[:, c]
    
    return Qbi_dep, Qbi_mob, Qbi_tr, Q_out    

def savedata(Q_out,timescale,port):
    outcum_tot = np.array([np.sum(Q_out[t,:,:], axis=0).sum() for t in range(timescale)])
    a = np.zeros((timescale)// 365)
    tot_sed_year = np.zeros((timescale) // 365)
    

    for i in range((timescale) // 365):
        a[i] = np.sum(outcum_tot[365 * (i):365 * (i+1)])
        tot_sed_year[i] = a[i] * 2.6 / 1000000
        
    
    
    file_path = f'../VICRes/RoutingSetupCores/RoutingSetup3S{port}/Results/totsedyear.csv'
    np.savetxt(file_path, tot_sed_year, delimiter=',')

def main():
    port = int(sys.argv[1])
    dtype = np.float32
    #Socket!
    server = InterfaceSocket(port=port)
    server.open()
    client, address = server.server.accept()
    #end socket

    # client.settimeout(server.timeout)
    driver = Driver(client)
    print(" @SOCKET:   Client asked for connection from " + str(address) + ". Now hand-shaking.")

    stat = driver.get_status()
    if stat == Status.Up | Status.NeedsInit:
        driver.initialise()
        print('Driver initialised.')

    #def DCASCADE_3S_sedyield(sed_input_param):
    # Load the .mat file
    data = scipy.io.loadmat('network_data_3Snew.mat')

    # Extract the variables from the loaded data
    network = data['Network'][0, 0]
    reach_data = data['ReachData']
    psi_3S = data['psi_3S'][0]
    dmi_finer = 2.0**(-psi_3S)
    PSI = np.where(dmi_finer < 0.016, 1600, np.where((dmi_finer >= 0.016) & (dmi_finer < 0.1), 650, 300)) #constant for flushing transport capacity
    rho_s = 2650 #sediment density [kg/m3]
    n_reaches = len(reach_data)
    numres=109
    sh = (1,n_reaches+numres)
    DownNode=list(network['Downstream']['Node'][0])
    DownNode=DownNode[0]
    DownNode=np.squeeze(DownNode)
    DownNode_float = np.zeros_like(DownNode, dtype=int)
    for i in range(len(DownNode)):
        if i != 74:
            DownNode_float[i] = DownNode[i][0][0]
    
    DownDist=list(network['Downstream']['Distance'][0])
    DownDist=DownDist[0]
    DownDist=np.squeeze(DownDist)
    num_reaches = len(DownDist)
    max_num_elements = max(len(item[0]) for item in DownDist)
    DownDist_float = np.zeros((num_reaches, max_num_elements))

    for i in range(num_reaches):
        DownDist_float[i, :len(DownDist[i][0])] = DownDist[i][0]

    path2out=list(network['Downstream']['Path'][0])
    path2out=path2out[0]
    path2out=np.squeeze(path2out)

    downdist_path=list(network['Downstream']['Distance'][0])
    downdist_path=downdist_path[0]
    downdist_path=np.squeeze(downdist_path)
    # Convert the MATLAB struct arrays to Python lists
    NH = network['NH']
    NH = NH.flatten()
    NH=NH-1
    Wa = [reach_data[i]['Wac'] for i in range(len(reach_data))]
    Wa = np.array([d[0] for d in Wa])
    Wa = Wa.flatten()
    Wac=np.copy(Wa)
    ResFlag = [reach_data[i]['Res'] for i in range(len(reach_data))]
    ResFlag = np.array([d[0] for d in ResFlag])
    ResFlag = ResFlag.flatten()

    ResDownReach = [reach_data[i]['ResDownReach'] for i in range(len(reach_data))]
    ResDownReach = np.array([d[0] for d in ResDownReach])
    ResDownReach = ResDownReach.flatten()
    ResDownReach = ResDownReach-1

    ResID = [reach_data[i]['ResID'] for i in range(len(reach_data))]
    ResID = np.array([d[0] for d in ResID])
    ResID = ResID.flatten()

    ResName = [reach_data[i]['ResName'] for i in range(len(reach_data))]
    ResName = np.array([d[0] for d in ResName])
    ResName = ResName.flatten()

    Slope = [reach_data[i]['Slope'] for i in range(len(reach_data))]
    Slope = np.array([d[0] for d in Slope])
    Slope = Slope.flatten()
    Lngt = [reach_data[i]['Length'] for i in range(len(reach_data))]
    Lngt = np.array([d[0] for d in Lngt])
    Lngt = Lngt.flatten()
    n_Man = [reach_data[i]['n'] for i in range(len(reach_data))]
    n_Man = np.array([d[0] for d in n_Man])
    n_Man = n_Man.flatten()
    fromN = [reach_data[i]['FromN'] for i in range(len(reach_data))]
    fromN = np.array([d[0] for d in fromN])
    fromN = fromN.flatten()
    toN = [reach_data[i]['ToN'] for i in range(len(reach_data))]
    toN = np.array([d[0] for d in toN])
    toN = toN.flatten()
    el_FN = [reach_data[i]['el_FN'] for i in range(len(reach_data))]
    el_FN = np.array([d[0] for d in el_FN])
    el_FN = el_FN.flatten()
    el_TN = [reach_data[i]['el_TN'] for i in range(len(reach_data))]
    el_TN = np.array([d[0] for d in el_TN])
    el_TN = el_TN.flatten()
    sed_yield = [reach_data[i]['sed_yield'] for i in range(len(reach_data))]
    sed_yield = np.array([d[0] for d in sed_yield])
    sed_yield = sed_yield.flatten()
    directAd = [reach_data[i]['directAd'] for i in range(len(reach_data))]
    directAd = np.array([d[0] for d in directAd])
    directAd = directAd.flatten()
    II = network['II']

    #Distance=[Upstream[i]['Distance'] for i in range(462)]

    # Load Q data
    Qall =scipy.io.loadmat('Q_Historical.mat')
    dates_Q =scipy.io.loadmat('dates.mat')
    Qall = Qall['Q']
    dates_Q=dates_Q['dates_Q']

    #dates_Q = q_data['data_structure']['Q_dates']

    # time initialization
    #timescale = len(Qall)
    #timescale=  365+366+365
    timescale=  5479
    # load old variables
    sed_yield = [55]
    sed_D50 = [0.25]


    scen=[1]

    # Initialize output elements
    sed_input_param_comb = np.array(np.meshgrid(sed_yield, sed_D50,scen)).T.reshape(-1, 3).T
    sed_input_param=sed_input_param_comb[:, 0]

    outlet = np.where(fromN == toN)[0]
    outlet=NH[-1]


    # define input sediment load
    # define magnitude of the input sed yield, in t/km2/y
    tot_yield = sed_input_param[0] * 1e6  # total yield for the 3S, in t/y
    yield_km_y = tot_yield / np.sum(sed_yield / np.max(sed_yield) * directAd)
    sed_yield = sed_yield / np.max(sed_yield) * yield_km_y

    # define GSD of the input sed yield
    Fi_yield = np.tile(Fi_extraction(sed_input_param[1] / 1000, 0.8, psi_3S), (len(NH), 1))


    phi = 0.4
    minvel = np.min(Slope)
    
    # calculate input sed. yield


    Qbi_input = initialize_sed_input(reach_data, Fi_yield, Qall, dates_Q, psi_3S, 10, sed_yield, directAd)
    last_matrix = Qbi_input[-1]

   
    Qbi_input.append(last_matrix)
    
    ResVolume = np.zeros((timescale + 1, numres))
    sed_Volume = np.zeros((timescale + 1, numres))
    sed_storage =  np.zeros((timescale + 1, numres))




    # Routing scheme
    # initialize new variables
    rows, col = II.shape
    Qbi_tr = np.zeros((2, rows, col, len(psi_3S)))
    Qbi_mob = np.zeros((1, rows, col, len(psi_3S)))
    Qbi_dep = np.zeros((2, rows, col, len(psi_3S)))
    Q_out = np.zeros((timescale + 1, len(II), len(psi_3S)))
    #Qbi_tr = [np.zeros((rows,col, len(psi_3S))) for _ in range(2)]
    #Qbi_mob = [np.zeros((rows,col, len(psi_3S))) for _ in range(1)]
    #Qbi_dep = [np.zeros((rows,col, len(psi_3S)))  for _ in range(2)]
    #Q_out = [np.zeros((len(II), len(psi_3S))) for _ in range(timescale+1)]
    #Q_original = np.copy(Q)

    EnergySlope = np.copy(Slope) 
    # Routing scheme
    #setout=[np.zeros((len(II), 5)) for _ in range(timescale+1)]
    max_index = 463
    max_size=463
    
    all_path2oute = np.zeros((max_index, max_size), dtype=np.int64)  # Dove max_size è la dimensione massima possibile per ogni array path2oute

    for n in range(463):
        e = DownNode_float[n]
        e = e - 1
        path2oute = np.squeeze(path2out[e])
        path2oute = np.squeeze(path2oute[outlet]).astype(np.int64)
        if np.isscalar(path2oute):  
            path2oute = np.array([path2oute])  

        size = min(path2oute.size, max_size)  
        all_path2oute[n, :size] = path2oute
    
    Q = np.zeros(sh, dtype)
    startflush = np.zeros(numres,dtype=int)
    startsluicing = np.zeros(numres,dtype=int)
    reserlength = np.zeros(numres,dtype=int)
    
    
    for t in range(0,timescale):
            stat = driver.get_status()
            
            if stat == Status.Up | Status.HasData:

                Q = driver.get_data(Q)
                
                
                
                
                Q1, vol, startflush, startsluicing, reserlength, ResVolume, h, v = invertQ(Q,n_reaches,numres,t,n_Man, Wac, Slope,startflush,startsluicing,reserlength,ResVolume)
                
                V_dep_numba =Qbi_dep[0,:,:,:]

                # 1) extract the deposit layer from the storage matrix and load the incoming cascades
                Qbi_incoming_numba = Qbi_tr[0,:,:,:]
                Qbi_input_t = Qbi_input[t]      
                Qbi_dep[1,:,:,:], Qbi_mob[0,:,:,:], Qbi_tr[1,:,:,:],  Q_out[t,:,:]= looping(NH,V_dep_numba,Qbi_incoming_numba,Qbi_input_t,EnergySlope,h,v,rho_s,Wac,ResFlag,ResID,vol,startflush,PSI,Q1,ResVolume,t,psi_3S,startsluicing,reserlength,phi,minvel,all_path2oute, Lngt, DownNode_float,DownDist_float,outlet, Qbi_tr[1,:,:,:], Qbi_dep[1,:,:,:], Qbi_mob[0,:,:,:], Q_out[t,:,:])
                        
                    
                # 5) Move the mobilized volumes to the destination reaches according to the sediment velocity
                #v_sed = velocity_EH(EnergySlope, Wac, v, h, minvel, phi, psi_3S)
                
                #If we have a reservoir, the velocity of the reach should be almost 0
                #v_sed = changevelocity(NH,ResFlag,vol,ResID,startflush,v_sed)
                
                Qbi_dep[0,:,:,:] = Qbi_dep[1,:,:,:]
                
                Qbi_tr[0,:,:,:] = Qbi_tr[1,:,:,:]
                Qbi_tr_t = np.squeeze(np.sum(Qbi_tr[1,:,:,:], axis=0))
                Qbi_dep_t = np.squeeze(np.sum(Qbi_dep[1,:,:,:], axis=0))

                sed_Volume[t, :], sed_storage[t, :] = changeFSLResVolume(Qbi_tr_t, Qbi_dep_t, NH,ResFlag,ResID,numres, phi, vol)
                sed_Volume[np.isnan(sed_Volume)] = 0
                sed_storage[np.isnan(sed_storage)] = 0
                
                
                        
            stat = driver.get_status()   
            if stat == Status.Up | Status.Ready:
                
                b = np.zeros(numres, dtype=vol.dtype)
                b = exchangedata(numres,ResID,b,t,sed_Volume)
                
                driver.send_data(b)
    
                
                
    # output processing
    
    savedata(Q_out,timescale,port)



    
    
if __name__ == "__main__":
    main()


