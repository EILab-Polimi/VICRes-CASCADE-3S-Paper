# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 09:48:53 2023

DCASCADE_3S runs an efficient version of D-CASCADE that does not trace
deposit layering to save computational time while still tracing sediment provenance. 

This version builds on DCASCADE_3S, released by Marco Tangi (github: mtangi/DCASCADE_3S), 
it has been translated in python and adds the coupling of the model with VICRes through a Socket interface.
It receives daily discharges for each river reach at each timestep from VICRes and returns 
the amount of sediment deposited in each reservoir.

The 16 acceptable scenarios (annual sediment yield of the basin and GSD of the sediment yield) 
identified in the first step of the work are saved in the input_scenarios array. 
By running the orchestratorscenario.py file, you will execute and save the results
of VICRes-DCASCADE 16 times, once for each input scenario (Results folder from 1 to 16).

@author: bruno
"""
import sys
import numpy as np
import scipy.io
from Fi_extraction import Fi_extraction
from initialize_sed_input import initialize_sed_input
from Engelund_Hansen_tr_cap import Engelund_Hansen_tr_cap
from velocity_EH import velocity_EH
from sed_transfer_fast import sed_transfer_fast
from D_finder import D_finder
from changeFSLResVolume import changeFSLResVolume
from sockets import Driver, InterfaceSocket, Status
## Read index of the run scenario from the orchestrator
if len(sys.argv) > 1:
        indexscenario = int(sys.argv[1])-1

##  Create Socket Connection
## CASCADE is acting as a server, and it is opening a 
## socket interface with a given address, to which VICRes will connect
dtype = np.float32
server = InterfaceSocket()
server.open()
client, address = server.server.accept()
# CASCADE is waiting for VICRes connection
driver = Driver(client)
print(" @SOCKET:   Client asked for connection from " + str(address) + ". Now hand-shaking.")
# The server (CASCADE) and the client (VICRes) are connected, and they can now exchane information

##  Read input

## All the first steps concern the reading of the input file from the .mat structure: 
## it's very inefficient but everything can easily adapted to read the same information from a shapefile
# Load the network data
data = scipy.io.loadmat('network_data_3Snew.mat')

# Extract the variables from the network structure
network = data['Network'][0, 0]
reach_data = data['ReachData']
psi_3S = data['psi_3S'][0]
n_reaches = len(reach_data)
n_classes = len(psi_3S)

# Convert the MATLAB struct arrays to Python lists
NH = network['NH']
NH = NH.flatten()
NH=NH-1 #These are the reaches ID but they are used as indices and in python everything start from zero
Wa = [reach_data[i]['Wac'] for i in range(len(reach_data))]
Wa = np.array([d[0] for d in Wa])
Wa = Wa.flatten()
Wac=np.copy(Wa)
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
sed_yieldreach = [reach_data[i]['sed_yield'] for i in range(len(reach_data))]    
sed_yieldreach = np.array([d[0] for d in sed_yieldreach])
sed_yieldreach = sed_yieldreach.flatten()
directAd = [reach_data[i]['directAd'] for i in range(len(reach_data))]
directAd = np.array([d[0] for d in directAd])
directAd = directAd.flatten()
II = network['II']
outlet = np.where(fromN == toN)[0]
el_node = np.concatenate((el_FN, el_TN[outlet]))
Node_el = el_node

# Computing constant (PSI) for flushing transport capacity
dmi_finer = 2.0**(-psi_3S)
PSI = np.where(dmi_finer < 0.016, 1600, np.where((dmi_finer >= 0.016) & (dmi_finer < 0.1), 650, 300)) 
fine_sand_indices = 2**(-psi_3S) < 0.125
coarse_sand_indices = 2**(-psi_3S) > 1
medium_indices = ~fine_sand_indices & ~coarse_sand_indices

# Defining sediment density [kg/m3]
rho_s = 2650 

# Defining Reservoir Parameters
numres=109
TEFlag = 1
sh = (1,n_reaches+numres)
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

## Load Historical Discharge Data
q_data =scipy.io.loadmat('Q_Historical.mat')
Qall = q_data['Q']
dates_Q =scipy.io.loadmat('dates.mat')
dates_Q=dates_Q['dates_Q']
# Time horizon Inizialization
timescale = len(Qall)

## define sediment budget 

# sed_yield and sed_D50 are the inputs for the sediment budget definition (in this case we have 16 acceptable scenarios, for each scenario we run VICRes-CASCADE and we save the results)
input_scenarios = [
        [55, 0.25],
        [60, 0.25],
        [30, 0.1],
        [35, 0.1],
        [40, 0.1],
        [45, 0.1],
        [50, 0.1],
        [30, 0.075],
        [35, 0.075],
        [40, 0.075],
        [45, 0.075],
        [25, 0.05],
        [30, 0.05],
        [35, 0.05],
        [20, 0.025],
        [25, 0.025]
    ]
sed_yield = input_scenarios[indexscenario][0]
sed_D50 = input_scenarios[indexscenario][1]

#sed_input_param = [sed_yield, sed_D50]

# Initialize output elements
sed_input_param_comb = np.array(np.meshgrid(sed_yield, sed_D50)).T.reshape(-1, 2).T

#select is of input sediment budget (in our case there is only one input)
sb=0
sed_input_param=sed_input_param_comb[:, sb]

# define input sediment load
# define magnitude of the input sed yield, in t/km2/y
tot_yield = sed_input_param[0] * 1e6  # total yield for the 3S, in t/y
yield_km_y = tot_yield / np.sum(sed_yieldreach / np.max(sed_yieldreach) * directAd)
sed_yieldreach = sed_yieldreach / np.max(sed_yieldreach) * yield_km_y

# define GSD of the input sed yield
Fi_yield = np.tile(Fi_extraction(sed_input_param[1] / 1000, 0.8, psi_3S), (len(NH), 1))


# calculate input sed. yield
Qbi_input = initialize_sed_input(reach_data, Fi_yield, Qall, dates_Q, psi_3S, 10, sed_yieldreach, directAd)
last_matrix = Qbi_input[-1]
Qbi_input.append(last_matrix)
#np.save('Qbi_input.npy', Qbi_input)
#Qbi_input = np.load('Qbi_input1.npy')

# initialize variables
# sediment velocity parameters
phi = 0.4
minvel = np.min(Slope)

## Routing scheme

# initialize new variables

rows, col = II.shape

# Reservoir Variables
ResVolume = np.zeros((timescale + 1, numres))
sed_Volume = np.zeros((timescale + 1, numres))
sed_storage =  np.zeros((timescale + 1, numres))

Qbi_tr = [np.zeros((rows,col, len(psi_3S))) for _ in range(2)] #Qbi_tr report the sediment mobilized present in the reach AFTER transfer. If you want to save daily information of Qbi_tr, put for _ in range(timescale), and replace all the Qbi_tr[0] with Qbi_tr[t-1] and all the Qbi_tr[1] with Qbi_tr[1] 
Qbi_mob = [np.zeros((rows,col, len(psi_3S))) for _ in range(1)] #Qbi_mob report the sediment mobilized present in the reach BEFORE transfer. If you want to save daily information of Qbi_mob, put for _ in range(timescale), and replace all the Qbi_mob[0] with Qbi_mob[t]
Qbi_dep = [np.zeros((rows,col, len(psi_3S))) for _ in range(2)] # If you want to save daily information of Qbi_dep, put for _ in range(timescale), and replace all the Qbi_dep[0] with Qbi_dep[t-1] and all the Qbi_dep[1] with Qbi_dep[1]
Q_out = [np.zeros((len(II), len(psi_3S))) for _ in range(timescale+1)]

# Routing scheme
setout=np.zeros((len(II), 5))

Q = np.zeros(sh, dtype)
startflush = np.zeros(numres,dtype=int)
startsluicing = np.zeros(numres,dtype=int)
reserlength = np.zeros(numres,dtype=int)

for t in range(0, 2*(timescale)):
        # For each timestep, it needs to receive dischrge data from VICRes and to send back reservoir deposition volumes
        stat = driver.get_status()
        if stat == Status.Up | Status.NeedsInit:
            driver.initialise()
            print('Driver initialised.')
        if stat == Status.Up | Status.HasData:
            # Reading Discharge (and Reservoir Volumes) sent by VICres
            Q = driver.get_data(Q)
            Q1=Q[0,0:n_reaches]    
            Q2=np.copy(Q1)
            Qfinal=Q2[n_reaches-1]
            for l in range (75,n_reaches):
                Q1[l]=Q2[l-1]
            Q1[74]=Qfinal
            
            vol=Q[0,n_reaches:n_reaches+numres]
            # From the received volumes, it can understand whether a reservoir in VICRes is starting/finishing Flushing or SLuicing
            vol[np.isclose(vol, -1, atol=0.001)] = -1
            vol[np.isclose(vol, -2, atol=0.001)] = -2
            vol[np.isclose(vol, 0.001, atol=0.0001)] = 0.001
            reserflushed = np.where(abs(vol + 1) < 0.0001)[0] #index of the reservoir
            resersluiced = np.where(vol < -3)[0] #index of the reservoir
            startflush[reserflushed] = 1 #Determine the reservoir that is starting sediment management
            startsluicing[resersluiced] = 1 
            reserlength[resersluiced] = -vol[resersluiced]   
            finishedflush = np.where(abs(vol + 2) < 0.0001)[0] #index of reservoir that finishes flushing
            startflush[finishedflush] = 0 #flushing is finished for that reservoir
            startsluicing[finishedflush] = 0 #Sluicing is finished
            #From VICRes the volume is in 1000 m3, here we need it in m3
            ResVolume[round(t/2+0.2),:]=vol*1000 #m3
            
            # At each new timestep, we put back to zero Qbi_tr[1] Qbi_dep[1] and Qbi_mob[0], while the information of
            # Qbi_tr and Qbi_dep of the previous step remains stored in Qbi_tr[0] and Qbi_dep[0]
            Qbi_tr[1] = np.zeros(Qbi_tr[0].shape)
            Qbi_dep[1] = np.zeros(Qbi_dep[0].shape)
            Qbi_mob[0] = np.zeros(Qbi_mob[0].shape)
            Q_out[round(t/2+0.2)] = np.zeros(Q_out[0].shape)
            
            # calculate new water depth for all reaches
            # Manning
            h = (Q1[:] * n_Man / (Wac * np.sqrt(Slope))) ** (3 / 5)
            v = 1 / n_Man * h ** (2 / 3) * np.sqrt(Slope)
                  
            for n in NH:
                # extract the deposit layer of the reach from the relative cell in the previous timestep
                V_dep_old = np.squeeze(Qbi_dep[0][:, n, :])

                # 1) extract the deposit layer from the storage matrix and load the incoming cascades
                Qbi_incoming = np.squeeze(Qbi_tr[0][:, n, :])
                Qbi_incoming[n, :] = Qbi_input[round(t/2+0.2)][n, :]
                

                # 2) find cascades to be included into the active layer according to the limit V_lim_tot,
                # and use the cumulative GSD to compute tr_cap
                # find total sediment volume GSD
                Fi_r_reach = np.sum((V_dep_old + Qbi_incoming), axis=0) / np.sum((V_dep_old + Qbi_incoming)) #i find the GSD of the active layer, for the transport capacity calculation
                Fi_r_reach[np.isnan(Fi_r_reach)] = 0 #if V_act is empty, i put Fi_r equal to 0 for all classes

                # calculate transport capacity using the Fi of the active layer, the resulting tr_cap is in m3/s and is converted in m3/day
                # in case the active layer is empty, the tr_cap is equal to 0
                tr_cap = Engelund_Hansen_tr_cap(Fi_r_reach, Slope[n], Wac[n], v[n], h[n], psi_3S) * 24 * 60 * 60

                #If in a reach we have a reservoir, and in these timestep is performing flushing, and there are the hydraulics conditions to have effective flushng in this particular timestep
                # We use atkinson formula to compute the flushing transport capacty
                if ResFlag[n] == 1 and abs(vol[ResID[n]-1] + 1) > 0.0001 and vol[ResID[n]-1]>0.0001 and startflush[ResID[n]-1]==1:
                    tr_cap = PSI*vol[ResID[n]-1] #t/s
                    tr_cap = tr_cap * 1000 / rho_s  #m3/s
                    tr_cap = tr_cap * Fi_r_reach* 24 * 60 * 60 #m3/d
                    
                    
                # If in the reach n we have a reservoirs, the incoming volume is reduced with the Brune Curve (if there is not a flushing event in these days)
                if ResFlag[n] == 1 and abs(vol[ResID[n]-1] - 0.001) > 0.0001 and startflush[ResID[n]-1]==0:   
                    if TEFlag == 1: 
                        # V_mob contains the volume mobilized in the reach, that is about to be transferred downstream, without the presence of the reservoir
                        #Q in m3/s, volume in m3, residence time in d^-1
                        Qday=Q1[n]*3600*24
                        if Qday == 0:
                            Qday = 0.0001
                        ResidenceTime = ResVolume[round(t/2+0.2),ResID[n]-1]/Qday
                        if ResidenceTime == 0:
                            ResidenceTime = 0.0001
                        Te=np.zeros(len(psi_3S))
                        # # # Fine Sand or finer: Gill equation for fine grained sediments (< 0.125 mm)
                        Te[fine_sand_indices] = (ResidenceTime**(3))/(1.0265*ResidenceTime**(3)+0.02621*ResidenceTime**(2)-0.133*10**(-3)*ResidenceTime+0.1*10**(-5))
                        
                        # # # Coarse Sand or Gravel: Jothiprakash and Grag (2008), coarse sediment (> 1 mm)
                        Te[coarse_sand_indices] = (8000 - 36 * (ResidenceTime**(-0.78))) / (78.85 + (ResidenceTime**(-0.78)))
                        Te[coarse_sand_indices] /= 100
                        
                        # # # Medium: Jothiprakash and Grag (2008), medium sediment (0.125 to 1 mm)
                        Te[medium_indices] = ResidenceTime / (0.00013 + 0.01 * ResidenceTime + 0.0000166 * (ResidenceTime**0.5))
                        Te[medium_indices] /= 100

                        Te[Te < 0] = 0
                        Te[Te > 1] = 1

                        # Handle NaN values
                        Te[np.isnan(Te)] = 0
                        
                        # If in the reservoirs we are having a sluicing event, we use the churcill curve to compute the TE
                        if startsluicing[ResID[n]-1]==1:
                            si=((ResVolume[round(t/2+0.2),ResID[n]-1]/Q1[n])**2)/reserlength[ResID[n]-1]
                            Te[:] = (112 - 800 * (0.3048 * si) ** (-0.2)) / 100
                            Te[Te < 0] = 0
                            Te[Te > 1] = 1
                            Te[np.isnan(Te)] = 0
                        # We reduce the incoming volume with the TE of the reservoir
                        V_mob = Qbi_incoming * (1-Te)
                        V_dep = np.maximum(0, (V_dep_old + Qbi_incoming - V_mob))
                        # Qbi_mob contains the volume mobilized in the reach, that is about to be transferred downstream
                        Qbi_mob[0][:, n, :] = V_mob
                        Qbi_dep[1][:, n, :] = V_dep
                    
                # 3) Deposit the cascades in the active layer until the volume mobilized for each class is equal to the tr_cap                                                
                else:
                    prc = np.minimum(tr_cap / np.sum((V_dep_old + Qbi_incoming), axis=0), 1)
                    prc[np.isnan(prc)] = 0
                    V_mob = (V_dep_old + Qbi_incoming) * prc
                    V_dep = np.maximum(0, (V_dep_old + Qbi_incoming - V_mob))
                    Qbi_dep[1][:, n, :] = V_dep
                    # Qbi_mob contains the volume mobilized in the reach, that is about to be transferred downstream
                    Qbi_mob[0][:, n, :] = V_mob

            # end of the reach loop
            # 5) Move the mobilized volumes to the destination reaches according to the sediment velocity
            v_sed = velocity_EH(Slope, Wac, v, h, minvel, phi, psi_3S)
            #If we have a reservoir, the velocity of the reach should be almost 0
            for n in NH:
                if ResFlag[n]==1 and abs(vol[ResID[n]-1] - 0.001) > 0.0001:
                    if startflush[ResID[n]-1]==1 and abs(vol[ResID[n]-1] + 1) > 0.0001: #flushing is active and drawdown is terminated
                        v_sed[:,n]=v_sed[:,n]
                    else:
                        v_sed[:,n]=(1-0.99)*v_sed[:,n]
    
            # loop for all reaches, now that we have the Fi_r and thus can compute transfer rates for all reaches
            for n in NH:
                V_mob = np.squeeze(Qbi_mob[0][:, n, :])

                if np.sum(V_mob) > 0:
                    reach_dest, setout[n,:] = sed_transfer_fast(n, v_sed * (60 * 60 * 24), Lngt, network, NH)
                    # i sum the volumes transported in reach n with all the other
                    # volumes mobilized by all the other reaches in time t
                    for c in range(len(psi_3S)):
                        if setout[n,c] == 0: #if the volume does not leave throught the outlet
                            Qbi_tr[1][:, reach_dest[c], c] += V_mob[:, c]
                        else:
                            Q_out[round(t/2+0.2)][:, c] =Q_out[round(t/2+0.2)][:, c] + V_mob[:, c]
            # The computation of Qbi_dep and Qbi_tr for this timestep is finished
            # Before proceeding to the next timestep, we save the current values of Qbi_dep and Qbi_tr on  Qbi_dep[0] and  Qbi_tr[0]
            # Remove these two lines if Qbi_dep and Qbi_tr have dimension timescale and not 2
            Qbi_dep[0] = Qbi_dep[1]
            Qbi_tr[0] = Qbi_tr[1]
            Qbi_tr_t = np.squeeze(np.sum(Qbi_tr[1], axis=0))
            Qbi_dep_t = np.squeeze(np.sum(Qbi_dep[1], axis=0))
            # As last step, we need to compute the deosited volumes in each reservoir, and give this information back to VICRes
            sed_Volume[round(t/2+0.2), :], sed_storage[round(t/2+0.2), :] = changeFSLResVolume(Qbi_tr_t, Qbi_dep_t, NH,ResFlag,ResID,numres, phi, vol)
            sed_Volume[np.isnan(sed_Volume)] = 0
            sed_storage[np.isnan(sed_storage)] = 0
                    
        #Each timestep is divided in two parts, firstly we read discharge data and simulate sediment transport, secondly we give reservoir deposition feedback to VICRes  
        elif stat == Status.Up | Status.Ready:
            b = np.zeros(numres, dtype=vol.dtype)
            for i in range(numres+1):
                if i in ResID:
                    b[i-1] = sed_Volume[round(t/2+0.2), i-1]#1000 m3
            #print(b)
            driver.send_data(b)

## output processing
outcum_tot = np.array([np.sum(Q_out[t], axis=0).sum() for t in range(timescale)])
a = np.zeros(timescale // 365)
tot_sed_year = np.zeros(timescale // 365)
D50_year = np.zeros(timescale // 365)
Fi_year = np.zeros((len(psi_3S), timescale // 365))

for i in range(timescale // 365):
    a[i] = np.sum(outcum_tot[365 * (i):365 * (i+1) ])
    tot_sed_year[i] = a[i] * 2.6 / 1000000
    D50_year[i] = D_finder(np.sum([np.sum(Q_out[j], axis=0) for j in range( 365 * (i), 365 * (i+1) )], axis=0) /
                            a[i], psi_3S, 50) * 1000
    Fi_year[:, i] = np.sum([np.sum(Q_out[j], axis=0) for j in range( 365 * (i), 365 *(i+1) )], axis=0) / a[i]

data_plot_out = Q_out
distanceupstream=network['Upstream']['Distance'][0]
distanceupstream=distanceupstream[0]
distanceupstream=np.squeeze(distanceupstream)
river_reach = {
    1: {'name': 'Se Kong', 'index': np.where(distanceupstream[73] != np.inf)[1].tolist() + [74]},
    2: {'name': 'Se San', 'index': np.where(distanceupstream[161] != np.inf)[1].tolist()},
    3: {'name': 'Sre Pok', 'index': np.where(distanceupstream[242] != np.inf)[1].tolist()},
    4: {'name': 'Se San - Sre Pok confluence', 'index': list(range(162, 167))}
}

outcum_provenance = np.zeros_like(data_plot_out[0])
#data_plot_out.append(data_plot_out[0])

for n in range(data_plot_out[0].shape[0]):
    outcum_provenance_n_year = np.zeros((timescale // 365, data_plot_out[0].shape[1]))
    for y in range(timescale // 365):
        for t in range(365):
            outcum_provenance_n_year[y, :] += data_plot_out[365 * y + t][n, :]
    outcum_provenance[n, :] = np.round(np.mean(outcum_provenance_n_year, axis=0))

outcum_provenance_river = np.zeros((len(river_reach) + 1, data_plot_out[0].shape[1]))

for i in range(len(river_reach)):
    outcum_provenance_river[i, :] = np.sum(outcum_provenance[river_reach[i+1]['index'], :], axis=0)

outcum_provenance_river[-1, :] = np.sum(outcum_provenance_river[0:len(river_reach), :], axis=0)

## save results

file_path = f"../VICRes/RoutingSetup3S/Results{indexscenario + 1}/outcumprovenance.csv"
np.savetxt(file_path, outcum_provenance, delimiter=',')

file_path = f"../VICRes/RoutingSetup3S/Results{indexscenario + 1}/outcumprovenanceyear.csv"
np.savetxt(file_path, outcum_provenance_n_year, delimiter=',')

file_path = f"../VICRes/RoutingSetup3S/Results{indexscenario + 1}/outcumprovenanceriver.csv"
np.savetxt(file_path, outcum_provenance_river, delimiter=',')

file_path = f"../VICRes/RoutingSetup3S/Results{indexscenario + 1}/totsedyear.csv"

np.savetxt(file_path, tot_sed_year, delimiter=',')

file_path = f"../VICRes/RoutingSetup3S/Results{indexscenario + 1}/D50year.csv"

np.savetxt(file_path, D50_year, delimiter=',')

file_path = f"../VICRes/RoutingSetup3S/Results{indexscenario + 1}/Fiyear.csv"

np.savetxt(file_path, Fi_year, delimiter=',')
# output struct definition
data_output = {
    'outcum_provenance_river': outcum_provenance_river,
    'tot_sed_year': tot_sed_year,
    'D50_year': D50_year,
    'Fi_year': Fi_year,
    'Fi_input': Fi_extraction(sed_input_param[1] / 1000, 0.8, psi_3S)
}

