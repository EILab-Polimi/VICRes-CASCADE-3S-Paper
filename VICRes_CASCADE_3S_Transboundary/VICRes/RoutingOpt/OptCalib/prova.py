from platypus import EpsNSGAII, Problem, Real, Integer, ProcessPoolEvaluator, Hypervolume, nondominated
import os
import csv
import time
import numpy as np
import multiprocessing
import datetime
import math
import subprocess
from datetime import timedelta
from concurrent.futures import ProcessPoolExecutor

if __name__ == '__main__':
    start_time = time.time()
    number_of_days = 5479
    spinning_period = 0
    maximum_no_reservoirs = 9
    number_of_cores = 20
    reservoirs = [[0 for x in range(3)] for y in range(maximum_no_reservoirs)]  
    
    reservoirs[0][0] = 7
    reservoirs[1][0] = 11
    reservoirs[2][0] = 8
    reservoirs[3][0] = 109
    reservoirs[4][0] = 2
    reservoirs[5][0] = 3
    reservoirs[6][0] = 6
    reservoirs[7][0] = 22
    reservoirs[8][0] = 13

    
    
    # Setup intial conditions
    constraint = []
    num_var = 0
    os.chdir('../../ReservoirsCompromizeFlushing')
    for i in range(maximum_no_reservoirs):
        text_file = open('res'+str(reservoirs[i][0])+'.txt','r')
        lines = text_file.read().split('\n')
        hmax = float(lines[1].split('\t')[0])
        hmin = float(lines[1].split('\t')[1])
        opt = int(lines[8].split('\t')[0])
        vcap = float(lines[1].split('\t')[2])							#vcapacity
        vd = float(lines[1].split('\t')[3])								#vdead
        text_file.close()
        if (opt==1):													#simplified rule curve (OP1)
           constraint.append(Real(hmin,hmax))							#h1
           constraint.append(Real(hmin,hmax))							#h2
           constraint.append(Real(1,365))								#t1
           constraint.append(Real(1,365))								#t2
           num_var+=4
           reservoirs[i][1] = 1
        elif (opt==2):													#rule curve (OP2)
            for j in range(12):
                constraint.append(Real(hmin,hmax))						#hi (i = 1 - 12)
            num_var+=12
            reservoirs[i][1] = 2
        elif (opt==3):													#operating rule (OP3)
            constraint.append(Real(0.0001,(math.pi/2-0.0001)))			#x1 (0-pi/2), but avoid 0 and pi (90 degree)
            constraint.append(Real(vd,vcap))								#x2 (range from vdead to vcapacity)
            constraint.append(Real(vd,vcap))								#x3 (range from vdead to vcapacity), check condition later
            constraint.append(Real(0.0001,(math.pi/2-0.0001)))			#x4 (0-pi/2), but avoid 0 and pi (90 degree)
            num_var+=4
            reservoirs[i][1] = 3
        elif (opt==5):													#operating rule (OP5 or OP3 12 demand values)
            for j in range (12):
               constraint.append(Real(0.0001,(math.pi/2-0.0001)))			#x1 (0-pi/2), but avoid 0 and pi (90 degree)
               constraint.append(Real(vd,vcap))								#x2 (range from vdead to vcapacity)
               constraint.append(Real(vd,vcap))								#x3 (range from vdead to vcapacity), check condition later
               constraint.append(Real(0.0001,(math.pi/2-0.0001)))			#x4 (0-pi/2), but avoid 0 and pi (90 degree)
            num_var+=48
            reservoirs[i][1] = 5
        sedopt = int(lines[11].split('\t')[0])
        if (sedopt==1):
           constraint.append(Integer(1,365))                           #Drawdown Date
           constraint.append(Real(1,5000))                         #Minimum Inflow to activate drawdown
           constraint.append(Integer(1,30))                            #Duration in days
           constraint.append(Integer(1,16))                             #Frequency in years
           num_var+=4
           reservoirs[i][2] = 1
        if (sedopt==2):
           constraint.append(Integer(1,365))                           #Drawdown Date
           constraint.append(Real(5,5000))                         #Minimum Inflow to activate drawdown
           constraint.append(Integer(30,150))                            #Duration in days
           constraint.append(Integer(1,16))                             #Frequency in years
           num_var+=4
           reservoirs[i][2] = 2
        if (sedopt==0):
           reservoirs[i][2] = 0

    # Reading Parametrization of optimal solutions
    # Reading optimization results for sluicing
    os.chdir('../OptimalSolutions')
    filename1 = 'optimization_variables_flush_final1.txt'
    
    data1 = np.loadtxt(filename1)
    
    all_data1 = data1

    
    varnew = all_data1

    
    varnew_unique = np.unique(varnew, axis=0)
    

    
    print("Lunghezza di varnew_unique:", len(varnew))

    
    varnew_unique = np.unique(varnew, axis=0)
    # Salvataggio di varnew_unique in un file di testo
    
    
    print("Lunghezza di varnew_unique:", len(varnew_unique))

    a = fnew_unique[:, 0].copy()
    b = fnew_unique[:, 1].copy()
    fnew_unique[:, 0] = fnew_unique[:, 2]
    fnew_unique[:, 1] = a
    fnew_unique[:, 2] = b

    # Cap values for each indicator (the ones used to normalize the values of each objective between 0 and 1 during the optimization)
    maxtotalpro = 10627501.095537942 #hp production
    maxfirm = 3.207057309002500e+04 #firm hp
    maxsed = 30

    fnew_unique_flush = fnew_unique
    fnew_unique_flush[:, 0] = (1 - fnew_unique_flush[:, 0]) * maxsed;
    fnew_unique_flush[:, 1] = (1 - fnew_unique_flush[:, 1]) * maxtotalpro*10**(-6) * 24 / 15
    fnew_unique_flush[:, 2] = (1 - fnew_unique_flush[:, 2]) * maxfirm*24*10**(-3)


    os.chdir('../RoutingSetup3S_Opt/Results')
    # Call optimization function
    
    
    start = datetime.datetime.now()
    #results = parallel_execution(varnew_unique, number_of_cores)
    print("Risultati finali:", results)
    end = datetime.datetime.now()


    
    



    # END OF FILE