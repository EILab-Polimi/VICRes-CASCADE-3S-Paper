# coding: utf-8
# Code: Environmental Intelligence Lab / Polimi
# author: bruno
# ------------------------------------------------------------------------------------------------------------------------------------------------
'''
This code builds on the optimisation.py code from the Resilient Water Systems Group / Singapore University of Technology and Design and adds 
the capability to call the coupled version of CASCADE and VICRes together, as well as optimize the parameters that determine the release of 
water and sediments in VICRes. To consider and include the sediment size in the optimizations, an additional objective (obtained from CASCADE) 
has been added, which represents the annual amount of sediment that reaches the outlet section.
If you run it using multiple cores, make sure that you have prepared the right number of folders for running the coupled framework in parallel 
(one folder for each core, to create different socket connections and save the different results).
If you need further information to run the script for optiomization contact: bruno.invernizzi@polimi.it
'''
# ------------------------------------------------------------------------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------------------------------------------------------------------------
# THE FOLLOWING CODE CALLS CASCADE-VICRes TO OPTIMIZE RESERVOIRS RELEASE PARAMETERS CONSIDERING DIFFERENT OBJECTIVE

def viccall(vars):
    # Declare decision variables (VIC)
    global maximum_no_reservoirs
    global number_of_days
    global spinning_period
    global objs
    global observedpeak
    global floodthreshold
    global pathwaterdemand
    global pathreswl
    global reservoirs
    totalinstallcapcity = 0
    os.chdir('../../RoutingSetup3S_Opt')
    text_file = open('reservoiroptimization.txt','r')
    lines = text_file.read().split('\n')
    number_of_days = int(lines[2].split('\t')[0])
    spinning_period = int(lines[3].split('\t')[0])
    maximum_no_reservoirs = len(lines[11].split(' '))
    objs = [0 for x in range(7)]
    reservoirs = [[0 for x in range(2)] for y in range(maximum_no_reservoirs)]
    pathwaterdemand = lines[12].split('\t')[0]
    observedpeak = float(lines[13].split('\t')[0])
    floodthreshold = float(lines[14].split('\t')[0])
    pathreswl = lines[15].split('\t')[0]
    totalobs = 0
    minvar = []
    maxvar = []
    eps = []
    for i in range(7):
        objs[i] = int(lines[9].split(' ')[i])
        if (objs[i]>=1):
            totalobs+=1
            minvar.append(0)
            maxvar.append(1)											# we normalize objective functions to the range from 0 - 1
            eps.append(0.001)											# epsilon 
    for i in range(maximum_no_reservoirs):
        reservoirs[i][0] = int(lines[11].split(' ')[i])
    text_file.close()
    os.chdir('../RoutingSetup3S_Opt/Results')
    # Starting the process
    port=33000
    rank = multiprocessing.current_process()._identity[0]
    #print("-----------------------------------------------------------------------------------")
    #print("Thread no: ",rank," out of ",multiprocessing.cpu_count())
    #maximum_no_reservoirs=1
    # Modify reservoir parameters
    port=port+rank
    countrow = 0
    directory = f"../../ReservoirsCores/ReservoirsOpt{port}"
    os.chdir(directory)
    D5demand = [0.0 for x in range(12)]
    #time.sleep(rank*2)
    for i in range(maximum_no_reservoirs):
        text_file = open('res'+str(reservoirs[i][0])+'.txt','r')
        lines = text_file.read().split('\n')
        hmax = float(lines[1].split('\t')[0])
        hmin = float(lines[1].split('\t')[1])
        volume = float(lines[1].split('\t')[2])
        dvolume = float(lines[1].split('\t')[3])
        height = float(lines[1].split('\t')[4])
        discharge = float(lines[1].split('\t')[5])
        yearopt = int(lines[1].split('\t')[6])
        initialvolume = float(lines[1].split('\t')[7])
        resname = lines[1].split('\t')[8]
        seepage = float(lines[3].split('\t')[0])
        infil = float(lines[3].split('\t')[1])
        opt = int(lines[8].split('\t')[0])
        if (opt==3):
            D4demand = float(lines[6].split('\t')[0])
        elif (opt==5):
            for j in range(12):
                D5demand[j] = float(lines[6+j].split('\t')[0])
        sedopt = int(lines[11].split('\t')[0])
        if (sedopt==1):
            maxrate=float(lines[13].split('\t')[2])
            maxwl=float(lines[13].split('\t')[5])
            mindisch=float(lines[13].split('\t')[6])
            width=float(lines[13].split('\t')[7])
            slope=float(lines[13].split('\t')[8])
            length=float(lines[13].split('\t')[9])
            gatelevel=float(lines[13].split('\t')[10])
            bottomcap=float(lines[13].split('\t')[11])
        elif (sedopt==2):
            maxrate=float(lines[13].split('\t')[2])
            maxwl=float(lines[13].split('\t')[5])
            mindisch=float(lines[13].split('\t')[6])
            width=float(lines[13].split('\t')[7])
            slope=float(lines[13].split('\t')[8])
            length=float(lines[13].split('\t')[9])
            gatelevel=float(lines[13].split('\t')[10])
            bottomcap=float(lines[13].split('\t')[11])
        text_file.close()
        with open('res'+str(reservoirs[i][0])+'.txt','w') as my_csv:
            my_csv.write("Hmax(M)	Hmin(M)	Scap(1000M3)	Sd(1000M3)	Hturbine(M)	Qdesign (M3/s)	Year	Sinitial(1000M3)	Name\n")
            my_csv.write("%f\t%f\t%f\t%f\t%f\t%f\t%i\t%f\t%s\n"%(hmax,hmin,volume,dvolume,height,discharge,yearopt,initialvolume,resname))
            my_csv.write("SEEPAGE	INFILTRATION\n")
            my_csv.write("%f\t%f\n"%(seepage,infil))
            my_csv.write("IRRIGATION\n")
            my_csv.write("%d\n"%0)
            my_csv.write("DUMMY\n")
            my_csv.write("OPERATION STRATEGY\n")
            my_csv.write(str(opt)+"\n")
            if (opt==1):
                my_csv.write("%f\t%f\t%i\t%i\n"%(vars[countrow],vars[countrow+1],vars[countrow+2],vars[countrow+3]))
                countrow+=4
            elif (opt==2):
                my_csv.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(vars[countrow],vars[countrow+1],vars[countrow+2],vars[countrow+3],vars[countrow+4],
                vars[countrow+5],vars[countrow+6],vars[countrow+7],vars[countrow+8],vars[countrow+9],vars[countrow+10],vars[countrow+11]))
                countrow+=12
            elif (opt==3):
                x1 = vars[countrow]
                x2 = vars[countrow+1]
                x3 = vars[countrow+2]
                x4 = vars[countrow+3]
                if (x2>x3):
                    temp = x2
                    x2 = x3
                    x3 = temp
                my_csv.write("%f\t%f\t%f\t%f\t%f\n"%(D4demand,x1,x2,x3,x4))
                countrow+=4
            elif (opt==5):
                for j in range(12):
                    x1 = vars[countrow]
                    x2 = vars[countrow+1]
                    x3 = vars[countrow+2]
                    x4 = vars[countrow+3]
                    if (x2>x3):
                        temp = x2
                        x2 = x3
                        x3 = temp
                    my_csv.write("%f\t%f\t%f\t%f\t%f\n"%(D5demand[j],x1,x2,x3,x4))
                    countrow+=4
            my_csv.write("SEDIMENT STRATEGY\n")
            my_csv.write(str(sedopt)+"\n")
            my_csv.write("DrawdownDate  MinInflow(m3/s) MaxRate(m)  Duration(d) Frequency(y)    MAXWL(m)    MINDISCH(m3/s)  WIDTH(m)    SLOPE(-)    LENGTH(m)   GATELEVEL(m)    QDESIGN(m3/s)\n")
            if (sedopt==1):
                my_csv.write("%i\t%f\t%f\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(vars[countrow],vars[countrow+1],maxrate,vars[countrow+2],vars[countrow+3],maxwl,mindisch,width,slope,length,gatelevel,bottomcap))
                countrow+=4
            elif (sedopt==2):
                my_csv.write("%i\t%f\t%f\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(vars[countrow],vars[countrow+1],maxrate,vars[countrow+2],vars[countrow+3],maxwl,mindisch,width,slope,length,gatelevel,bottomcap))
                countrow+=4
        totalinstallcapcity+=0.9*9.81*height*discharge
    # 9.81 is the gravitational acceleration; 0.9 is turbine coefficient (change if needed but remember to change in routing model - reservoir.f)

    # Run routing model
    os.chdir('../../..')
    orche = subprocess.Popen(['python3', 'orchestratorOpt.py', str(port)])
    orche.wait()
    directory = f"VICRes/RoutingSetupCores/RoutingSetup3S{port}/Results"
    os.chdir(directory)
    
    text_file = open('date.day','r')
    lines = text_file.read().split('\n')
    countrow = 0
    violation = 0.0
    waterdeficit = 0.0
    peakds = 0.0
    waterdeviation = 0.0
    hydrofirm = totalinstallcapcity * 30 * 24 * 1000		 # In a month: 30 days x 24h x 1000 (MWh) - just for a cap value; not important
    VICdata = [[0 for x in range(4)] for y in range(number_of_days)]
    for line in lines:
        try:
            year,month,day,flows = filter(None,line.split(' '))
            VICdata[countrow][0] = float(year)
            VICdata[countrow][1] = float(month)
            VICdata[countrow][2] = float(day)
            VICdata[countrow][3] = float(flows)
    #         if (countrow>spinning_period):
    #             if (VICdata[countrow][3]>=floodthreshold):
    #                 violation+=1
    #             if (VICdata[countrow][3]>peakds):
    #                 peakds = VICdata[countrow][3]
            countrow+=1
        except:
            print("FINISH READING MODELED FILE ...")
    text_file.close()

    # Read water demand file
    #waterdemand= [[0 for x in range(maximum_no_reservoirs+1)] for y in range(number_of_days)]	# Read water demand file for water deficit calculation
    #text_file = open(pathwaterdemand,'r')
    #lines = text_file.read().split('\n')
    #countrow = 0
    #maxirr = 0.0
    #for line in lines:
    #    try:
    #        for i in range(maximum_no_reservoirs+1):								# water demand for n reservoir + basin outlet
    #            waterdemand[countrow][i] = float(filter(None,line.split('\t'))[i])
    #            if (countrow>spinning_period):
    #                maxirr+=waterdemand[countrow][i]
    #        countrow+=1
    #    except:
    #        print("...")
    #text_file.close()
    # reservoirwaterlevel = [[0 for x in range(maximum_no_reservoirs)] for y in range(365)] # Read predefined water levels for reservoirs
    # text_file = open(pathreswl,'r')
    # lines = text_file.read().split('\n')
    # maxreswl = 0
    # for line in lines:
    #     try:
    #        for k in range(365):
    #            for i in range(maximum_no_reservoirs):
    #                 reservoirwaterlevel[k][i] = float(filter(None,line.split('\t'))[i])
    #                 if (maxreswl<reservoirwaterlevel[k][i]):
    #                     maxreswl = reservoirwaterlevel[k][i]
    #     except:
    #        print("...")
    # text_file.close()

    # We normalize the flood peak by dividing the maximum flood peak value with a historical value (provided in the configuration file)
    # We want to minimize the flood peak
    # Unit: m
    peakds = peakds / observedpeak
    # We normalize the number of violation; this value can be normalized by dividing with the number of simulation days
    #but, in most of the region, there are rainny and dry seasons, so that we divide the number of violation days by (simulation days/2)
    # Unit: days
    violation = 1.0*violation/(countrow/2+1)
    totalproduction = 0.0
    temphydrofirm = [0 for x in range(int(number_of_days/365)*12)]
    countrow = 0
    occurences = [0 for x in range(int(number_of_days/365)*12)]
    # Read modeling results
    # - Calculate hydro power total production
    dischargers = [[0 for x in range(maximum_no_reservoirs)] for y in range(number_of_days)]
    for i in range(maximum_no_reservoirs):
        monthth = 0
        countrow = 0
        with open('reservoir_'+str(reservoirs[i][0])+'_75STA.day','r') as my_csv: 
            lines = my_csv.read().split('\n')
            for line in lines:
                try:
                    produce = 0
                    if (countrow>spinning_period):
                        list_so = list(filter(None,line.split()))
                        produce = float(list_so[6])
                        dischargers[countrow][0]=float(list_so[4]) + float(list_so[5])
                        #waterdeviation+=abs(float(list_so[1])-reservoirwaterlevel[int(countrow%365)][i])
                    if (countrow>92):
                        if VICdata[countrow][1] in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
                            if (VICdata[countrow][1]!=VICdata[countrow-1][1]):
                                monthth+=1
                                occurences[monthth]=0
                            temphydrofirm[monthth]+=produce
                            occurences[monthth]+=1
                    totalproduction+=produce
                    countrow+=1
                except:
                    countrow = countrow
    # We now normalize to the range from 0 - 1
    # Unit: MWh ---> we maximize hydropower production or minimize 1 - its normalized value
    #waterdeviation = waterdeviation/maxreswl/countrow      # small number, but fine
    waterdeviation=0.3
    for i in range(monthth):                            # ignore the last month which is normally in the wet season
        if  (temphydrofirm[i]<hydrofirm):
            hydrofirm=temphydrofirm[i]
            occ=occurences[i]
    #hydrofirm=1-hydrofirm/30/totalinstallcapcity*1000          # normalize the value smallest montlhy value of the firm hydropower
    hydrofirm = 1 - hydrofirm/32070.5731
    #totalproduction = 1 - totalproduction/totalinstallcapcity *1000/countrow
    totalproduction = 1 - totalproduction/10627501.095537942
    # - Calculate water deficit for the whole system
    countrow2 = 0
    # for i in range(spinning_period,countrow):
    #     try:
    #         releaseatoutlet = 0
    #         for k in range(maximum_no_reservoirs):
    #             waterdeficit+=max(0.0,waterdemand[i][k]-dischargers[i][k])
    #             releaseatoutlet+=min(waterdemand[i][i],dischargers[i][i])
    #         waterdeficit+=max(0.0,waterdemand[i][maximum_no_reservoirs]-VICdata[i][3]+releaseatoutlet)
    #         countrow2+=1
    #     except:
    #         print("...")
    # We minimize water deficit and normalize (Unit: m3/s)
    #waterdeficit=waterdeficit/maxirr							# annual water deficit + normalize. note: no need to devide by number of years because the normalization
    #Analyze result and return objective function values
    #New Sediment Indicator
    totalsed = 0
    count = 0
    pristinesed=30
    # file name with total sediment per year CSV
    file_name = "totsedyear.csv"

    with open(file_name, "r") as csv_file:
        csv_reader = csv.reader(csv_file)
    
        for row in csv_reader:
            try:
                sediment = float(row[0])
                totalsed += sediment
                count += 1
            except: 
                count=count
    # Compute the annual average and then normalize
    if count > 0:
        annualsedavg = totalsed / count
        sedload=1-annualsedavg/pristinesed
    else:
        print("The sediment file is empty.")

    os.chdir('..')
    
    
    
    if (violation>1):
        violation = 1
    elif (violation<0):
        violation = 0
    if (totalproduction>1):
        totalproduction = 1
    elif (totalproduction<0):
        totalproduction = 0
    if (waterdeficit>1):
        waterdeficit = 1
    elif (waterdeficit<0):
        waterdeficit = 0
    if (peakds>1):
        peakds = 1
    elif (peakds<0):
        peakds = 0
    if (hydrofirm>1):
        hydrofirm = 1
    elif (hydrofirm<0):
        hydrofirm = 0
    if (waterdeviation>1):
        waterdeviation = 1
    elif (waterdeviation<0):
        waterdeviation = 0
    if (sedload>1):
        sedload = 1
    elif (sedload<0):
        sedload = 0
    # with open("objectives.txt", "a") as file:
    #     file.write(f"{sedload} {totalproduction} {hydrofirm}\n")
    # vars_str = " ".join(map(str, vars))
    # # Scrivi la stringa nel file
    # with open("variables.txt", "a") as file:
    #     file.write(vars_str + "\n")
    
    returnobjs = []
    if (objs[0]>0):
        returnobjs.append(waterdeficit)
    if (objs[1]>0):
        returnobjs.append(totalproduction)
    if (objs[2]>0):
        returnobjs.append(hydrofirm)
    if (objs[3]>0):
        returnobjs.append(peakds)
    if (objs[4]>0):
        returnobjs.append(violation)
    if (objs[5]>0):
        returnobjs.append(waterdeviation)
    if (objs[6]>0):
        returnobjs.append(sedload)
    return returnobjs

# -----------------------------------------------------------------------------------------------------------------------------------
# Open operation file
if __name__ == '__main__':
    start_time = time.time()
    os.chdir('../../RoutingSetup3S_Opt')
    text_file = open('reservoiroptimization.txt','r')
    lines = text_file.read().split('\n')
    number_of_days = int(lines[2].split('\t')[0])
    spinning_period = int(lines[3].split('\t')[0])
    maximum_no_reservoirs = len(lines[11].split(' '))
    number_of_functions = int(lines[4].split('\t')[0])
    population = int(lines[5].split('\t')[0])
    number_of_cores = int(lines[6].split('\t')[0])
    objs = [0 for x in range(7)]
    reservoirs = [[0 for x in range(3)] for y in range(maximum_no_reservoirs)]
    pathwaterdemand = lines[12].split('\t')[0]
    observedpeak = float(lines[13].split('\t')[0])
    floodthreshold = float(lines[14].split('\t')[0])
    pathreswl = lines[15].split('\t')[0]
    totalobs = 0
    minvar = []
    maxvar = []
    eps = []
    for i in range(7):
        objs[i] = int(lines[9].split(' ')[i])
        if (objs[i]>=1):
            totalobs+=1
            minvar.append(0)
            maxvar.append(1)											# we normalize objective functions to the range from 0 - 1
            eps.append(0.01)											# epsilon 
    for i in range(maximum_no_reservoirs):
        reservoirs[i][0] = int(lines[11].split(' ')[i])
    text_file.close()
    
    # Setup intial conditions
    constraint = []
    num_var = 0
    os.chdir('../ReservoirsOptSluic')
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
    os.chdir('../RoutingSetup3S_Opt/Results')
    # Call optimization function
    
    problem = Problem(num_var, totalobs)
    problem.types[:] = constraint
    problem.function = viccall
    problem.constraints[:] = ">=0"
    start = datetime.datetime.now()
    hyp = Hypervolume(minimum=minvar, maximum=maxvar)
    allresults = []
    prova=[]
    provavar=[]
    alltrial=[]
    with ProcessPoolEvaluator(number_of_cores) as evaluator:
        algorithm = EpsNSGAII(problem, eps, population_size=population, evaluator = evaluator)
        while algorithm.nfe<number_of_functions:
            algorithm.step()
            one_step_result = hyp.calculate(algorithm.result)  # this param stores information of hypervolume indicator, save as a file if needed
            allresults.append(one_step_result)
            alltrial.extend(algorithm.result)
            print('NUMBER OF FUNCTION EVALUATIONS == ', algorithm.nfe)
    end = datetime.datetime.now()
    os.chdir('../Results')
    np.savetxt("allresults.txt",[s.objectives[:] for s in alltrial],fmt="%s")
    np.savetxt("allvariables.txt",[s.variables[:] for s in alltrial],fmt="%s")    
    nondominated_solutions = nondominated(alltrial)
    
    
    # Finish and save results to files
    print("Start",start)
    print("End",end)
    print("Finish running simulations! See opt_objectives.txt and opt_variables.txt for results.")
    np.savetxt("optimization_objectives.txt",[s.objectives[:] for s in nondominated_solutions],fmt="%s")
    
    np.savetxt("optimization_variables.txt",[s.variables[:] for s in nondominated_solutions],fmt="%s")
    np.savetxt("hypervolume.txt", [allresults],fmt="%s")
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Total Computational Time: {total_time} secondi")
    total_time_str = str(timedelta(seconds=total_time))

    # Stampa il tempo totale
    print(f"Total Computational Time: {total_time_str}")
   

    
    


    # END OF FILE
    # ------------------------------------------------------------------------------------------------------------------------------------------------ 