import os
import subprocess
import time
import sys

def run_main_and_rout():
    # Run VICRes_CASCADE. If the optimization is run on multiple cores, it will create a socket connection for each of the model run in parallel.
    port=int(sys.argv[1])
    cascadeprocess={}
    rout_process={}
    directory="CASCADE"
    os.chdir(directory)
    cascadeprocess[port]=subprocess.Popen(['python3', 'DCASCADE_3SOpt.py', str(port)])
    time.sleep(10)
    os.chdir("..")
    directory=f"VICRes/RoutingSetupCores/RoutingSetup3S{port}"
    os.chdir(directory)
    new_config_file_name = f'configuration{port}.txt'

    ## Read the contents of the file configuration.txt
    with open('configuration.txt', 'r') as config_file:
        lines = config_file.readlines()

    # Add the value of port to the last line
    lines.append(f'{port}\n')

    # Create a new file with the name based on the variable port
    with open(new_config_file_name, 'w') as new_config_file:
        new_config_file.writelines(lines) 
    os.chdir("../..")
    os.chdir("RoutingOpt/SourceCode")   
    #os.system('./rout ../../RoutingSetup3S/configuration{}.txt'.format(port))
    rout_process[port] = subprocess.Popen(['./rout', f'../../RoutingSetupCores/RoutingSetup3S{port}/configuration{port}.txt'])
    
    # Wait for both processes to complete
    cascadeprocess[port].wait()
    rout_process[port].wait()
# Call the function to run main.py and rout.exe sequentially
run_main_and_rout()

