import os
import subprocess
import time
import sys
from datetime import timedelta

def run_main_and_rout():
    # Run VICRes_CASCADE 16 times, with the 16 different defined sediment budget 
    for i in range(1,17):
        os.chdir("CASCADE")
        # Run DCASCADE_3S from terminal
        cascadeprocess=subprocess.Popen(['python3', f'DCASCADE_3S.py', str(i)])
        time.sleep(10)
        
        # Run rout.exe from terminal (VICRes)
        os.chdir("../VICRes/Routing/SourceCode")
        rout_process=subprocess.Popen(['./rout', f'../../RoutingSetup3S/configuration{i}.txt'])
        
        cascadeprocess.wait()
        rout_process.wait()
        os.chdir("../../..")

# Call the function to run DCASCADE_3S.py and rout.exe sequentially
run_main_and_rout()

