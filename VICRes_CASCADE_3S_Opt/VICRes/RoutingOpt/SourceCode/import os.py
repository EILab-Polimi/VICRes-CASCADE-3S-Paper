#import os
#import tensorflow as tf
#os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


#working_directory = "'C:/Users/bruno/OneDrive - Politecnico di Milano/Tesi/VICRES-CASCADE/integration/vicres-socket/VICResOpt-master (2)/Routing/SourceCode'"
#command = "C:/cygwin64/bin/bash.exe -l -c cd 'C:/Users/bruno/OneDrive - Politecnico di Milano/Tesi/VICRES-CASCADE/integration/vicres-socket/VICResOpt-master (2)/Routing/SourceCode'  &&  provaa ../../routingsetup/configuration.txt"
#os.system(command)

import os
import subprocess
import time

def run_main_and_rout():
    # Run main.py
    subprocess.Popen(['python', 'python/main.py'])
    time.sleep(1)
    # Run rout.exe in Cygwin terminal
    #os.system('rout  ../../routingsetup/configuration.txt')

    cygwin_command = '''C:/cygwin64/bin/bash.exe -l -c "cd 'C:/Users/bruno/OneDrive - Politecnico di Milano/Tesi/Github fork/VICResOpt-master (2)/Routing/SourceCode'  &&  ./rout ../../routingsetup3S/configuration.txt"'''
    subprocess.call(cygwin_command)
    #time.sleep(30)

# Call the function to run main.py and rout.exe sequentially
run_main_and_rout()

