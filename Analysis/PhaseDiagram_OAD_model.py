# 

import numpy as np
from lib.gillespie import (gillespie, gillespie_fraction)
from lib.misc import (load_topology, sample_graph_configuration_model,make_dirs)

import multiprocessing as mp
import os

# PAckages for profiling
import time 
import cProfile
import pstats




def analyze(r0,r1_list,filepath_output):
    """ Simulate the reaction-diffusion over the graph G 
        and save the results in a text file
    """
    FNAME_OUTPUT = filepath_output + f"{name_topology}_r0_{r0:.5f}.dat"
    with open(FNAME_OUTPUT,"w+") as file:
        file.write(f"topology: {name_topology}\nnumber of nodes: {NUM_NODES}\nnumber of samples: {NUM_SAMPLES}\nmax time step: {MAX_TSTEP}\ninitial fraction disrupted nodes: {dynp_pINI}\n")

    for r1 in r1_list:
        print(f"analyzing R0: {r0:2.3f}, R1: {r1:2.3f}")
        O = gillespie_fraction([G,r0,r1,MAX_TSTEP,NUM_SAMPLES,dynp_pINI])
        with open(FNAME_OUTPUT,"a+") as file:
            # file.write(str(r1) + "\t" + str(O) + "\n")
            file.write(f"{r1:.6f}"+"\t"+f"{O:.6f}"+"\n") 
        if O < 1e-3: # 
            break   


#########################################################
#########################################################
#########################################################

name_topology   = "europe" # america europe airports random


NUM_NODES    = 5000   # 1000
NUM_SAMPLES  = 150     # 100
MAX_TSTEP    = 2000   # 1000
dynp_pINI    = 0.0005 # Fraction of disrupted nodes on the network as initial condition (default 0.05)

G, NUM_NODES = load_topology(name_topology,NUM_NODES)

filepath_output = f"./Output_OAD/{name_topology}_nodes_{NUM_NODES}_samples_{NUM_SAMPLES}_maxtime_{MAX_TSTEP}/"

make_dirs(filepath_output)



r0_list = np.linspace(0,0.3,50)
r1_list = np.linspace(0,1.1*(1/dynp_pINI),50)
r1_list = np.linspace(0,4,50)



def apply_async_with_callback():
    """ Parallelize the execution of the function analyze """
    pool = mp.Pool(8)
    for i in range(len(r0_list)):
        pool.apply_async(analyze, args = (r0_list[i], r1_list, filepath_output, ))
        # analyze(r0_list[i], r1_list, filepath_output)
    pool.close()
    pool.join()

if __name__ == '__main__':
    apply_async_with_callback()


# if __name__ == '__main__':
#     with cProfile.Profile() as profile:
#         apply_async_with_callback()
#     results = pstats.Stats(profile)
#     results.sort_stats(pstats.SortKey.TIME)
#     results.print_stats()
#     results.dump_stats("results.profile")
