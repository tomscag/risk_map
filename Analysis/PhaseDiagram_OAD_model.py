# 

import numpy as np
from lib.gillespie import (gillespie_optimized, gillespie_optimized_fraction)
from lib.misc import (load_topology, make_dirs)
from lib.input import Inputs

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
        O = gillespie_optimized_fraction([G,r0,r1,MAX_TSTEP,NUM_SAMPLES,dynp_pINI])
        with open(FNAME_OUTPUT,"a+") as file:
            # file.write(str(r1) + "\t" + str(O) + "\n")
            file.write(f"{r1:.6f}"+"\t"+f"{O:.6f}"+"\n") 
        if O < 1e-2: # 
            break   


def apply_async_with_callback():
    """ Parallelize the execution of the function analyze """
    pool = mp.Pool(8)
    for i in range(len(r0_list)):
        pool.apply_async(analyze, args = (r0_list[i], r1_list, filepath_output, ))
        # analyze(r0_list[i], r1_list, filepath_output)
    pool.close()
    pool.join()


#########################################################
#########################################################
#########################################################





if __name__ == '__main__':

    name_topology = "america" # america europe airports random
    inputs_dct = getattr(Inputs(), name_topology)
    NUM_NODES    = inputs_dct['num_nodes']   # 1000
    NUM_SAMPLES  = inputs_dct['num_samples'] # 100
    MAX_TSTEP    = inputs_dct['max_tstep']   # 1000
    dynp_pINI    = inputs_dct['init0']  # 0.0005 # Fraction of disrupted nodes on the network as initial condition (default 0.05)

    G  = load_topology(name_topology,NUM_NODES)

    
    filepath_output = make_dirs(f"./Output_OAD/{name_topology}_nodes_{NUM_NODES}_samples_{NUM_SAMPLES}_maxtime_{MAX_TSTEP}")


    r0_list = inputs_dct["r0_list"]
    r1_list = inputs_dct["r1_list"]

    apply_async_with_callback()


# if __name__ == '__main__':
#     with cProfile.Profile() as profile:
#         apply_async_with_callback()
#     results = pstats.Stats(profile)
#     results.sort_stats(pstats.SortKey.TIME)
#     results.print_stats()
#     results.dump_stats("results.profile")
