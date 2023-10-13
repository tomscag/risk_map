#

import numpy as np
import networkx as nx
from Functions import simulate_reaction_diffusion_gillespie_fraction

import multiprocessing as mp

# PAckages for profiling
import time 
import cProfile
import pstats



def load_fully_connected_graph(NUM_NODES):
    return nx.complete_graph(NUM_NODES)


def analyze_full(r0,r1_list,filepath_output):
    """ Simulate the reaction diffusion over the graph G 
        and save the results in a text file
    """
    FNAME_OUTPUT = filepath_output + f"{name_topology}_r0_{r0*NUM_NODES:.5f}.dat"
    with open(FNAME_OUTPUT,"w+") as file:
        file.write(f"topology: {name_topology}\nnumber of nodes: {NUM_NODES}\nnumber of samples: {NUM_SAMPLES}\nmax time step: {MAX_TSTEP}\ninitial fraction disrupted nodes: {dynp_pINI}\n")

    for r1 in r1_list:
        print(f"analyzing R0: {r0*NUM_NODES:2.3f}, R1: {r1*NUM_NODES:2.3f}")
        O = simulate_reaction_diffusion_gillespie_fraction.main([G,r0,r1,MAX_TSTEP,NUM_SAMPLES,dynp_pINI])
        with open(FNAME_OUTPUT,"a+") as file:
            # file.write(str(r1) + "\t" + str(O) + "\n")
            file.write(f"{r1*NUM_NODES:.6f}"+"\t"+f"{O:.6f}"+"\n") 


########################################
########################################
########################################


name_topology   = "full" # america europe airports random
filepath_output = "./Output_OAD/"

NUM_NODES    = 5000   # 1000
NUM_SAMPLES  = 50     # 100
MAX_TSTEP    = 1000   # 1000
dynp_pINI    = 0.05   # Fraction of disrupted nodes on the network as initial condition

G = load_fully_connected_graph(NUM_NODES)

r0_list = np.linspace(0,2,40)/NUM_NODES
r1_list = np.linspace(0,1.3*(1/dynp_pINI),40)/NUM_NODES





def apply_async_with_callback():
    """ Parallelize the execution of the function analyze """
    pool = mp.Pool(8)
    for i in range(len(r0_list)):
        pool.apply_async(analyze_full, args = (r0_list[i], r1_list, filepath_output, ))
        # analyze(r0_list[i], r1_list, filepath_output)
    pool.close()
    pool.join()

if __name__ == '__main__':
    apply_async_with_callback()