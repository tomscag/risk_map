# Single removal OAD

import numpy as np
from lib.gillespie import (gillespie_optimized)
from lib.misc import (load_topology, make_dirs)
from lib.input import Inputs

import networkx as nx
import multiprocessing as mp
import os
from collections import Counter

#################################
#################################
#################################

def read_network():
    G = nx.read_edgelist(fpathelist,nodetype=int)
    # Compute the largest connected component
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G0 = G.subgraph(Gcc[0]) 
    G0 = nx.convert_node_labels_to_integers(G0,first_label=0)
    return G0, Gcc[0]


def analyze(r0,r1_list,filepath_output):
    """ Simulate the reaction-diffusion over the graph G 
        and save the results in a text file
    """
    FNAME_OUTPUT = filepath_output + f"{name_topology}_r0_{r0:.5f}.dat"
    with open(FNAME_OUTPUT,"w+") as file:
        file.write(f"topology: {name_topology}\nnumber of nodes: {NUM_NODES}\nnumber of samples: {NUM_SAMPLES}\nmax time step: {MAX_TSTEP}\ninitial fraction disrupted nodes: {init0}\n")

    for r1 in r1_list:
        print(f"analyzing R0: {r0:2.3f}, R1: {r1:2.3f}")
        O = gillespie_optimized([G,r0,r1,MAX_TSTEP,NUM_SAMPLES,init0])
        with open(FNAME_OUTPUT,"a+") as file:
            # file.write(str(r1) + "\t" + str(O) + "\n")
            file.write(f"{r1:.6f}"+"\t"+f"{O:.6f}"+"\n") 
        if O < 1e-2: # 
            break   


def run_parallel(numcpu):
    """ Parallelize the execution of the function analyze """
    pool = mp.Pool(numcpu)
    for i in range(len(r0_list)):
        pool.apply_async(analyze, args = (r0_list[i], r1_list, filepath_output, ))
        # analyze(r0_list[i], r1_list, filepath_output)
    pool.close()
    pool.join()



#################################
#################################
#################################


if __name__ == "__main__":
    numcpu = 1
    num_nodes = 100  # Nodes 
    p_america = 0.00015 
    fpathelist = "./lib/motter/NetworkTopologies/network_america.edgelist"
    G,_  = read_network()
    # filepath_output = make_dirs(f"./Output_OAD/results_cluster_america_motter")

    N = 14990 
    num_run = int(1e2)
    
    MAX_TSTEP = 10^3
    NUM_SAMPLES = 1 # Samples per node
    
    r0 = 1e5*1/(N*p_america)
    r1 = 0e1
    dict_comp = Counter({})
    for run in range(num_run):
        init_node = np.random.choice(range(N),size=10)
        nodelist = [init_node]
        O,dyn_VS = gillespie_optimized([G,nodelist,r0,r1,MAX_TSTEP,NUM_SAMPLES])
        G1  = G.subgraph(dyn_VS)
        cc = sorted(nx.connected_components(G1), key=len, reverse=True)
        dict_comp += Counter([len(item) for item in cc])
        print(f"run {run}\t initial node {init_node}\t GCC {O}")
    # run_parallel(numcpu)
        
    file_path = 'output_OAD.txt'
    # Save the dictionary to a text file
    with open(file_path, 'w') as file:
        for key, value in dict_comp.items():
            file.write(f'{key}: {value}\n')