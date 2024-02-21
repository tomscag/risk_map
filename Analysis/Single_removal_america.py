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


def analyze(r0,r1):
    """ Simulate the reaction-diffusion over the graph G 
        and save the results in a text file
    """
    # print(f"Parameters r0 {r0:3.0f} r1 {r1:3.2f}")
    size = 20
    dict_comp = Counter({})

    for run in range(num_run):
        nodelist = np.random.choice(range(N),size=size)
        O,dyn_VS = gillespie_optimized([G,nodelist,r0,r1,MAX_TSTEP,NUM_SAMPLES])
        G1  = G.subgraph(dyn_VS)
        cc = sorted(nx.connected_components(G1), key=len, reverse=True)
        dict_comp += Counter([len(item) for item in cc])
        print(f"r0 {r0:3.0f} r1 {r1:3.2f} run {run}\t GCC {O:1.4f}")

    file_path = f"./OAD_[{r0:3.0f},{r1:3.2f}]_nremovals_{size}_nrun_{num_run}.txt"
    # Save the dictionary to a text file
    with open(file_path, 'w') as file:
        for key, value in dict_comp.items():
            file.write(f'{key}: {value}\n')


def run_parallel(par_list,numcpu):
    """ Parallelize the execution of the function analyze """
    pool = mp.Pool(numcpu)
    for i in range(len(par_list)):
        r0,r1 = par_list[i]
        pool.apply_async(analyze, args = (r0, r1, ))
        # analyze(r0,r1)
    pool.close()
    pool.join()



#################################
#################################
#################################


if __name__ == "__main__":
    numcpu = 20
    p_america = 0.00015 
    fpathelist = "./lib/motter/NetworkTopologies/network_america.edgelist"
    G,_  = read_network()
    # filepath_output = make_dirs(f"./Output_OAD/results_cluster_america_motter")

    N = 14990 
    num_run = int(1e2)
    
    MAX_TSTEP = 10^3
    NUM_SAMPLES = 1 # Samples per node
    
    size=20
    r0 = 0.01*1*(N/size)/(p_america) # Circa 5e4
    r1 = 1e0

    r0_list = np.arange(4e4,5e4,int(0.5e3))
    r1_list = np.arange(0,1.1,0.2)
    par_list = [(r0,r1) for r0 in r0_list for r1 in r1_list]
    run_parallel(par_list,numcpu)


            