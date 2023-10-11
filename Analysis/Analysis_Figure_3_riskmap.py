# Analysis risk map


import numpy as np
import networkx as nx
from Functions import simulate_reaction_diffusion_gillespie

import multiprocessing as mp
import os



def load_topology(fname):
    '''
        INPUT:
            fname (str)
                "america" - "europe" - "airports"
    '''

    if fname=="america":
        _fpath   = "../Data/Processed/Topologies/america/powergrid_north_america.el"
        return nx.read_edgelist(_fpath,nodetype=int)

    elif fname=="europe":
        _fpath   = "../Data/Processed/Topologies/europe/powergrid_europe.el"
        return nx.read_edgelist(_fpath,nodetype=int)    

    elif fname=="airports":
        _fpath = "../Data/Processed/Topologies/airports/airports.edgelist"
        return nx.read_edgelist(_fpath,nodetype=int)
    else:
        print("Topology not recognized \n EXIT")
        return
    

def write_results(foutname,O,node):

    with open(foutname,"a+") as file:
        file.write(f"{node}\t{O:.6f}\n")    
    
    


############################################
############################################
############################################


NUM_SAMPLES  = 10     # 100
MAX_TSTEP    = 2000   # 1000

name_topology   = "airports" # america europe airports random


# if not os.path.exists(foutname):
#     os.makedirs(foutname)
# else:
#     print("Directory already present: delete it or change parameters")
#     os._exit(0)


G = load_topology(name_topology)

NUM_NODES = len(G.nodes())


r0    = 0.05
r1    = 1
foutname = f"./Output_OAD/{name_topology}_r0_{r0}_r1_{r1}_samples_{NUM_SAMPLES}_maxtime_{MAX_TSTEP}.dat"
if os.path.exists(foutname):
    print("Existing file")
    os._exit(0)

# for node in range(NUM_NODES):
#     print(f"Analyzing node {node}")
#     O = simulate_reaction_diffusion_gillespie.main([G,[node],r0,r1,MAX_TSTEP,NUM_SAMPLES])
#     write_results(foutname,O,node)


def analyze(G,node,r0,r1,MAX_TSTEP,NUM_SAMPLES,foutname):
    print(f"Analyzing node {node}")
    O = simulate_reaction_diffusion_gillespie.main([G,[node],r0,r1,MAX_TSTEP,NUM_SAMPLES])
    write_results(foutname,O,node)

def apply_async_with_callback():
    """ Parallelize the execution of the function analyze """
    pool = mp.Pool(10)
    for node in range(NUM_NODES):
        pool.apply_async(analyze, args = (G,node,r0,r1,MAX_TSTEP,NUM_SAMPLES,foutname, ))
        # analyze(G,node,r0,r1,MAX_TSTEP,NUM_SAMPLES,foutname)
    pool.close()
    pool.join()

if __name__ == '__main__':
    apply_async_with_callback()