# 

import numpy as np
import networkx as nx
from Functions import simulate_reaction_diffusion_gillespie
from Functions import simulate_reaction_diffusion_gillespie_fraction

import multiprocessing as mp
import os

# PAckages for profiling
import time 
import cProfile
import pstats

def load_topology(fname,NUM_NODES):
    '''
        if fname=="full"
            Return a fully connected graph

            num_noded: int
                number of nodes
        else
            Return a sampled graph from an existing topology
    '''
    if fname=="random":
        return nx.fast_gnp_random_graph(NUM_NODES, 0.005, seed=None, directed=False)
    elif fname=="full":
        return nx.complete_graph(NUM_NODES), NUM_NODES
    elif fname=="america":
        _fpath   = "../Data/Processed/Topologies/america/powergrid_north_america.el"
        edgelist = nx.read_edgelist(_fpath,nodetype=int)
        if NUM_NODES < 16167:    # If less than the number, sample with configuration model
            return sample_graph_configuration_model(edgelist,NUM_NODES), NUM_NODES
        else:
            NUM_NODES = 16167
            return edgelist, NUM_NODES
    elif fname=="europe":
        _fpath   = "../Data/Processed/Topologies/europe/powergrid_europe.el"
        edgelist = nx.read_edgelist(_fpath,nodetype=int)    
        if NUM_NODES < 1467: # If less than the number, sample with configuration model
            return sample_graph_configuration_model(edgelist,NUM_NODES), NUM_NODES
        else:
            NUM_NODES = 1467
            return edgelist, NUM_NODES
    elif fname=="airports":
        _fpath = "../Data/Processed/Topologies/airports/airports_world.edgelist"
        return nx.read_edgelist(_fpath,nodetype=int), NUM_NODES
    else:
        print("Topology not recognized \n EXIT")
        return


def sample_graph_configuration_model(P,num_nodes):
    """ Generate a subgraph of P of with a smaller number of nodes """
    deg_list         = P.degree()
    nbins            = max(dict(deg_list).values())    # Bins range from 1 to 22
    bins             = list(range(0,nbins+1)) # Increasing array of bin edges, including the rightmost edge
    count, bin_edges = np.histogram(deg_list,bins=bins,range=(1,22))
    prob_deg         = count/sum(count)
    print(f"Number of bins: {nbins}")

    deg_list =  np.random.choice(range(0,nbins), size=num_nodes, p= prob_deg)
    while sum(deg_list)%2 != 0:
        deg_list =  np.random.choice(range(0,nbins), size=num_nodes, p= prob_deg)

    G = nx.configuration_model(deg_list,seed=1234)
    G = nx.Graph(G)     
    G.remove_edges_from(nx.selfloop_edges(G))

    gcc = max([len(item) for item in nx.connected_components(P)])/len(P)
    print(f"Relative size of GCC - US power grid: {gcc:2.2f}")

    gcc = max([len(item) for item in nx.connected_components(G)])/len(G)
    print(f"Relative size of GCC: {gcc:2.2f}")
    return G


def analyze(r0,r1_list,filepath_output):
    """ Simulate the reaction diffusion over the graph G 
        and save the results in a text file
    """
    FNAME_OUTPUT = filepath_output + f"{name_topology}_r0_{r0:.5f}.dat"
    with open(FNAME_OUTPUT,"w+") as file:
        file.write(f"topology: {name_topology}\nnumber of nodes: {NUM_NODES}\nnumber of samples: {NUM_SAMPLES}\nmax time step: {MAX_TSTEP}\ninitial fraction disrupted nodes: {dynp_pINI}\n")

    for r1 in r1_list:
        print(f"analyzing R0: {r0:2.3f}, R1: {r1:2.3f}")
        O = simulate_reaction_diffusion_gillespie_fraction.main([G,r0,r1,MAX_TSTEP,NUM_SAMPLES,dynp_pINI])
        with open(FNAME_OUTPUT,"a+") as file:
            # file.write(str(r1) + "\t" + str(O) + "\n")
            file.write(f"{r1:.6f}"+"\t"+f"{O:.6f}"+"\n") 
        if O < 1e-5: # 
            break   


#########################################################
#########################################################
#########################################################

name_topology   = "europe" # america europe airports random


NUM_NODES    = 5000   # 1000
NUM_SAMPLES  = 75     # 100
MAX_TSTEP    = 2000   # 1000
dynp_pINI    = 0.05   # Fraction of disrupted nodes on the network as initial condition

G, NUM_NODES = load_topology(name_topology,NUM_NODES)

filepath_output = f"./Output_OAD/{name_topology}_nodes_{NUM_NODES}_samples_{NUM_SAMPLES}_maxtime_{MAX_TSTEP}/"

if not os.path.exists(filepath_output):
    os.makedirs(filepath_output)
else:
    print("Directory already present: delete it or change parameters")
    os._exit(0)



r0_list = np.linspace(0,2,50)
r1_list = np.linspace(0,1.1*(1/dynp_pINI),50)

## wip- non equal binning


# np.logspace(0,1/(1-dynp_pINI),base=25,dtype=float)
# # np.interp
# r0_arr = []

# r1_crit = 1
# r0_crit = 1/(1-dynp_pINI) - r1_crit/dynp_pINI
# A = np.linspace(0,r0_crit,25)
# np.concatenate([A, np.logspace(np.log10(r0_crit),np.log10(2),num=25,base=10,dtype=float)])


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
