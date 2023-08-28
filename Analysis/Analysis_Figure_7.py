# 

import numpy as np
import networkx as nx
from Functions import simulate_reaction_diffusion_gillespie 

import multiprocessing as mp

# PAckages for profiling
import time 
import cProfile
import pstats

def load_topology(fname):
    """ Load the edgelist """
    return nx.read_edgelist(fname,nodetype=int)


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


def fully_connected_graph(num_nodes):
    '''
        Return a fully connected graph

            num_noded: int
                number of nodes
    '''
    return nx.complete_graph(num_nodes)


def analyze(r0,r1_list,filepath_output):
    """ Simulate the reaction diffusion over the graph G 
        and save the results in a text file
    """
    FNAME_OUTPUT = filepath_output + f"powergrid_r0_{r0:.5f}.dat"
    with open(FNAME_OUTPUT,"w+") as file:
        file.write(f"number of nodes: {NUM_NODES}\nnumber of samples: {NUM_SAMPLES}\nmax time step: {MAX_TSTEP}\n")

    for r1 in r1_list:
        print(f"analyzing R0: {r0:2.3f}, R1: {r1:2.3f}")
        O = simulate_reaction_diffusion_gillespie.main([G,r0,r1,MAX_TSTEP,NUM_SAMPLES])
        with open(FNAME_OUTPUT,"a+") as file:
            # file.write(str(r1) + "\t" + str(O) + "\n")
            file.write(f"{r1:.6f}"+"\t"+f"{O:.6f}"+"\n")    


#########################################################
#########################################################
#########################################################

filepath_input  = "./Data/Processed/Topologies/Powergrid_Europe/powergrid_europe.el"
filepath_output = "./Analysis/Output_OAD/"

NUM_NODES    = 51   # 500
NUM_SAMPLES  = 10    # 25
MAX_TSTEP    = 10    # 30

P = load_topology(filepath_input)
# G = sample_graph_configuration_model(P,NUM_NODES)
G = fully_connected_graph(NUM_NODES)

r0_list = np.linspace(0,5,10)
r1_list = np.linspace(0,1,10)






def apply_async_with_callback():
    """ Parallelize the execution of the function analyze """
    # pool = mp.Pool(1)
    for i in range(len(r0_list)):
        # pool.apply_async(analyze, args = (r0_list[i], r1_list, filepath_output, ))
        analyze(r0_list[i], r1_list, filepath_output)
    # pool.close()
    # pool.join()

if __name__ == '__main__':
    apply_async_with_callback()


# if __name__ == '__main__':
#     with cProfile.Profile() as profile:
#         apply_async_with_callback()
#     results = pstats.Stats(profile)
#     results.sort_stats(pstats.SortKey.TIME)
#     results.print_stats()
#     results.dump_stats("results.profile")
