import numpy as np
import networkx as nx
import os



def load_topology(name_topology,num_nodes):
    '''
        if name_topology=="full"
            Return a fully connected graph

            num_nodes: int
                number of nodes
        else
            Return a sampled graph from an existing topology
    '''
    if name_topology=="random":
        return nx.fast_gnp_random_graph(num_nodes, 0.005, seed=None, directed=False)
    elif name_topology=="full":
        return nx.complete_graph(num_nodes), num_nodes
    elif name_topology=="america":
        _fpath   = "../Data/Processed/Topologies/america/powergrid_north_america.el"
        edgelist = nx.read_edgelist(_fpath,nodetype=int)
        if num_nodes < 16167:    # If less than the number, sample with configuration model
            return sample_graph_configuration_model(edgelist,num_nodes), num_nodes
        else:
            num_nodes = 16167
            return edgelist, num_nodes
    elif name_topology=="europe":
        _fpath   = "../Data/Processed/Topologies/europe/powergrid_europe.el"
        edgelist = nx.read_edgelist(_fpath,nodetype=int)    
        if num_nodes < 1467: # If less than the number, sample with configuration model
            return sample_graph_configuration_model(edgelist,num_nodes), num_nodes
        else:
            num_nodes = 1467
            return edgelist, num_nodes
    elif name_topology=="airports":
        _fpath = "../Data/Processed/Topologies/airports/airports.edgelist"
        if num_nodes < 3182:
            return sample_graph_configuration_model(edgelist,num_nodes), num_nodes
        else:
            num_nodes = 3182
            return nx.read_edgelist(_fpath,nodetype=int), num_nodes
    else:
        print("Topology not recognized \n EXIT")
        return


def sample_graph_configuration_model(P,num_nodes):
    """ Sampling from P a subgraph with num_nodes """

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



def make_dirs(filepath_output):
    if not os.path.exists(filepath_output):
        os.makedirs(filepath_output)
    else:
        print("Directory already present: delete it or change parameters")
        os._exit(0)