import numpy as np
import pandas as pd
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



        
def load_data_stressor(event_name=None,name_topology="america"):


    if name_topology == "america" or name_topology == "europe":
        _dir  = f'../Data/Processed/Storms/{name_topology}/storm_'

        data_storm = pd.read_csv(
                    _dir+event_name+'.csv', delimiter=' ',
                    names=["Latitude", "Longitude", "wmo_wind.x", "PDM"],
                    header=0,dtype=float)        

        return data_storm
    elif name_topology == "airports":
        _dir  = f'../Data/Processed/Earthquakes/earthquakes_World_2000-2023_M7.csv'
        data_quakes = pd.read_csv(
                    _dir, delimiter=',', usecols=[1,2,3,4],
                    names=["Latitude", "Longitude", "depth", "magnitude"],
                    header=0,dtype=float) 
        return data_quakes



def load_topology_parameters(name_topology):
    '''
        OUTPUT
            path_edgelist
            path_nodelist
            path_risk
    '''
    if name_topology=="america":
        r0,r1 = (10,0.3)
        return  "../Data/Processed/Topologies/america/powergrid_north_america.el", \
                "../Data/Processed/Topologies/america/america.nodelist", \
                f"./Output_OAD/{name_topology}_r0_{r0}_r1_{r1}_samples_10_maxtime_2000.dat"

    elif name_topology=="europe":
        r0,r1 = (8,0.3)
        return "../Data/Processed/Topologies/europe/powergrid_europe.el", \
                "../Data/Processed/Topologies/europe/europe.nodelist", \
                f"./Output_OAD/{name_topology}_r0_{r0}_r1_{r1}_samples_10_maxtime_2000.dat"
    
    elif name_topology=="airports":
        r0,r1 = (0.3,2)
        return "../Data/Processed/Topologies/airports/airports.edgelist",\
                "../Data/Processed/Topologies/airports/airports.nodelist",\
                f"./Output_OAD/{name_topology}_r0_{r0}_r1_{r1}_samples_10_maxtime_2000.dat"
    
    else:
        print("Topology not recognized \n EXIT")
        return
    

        
def load_topology_geodata(name_topology):
    fpath = f"../Data/Processed/Topologies/{name_topology}/{name_topology}.nodelist"
    return pd.read_csv(fpath,delimiter=" ",index_col=0,names=["label","lng","lat","geoid"],usecols=[0,1,2,3],dtype={"geoid":"string"})



def fragility_model_storm(dist,force,const=1.092e-3,dist0=10, K=1300):
    """ 
    Return the probability of a physical damage
    at distance (dist) and with wind (force)
    """
    if dist > K:
        return 0
    else:
        return  const*((force/65)**(8.02))/(dist+dist0)**2

def fragility_model_earthquake(dist,magnitude,const=1.092e-12,dist0=10, K=1300):
    """ 
    Return the probability of a physical damage
    at distance (dist) and with magnitude
    """
    if dist > K:
        return 0
    else:
        return  const*(10**(1.5*(magnitude-1) )  )/(dist+dist0)**2