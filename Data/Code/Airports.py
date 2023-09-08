#!/bin/bash
import networkx as nx

def create_edgelist():
    _fnamepath = "../Raw/Airports/airlines_multiplex_layer-node-node.txt"
    _fnameoutp = "../Processed/Airports/airports_network_extracted.edgelist"

    with open(_fnameoutp,'w+') as fout:
        with open(_fnamepath,'r') as file:
            for item in  file.readlines():
                # print(item)
                node1 = item.split()[1]
                node2 = item.split()[2]
                fout.write(f"{node1} {node2}\n")


import os 
os.chdir("/home/tomsc/Projects/RiskMap/Data/Code")



####################################
####################################
####################################

_fnameoutp = "../Processed/Airports/airports_network_extracted.edgelist"
G = nx.read_edgelist(_fnameoutp,nodetype=int)
print(G.nodes())


G = nx.convert_node_labels_to_integers(G, first_label=0, ordering='default', label_attribute=None)

nx.write_edgelist(G,"./airports_world.edgelist",delimiter=' ',data=False)
