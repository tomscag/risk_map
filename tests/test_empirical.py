# Test OAD model on empirical networks
from oad.main import OAD
import networkx as nx

if __name__ == "__main__":
    print("Loading airports data...")
    folderpath = "../Data/Processed/Topologies/airports"
    nodelist = f"{folderpath}/airports.nodelist"
    edgelist = f"{folderpath}/airports.edgelist"
    G = nx.read_edgelist(edgelist, nodetype=int)
    
    num_samples = 5
    lam, gam = (1e3, 0.01)
    tmax, init0 = (100, 0.05)
    model = OAD(G, init0)
    s = model.run(lam, gam, tmax, num_samp=num_samples)
    print(f"Average susceptibles: {s:1.2f}")
