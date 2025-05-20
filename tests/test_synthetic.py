# Test model on a synthetic network

from oad.main import OAD
from oad.utils import load_network, Cascade


if __name__ == "__main__":
    print("Initializing...")
    network_type = 'barabasi_albert' # erdos_renyi barabasi_albert fully_connected
    num_nodes = 100
    num_samples = 5
    
    if network_type == 'erdos_renyi':
        params = {'p':0.3}
    elif network_type == 'barabasi_albert':
        params = {'m':30}
        
    G = load_network(network_type, num_nodes, params=params)

    lam, gam = (2, 0.01)
    tmax, init0 = (100, 0.05)
    model = OAD(G, init0)
    s = model.run(lam, gam, tmax, num_samp=num_samples)
    print(f"Average susceptibles: {s:1.2f}")
