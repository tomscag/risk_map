from dataclasses import dataclass
import networkx as nx
import numpy as np

def load_network(
        network_type:str,
        num_nodes:int, 
        params:dict={}):
    """
    
    Parameters
    ----------
    network_type : str
        Can be 'erdos_renyi' 'fully_connected'
    num_nodes : int, optional
        The number of nodes
    params: dict 
        Custom parameters

    Returns
    -------
    nx.Graph

    """

    
    if network_type=="erdos_renyi":
        return nx.fast_gnp_random_graph(n=num_nodes, p=params['p'], seed=None, directed=False)
    elif network_type=="barabasi_albert":
        return nx.barabasi_albert_graph(n=num_nodes, m=params['m'])
    
    elif network_type=="fully_connected":
        return nx.complete_graph(num_nodes)
    
    else:
        raise Exception("Topology not recognized \n")



@dataclass
class Cascade:
    """A small dataclass to collect results from cascades. """

    def to_json(self) -> dict[str, ...]:
        """Convert to JSON serializable types."""
        return {
            "time": np.datetime_as_string(self.time).tolist() if self.time is not None else None,
            "failing_0": self.failing_0,
            "failing_1": self.failing_1,
            "gcc": self.gcc,
        }


    def __len__(self) -> int:
        """Return the size of the cascade."""
        return len(self.failing_0) + len(self.failing_1)
