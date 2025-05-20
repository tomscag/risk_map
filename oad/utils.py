import networkx as nx
import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
import json
from typing import Any

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
    """ A small dataclass to collect results from simulation (cascades). 
    
    Parameters
    ----------
    survived :
        The nodes not affected by the cascading failure.
        It may be a list of node ids (int) or labels (str)
    removed :
        The nodes failed as consequence of the cascade.
        Same type as `survived`
    gcc :
        The fraction of nodes left in the greatest connected component.
    
    """

    survived: list[int | str] = field(default_factory=list)
    removed: list[int | str] = field(default_factory=list)
    gcc: float | None = None


    def to_json(self) -> dict[str, ...]:
        """Convert to JSON serializable types."""
        return {
            "survived": self.survived,
            "removed": self.removed,
            "gcc": self.gcc,
        }

    def sum_survived(self, key: str|None = None) -> int:
        """Compute the sum (or number) of nodes survived."""
        return len(self.survived)

    def set_gcc(self, gcc: float):
        """Set the great connected component size."""
        self.gcc = gcc
        return self
    
    
    def write(self, filepath: str) -> None:
        """Write the Cascade result to json file.
    
        Parameters
        ----------
        filepath : str
            filename output
    
        """
        with Path(filepath).open("w") as fout:
            json.dump(self.to_json(), fout, indent=4)
    

    def __len__(self) -> int:
        """Return the size of the cascade."""
        return len(self.survived) 
