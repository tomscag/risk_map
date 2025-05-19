import networkx as nx

def create_topology(num_nodes: int = 100,
                    p: float = 0.01
                    ) -> nx.Graph:

    return nx.erdos_renyi_graph(n=num_nodes, p=p)

