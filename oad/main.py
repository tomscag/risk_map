import networkx as nx
import numpy as np

from oad.utils import Cascade

class OAD:
    """
    
    Run the OAD model on a network [1] 
    using an optimized version of Gillespie algorithm [2]
    [1] https://www.sciencedirect.com/science/article/pii/S0960077924013651#GS1
    [2] https://www.sciencedirect.com/science/article/pii/S0010465517301893
    
    """

    def __init__(self, 
                 G:nx.Graph
                 ) -> None:
        '''
        G: nx.Graph
            The network topology
        '''
        self.num_nodes = G.number_of_nodes()
        self.G = G

        self.net_kmax = max([v for k, v in G.degree()]) # max degree
        self._results = []  

    def run(self,
            lam: float,
            gam: float,
            tmax: float,
            init0:float,
            num_samp: int = 10) -> None:
        """

        Parameters
        ----------
        lam: float
            Local spreading parameter
        gam: float
            Non-local spreading parameter
        tmax: float, optional
            Maximum time to find the solution
        init0: float
            Initial fraction of infected nodes
        num_samp : int
            Number of sample to calculate the average

        Returns
        -------
        float
            The average of the greatest connected component after the 
            dynamics reaches stationarity

        """
        print("Running OAD model ...")
        
        # Rescale the parameters
        lam = lam/self.num_nodes
        gam = gam/(self.num_nodes**2) 

        for _ in range(num_samp):
            
            
            self._results.append(
                Cascade(
                    **self._simulate_gillespie(lam, 
                                               gam, 
                                               tmax, 
                                               init0
                                               )
                    )
                )

            

    def _simulate_gillespie(self,
                            lam: float,
                            gam: float,
                            tmax: float,
                            init0:float) -> dict:
        
        
        nodes_state = {i : 0 for i in self.G.nodes()}  # 0:operational, 1:affected
        list_I = [None for item in range(self.num_nodes)] 
        list_S = [None for item in range(self.num_nodes)]
        num_I  = 0
        num_R  = 0
        tot_deg_I = 0 # Total degree of infected
        mu = 1     # Rescale parameters dividing by mu (see paper [1])

        # Initialization
        # Sort vertices and apply the initial condition
        ver_i_list = [] 
        for ver in np.random.permutation(self.G.nodes()): 
            ver_i_list.append(ver)
            list_I[num_I] = ver
            num_I += 1 # total infected
            nodes_state[ver]= 1 # 
            tot_deg_I += self.G.degree(ver) # 
            if num_I == int(self.num_nodes*init0):
                break
 
        ver_s_list =  list(set(range(self.num_nodes)) - set(ver_i_list))
        num_S = 0
        dyn_NSk = 0
        for ver in ver_s_list:
            list_S[num_S] = ver
            num_S += 1 # total susceptible
            dyn_NSk += self.G.degree(ver)

        # Run the dynamics until convergence
        while True:
            t = 0.0
            dt = 0.0      
            # Calculate the total rate
            M = mu * num_I
            L = lam* tot_deg_I   # This include phantom process
            A = (gam * num_I**2) * num_S # This does not include phantom process (not needed)
            dyn_R = M + L + A

            # Select the time step
            rnd = max(np.random.uniform(),1e-12) # Avoid u = 0
            dt = -np.log(rnd) / dyn_R
            t += dt

            dyn_m = M / dyn_R    # probability of a single event of healing
            dyn_a  = A / dyn_R      # probablity of a field infection
            dyn_lm = dyn_m + L / dyn_R    # probability of a normal infection
            rand   = np.random.uniform()  

            if rand < dyn_m: # If heal (infected in recovered)
                # Select a random occupied vertex and heal.
                pos_inf = np.random.randint(0,num_I) # int random number
                ver = list_I[pos_inf]

                # Then, heal it (put it in Recovered)
                nodes_state[ver] = 2   # 2: disrupted
                tot_deg_I -= self.G.degree(ver)
                num_I -= 1
                num_R += 1 # new
                list_I[pos_inf] = list_I[num_I]
                list_I[num_I] = None 
                
            elif rand >= dyn_m and rand < dyn_lm:  # Infect with the normal mechanism
                # Select the infected vertex i with prob. proportional to k_i
                while True:
                    pos_inf = np.random.randint(0,num_I)
                    ver = list_I[pos_inf]
                    if np.random.uniform() < 1.0*self.G.degree(ver) / (1.0*self.net_kmax):
                        #print('self.G.degree(ver) = ', str(self.G.degree(ver)), ',*net_kmax=', str(net_kmax),  sep= ' ')
                        break
                # Select one of its neighbors
                neighbors_values = [n for n in self.G.neighbors(ver)]
                ver = np.random.choice(neighbors_values)
                if nodes_state[ver] == 0: # if not a phantom process, infect
                    nodes_state[ver] = 1
                    tot_deg_I += self.G.degree(ver)
                    list_I[num_I] = ver    # Add one element to list of I
                    num_I += 1             # Increase by 1 the list
                    num_S -= 1 
                    list_S.remove(ver)  

            else: # Infect with the field mediating mechanism:
                
                # ver = np.random.choice(neighbors_values)
                ver = np.random.choice([item for item in list_S if item is not None])  # [ATTEMPT] here we don't need phantom
                if nodes_state[ver] == 0: # we don't need to test for phantom process, actually
                    nodes_state[ver] = 1

                    tot_deg_I += self.G.degree(ver)
                    list_I[num_I] = ver    # Add one element to list of I
                    num_I += 1             # Increase by 1 the list
                    num_S -= 1 

                    list_S.remove(ver) 

            # Termination.
            if t >= tmax or num_I == 0 or dt < 0:
                # self._results.append(num_S/self.num_nodes)
                survived = [key for key, value in nodes_state.items()
                            if value == 0]
                gcc = len(
                    max(
                        list(
                            nx.connected_components(self.G.subgraph(survived)
                                                    )
                        )
                    )
                )
                return {"survived": survived,
                        "removed": [key for key, value in nodes_state.items() if value == 2],
                        "gcc": gcc}
        
    
    
    
