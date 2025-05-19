import networkx as nx
import numpy as np

class OAD:
    """
        Simulate the OAD model on a network
    """

    def __init__(self, 
                 G:nx.Graph, 
                 lam:float, 
                 gam:float, 
                 tmax:float, 
                 init0:float
                 ) -> None:
        '''
        G: nx.Graph
            The graph topology
        lam: float
            Local spreading parameter
        gam: float
            Non-local spreading parameter
        tmax: float, optional
            Maximum time to find the solution
        init0: float
            Initial fraction of infected nodes
        '''
        self.num_nodes = G.number_of_nodes()
        self.G = G
        self.lam = lam/self.num_nodes
        self.gam = gam/(self.num_nodes**2) 
        self.mu = 1     # To rescale the parameters
        self.tmax = tmax
        self.init0 = init0

        my_degrees = G.degree()
        degree_values = [v for k, v in my_degrees]
        self.net_kmax = max(degree_values)

    def run(self, num_sam: int=10) -> float:
        """
        Run OAD model using an optimized Gillespie algorithm
        (ref. https://www.sciencedirect.com/science/article/pii/S0010465517301893)

        Parameters
        ----------
        num_sam : int
            Number of sample to calculate the average

        Returns
        -------
        float
            The average of the greatest connected component after the 
            dynamics reaches stationarity

        """
        print("Running OAD model ...")
        self._results = []  

        for _ in range(num_sam):

            dyn_sig = {i : 0 for i in self.G.nodes()}  # 0:operational, 1:affected
            dyn_sig = dict.fromkeys(dyn_sig, 0) # {1: 0, 136:0, 992: 0...}
            list_I = [None for item in range(self.num_nodes)] 
            list_S = [None for item in range(self.num_nodes)]
            num_I  = 0
            num_R  = 0
            tot_deg_I = 0 # Total degree of infected

            # Initialization
            # Sort vertices and apply the initial condition
            ver_i_list = [] 
            for ver in np.random.permutation(self.G.nodes()): 
                ver_i_list.append(ver)
                list_I[num_I] = ver
                num_I += 1 # total infected
                dyn_sig[ver]= 1 # 
                tot_deg_I += self.G.degree(ver) # 
                if num_I == int(self.num_nodes*self.init0):
                    break

            ver_TOT = range(self.num_nodes)  
            ver_s_list =  list(set(ver_TOT) - set(ver_i_list)) # set: unordered collection on unique elements
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
                M = self.mu * num_I
                L = self.lam * tot_deg_I   # This include phantom process
                A = (self.gam * num_I**2) * num_S # This does not include phantom process (not needed)
                dyn_R = M + L + A

                # Select the time step
                rnd = max(np.random.uniform(),1e-12) # Avoid u = 0
                dt = -np.log(rnd) / dyn_R
                t += dt

                # Probability m to heal
                dyn_m = M / dyn_R # m = M/R, probability of a single event of healing
            
                dyn_a  = A / dyn_R            # probablity of a field mediated infection
                dyn_lm = dyn_m + L / dyn_R 
                rand   = np.random.uniform()  # 

                if rand < dyn_m: # If heal (infected in recovered)
                    # Select a random occupied vertex and heal.
                    pos_inf = np.random.randint(0,num_I) # int random number
                    ver = list_I[pos_inf]

                    # Then, heal it (put it in Recovered)
                    dyn_sig[ver] = 2   # 2: disrupted
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
                    if dyn_sig[ver] == 0: # if not a phantom process, infect
                        dyn_sig[ver] = 1
                        tot_deg_I += self.G.degree(ver)
                        list_I[num_I] = ver    # Add one element to list of I
                        num_I += 1             # Increase by 1 the list
                        num_S -= 1 
                        list_S.remove(ver)  

                else: # Infect with the field mediating mechanism:
                    
                    # ver = np.random.choice(neighbors_values)
                    ver = np.random.choice([item for item in list_S if item is not None])  # [ATTEMPT] here we don't need phantom
                    if dyn_sig[ver] == 0: # we don't need to test for phantom process, actually
                        dyn_sig[ver] = 1

                        tot_deg_I += self.G.degree(ver)
                        list_I[num_I] = ver    # Add one element to list of I
                        num_I += 1             # Increase by 1 the list
                        num_S -= 1 

                        list_S.remove(ver) 

                # Termination.
                if t >= self.tmax or num_I==0 or dt <0:
                    self._results.append(num_S/self.num_nodes)  
                    break
        
        avg = sum(self._results)/len(self._results)
        return avg
    
    
    
