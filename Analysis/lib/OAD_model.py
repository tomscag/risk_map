import networkx as nx
import numpy as np



class OAD_model():

    def __init__(self,G,lam,gam,tmax,init0) -> None:
        '''
        G: networkx object
            The graph topology
        lam: float
            Local spreading parameter
        gam: float
            Non-local spreading parameter
        tmax: float, optional
            Maximum time to find the solution
        num_sam: int, optional
            Number of samples to run the dynamics
        init0: float
            Initial fraction of infected nodes
        '''
        self.net_N = G.number_of_nodes()
        self.G = G
        self.lam = lam/self.net_N
        self.gam = gam/(self.net_N**2) 
        self.mu = 1     # To rescale the parameters
        self.tmax = tmax
        self.init0 = init0

        my_degrees    = G.degree()
        degree_values = [v for k, v in my_degrees]
        self.net_kmax      = max(degree_values)


    def RunSimulation(self,num_sam=10):
        '''
        num_sam: int, optional
            Number of samples to run the dynamics
        '''
        self.res_sam = [] # Save results here

        for _ in range(num_sam):

            dyn_sig = { i : 0 for i in self.G.nodes()}  # 0:operational, 1:affected
            dyn_sig = dict.fromkeys(dyn_sig, 0) # {1: 0, 136:0, 992: 0...}
            list_I = [None for item in range(self.net_N)] # [None, None, None...]
            list_S = [None for item in range(self.net_N)]
            num_I  = 0
            num_R  = 0
            tot_deg_I = 0 # Total degree of infected

            # Initialization
            # Sort vertices and apply the initial condition
            ver_i_list = [] 
            for ver in np.random.permutation(self.G.nodes()): # np.random.permutation(10) = array([1, 7, 4, 3, 0, 9, 2, 5, 8, 6])
                ver_i_list.append(ver)
                list_I[num_I] = ver
                num_I += 1 # total infected
                dyn_sig[ver]= 1 # 
                tot_deg_I += self.G.degree(ver) # 
                if num_I == int(self.net_N*self.init0):
                    break

            ver_TOT = range(self.net_N)  
            ver_s_list =  list(set(ver_TOT) - set(ver_i_list)) # set: unordered collection on unique elements
            num_S = 0
            dyn_NSk = 0
            for ver in ver_s_list:
                list_S[num_S] = ver
                num_S += 1 # total susceptible
                dyn_NSk += self.G.degree(ver)


            # Run the dynamics until convergence
            running = True
            while running == True:
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
                    #print('sto guarendo ver= '+ str(ver), ' ,pos_inf=', str(pos_inf))
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
                    self.res_sam.append(num_S/self.net_N)  # Save result here
                    running = False
                # print(num_S)
        
        avg = sum(self.res_sam)/len(self.res_sam)
        return avg