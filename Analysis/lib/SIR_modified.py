import networkx as nx
import numpy as np

# TO-DO
# Obtain the same results in a controlled case using old and new code
# Generalize to the modified model (alpha parameter)
# Launch the simulation in the cloud     


class SIR_Simulation():

    def __init__(self, A, lam=0.1, gam=0.001, alp=0.0, i0=0.1) -> None:

        # Network setup.
        if type(A)==np.ndarray:
            self.A = A
        elif type(A)==nx.classes.graph.Graph:
            self.A = nx.adjacency_matrix(A).todense()
        else:
            raise BaseException("Input an adjacency matrix or networkx object only.")

        # Model Parameters.
        self.N = self.A.shape[0]
        self.lam = lam # Rate is prop to (1/N)
        self.gam = gam
        self.alp = alp

        # Time-keeping.
        self.t = 0
        self.times = [0]
        
        # Node numbers.
        self.I = [int(i0*self.N)]
        self.S = [self.N-self.I[0]]
        self.R = [0]

        # Node states.
        self.X = np.array([1]*self.S[0] +[2]*self.I[0]).reshape((self.N,1))
        np.random.shuffle(self.X)
         
        # Initial propensity setup.        
        self.UpdatePropensityALL() 



    def UpdatePropensityALL(self, n_nodes=None):
        self.IPL = (self.lam*self.A.dot(self.X==2))*(self.X==1)
        self.IPG =  self.gam*(self.X==2)
        return None

    def UpdatePropensityFANCY(self,n_nodes):    # Use fancy indexing        
        self.IPL[n_nodes] = (self.lam*self.A[n_nodes].dot(self.X==2))*(self.X[n_nodes]==1)
        self.IPG =  self.gam*(self.X==2)
        return None
    
    def UpdatePropensityTAKE(self,n_nodes):   
        self.IPL[n_nodes] = (self.lam*np.take(self.A,n_nodes, axis=0).dot(self.X==2))*((np.take(self.X,n_nodes, axis=0)==1))
        self.IPG =  self.gam*(self.X==2)
        return None  

    def UpdatePropensityLOOP(self,n_nodes):        
        for node in n_nodes:
            self.IPL[node] = (self.lam*self.A[node].dot(self.X==2))*(self.X[node]==1)
            self.IPG[node] =  self.gam*(self.X[node]==2)
        return None




    def RunIteration(self):
        
        # Termination.
        if self.I[-1] == 0:
            self.S = np.array(self.S, dtype=float)
            self.I = np.array(self.I, dtype=float)
            return False

        # 1. Generate random numbers r1,r2 uniformly distributed in (0,1)
        r1 = np.random.rand()
        r2 = np.random.rand()
        
        # 2. Calculate alpha.
        cumsum = np.concatenate((self.IPL.cumsum(),self.IPL.cumsum()[-1]+self.IPG.cumsum()),dtype=float)
        self.alpha = cumsum[-1]

        # 3. Compute the time until the next reaction takes place.
        tau = (1.0/self.alpha)*np.log(float(1.0/r1))
        self.t += tau
        self.times.append(self.t)

        # 4. Compute which reaction takes place.
        index = np.searchsorted(cumsum,r2*self.alpha)

        # 5. Update node states. 
        if (index//self.N == 0):    # Infection
            self.X[index%self.N] = 2
            self.S.append(self.S[-1] - 1)
            self.I.append(self.I[-1] + 1)
            self.R.append(self.R[-1])
        elif (index//self.N == 1):   # Recover
            self.X[index%self.N] = 3
            self.R.append(self.R[-1]+1)
            self.I.append(self.I[-1]-1)
            self.S.append(self.S[-1])
        else:
            raise Exception("Node state to update is wrong")
            
        # 6. Update propensities (FANCY or TAKE is the most efficient)
        n1 = np.nonzero(self.A[index%self.N])[0]
        n1 = np.append(n1,index%self.N)
        
        # self.UpdatePropensityALL()
        self.UpdatePropensityFANCY(n1)
        # self.UpdatePropensityTAKE(n1)
        # self.UpdatePropensityLOOP(n1)
        return True

    def RunToConvergence(self):
        running = True
        while running:
            running = self.RunIteration()
        return None




    
    def ParametricSweep(self,G,
                        lam_list = np.arange(0.01,0.5,0.05),
                        gam_list = np.arange(0.01,0.5,0.05),
                        filename = 'test.npz'
                        ):
        
        results = np.empty(shape=(lam_list.shape[0],gam_list.shape[0]))
        for idl, lam in enumerate(lam_list):
            for idg,gam in enumerate(gam_list):
                print(f"Analyzing {lam:.3f} and {gam:.3f}")
                model = SIR_Simulation(G, lam=lam, gam=gam, alp=0.0, i0=0.01)
                model.RunToConvergence()
                results[idl,idg] = model.S[-1]
        
        np.savez(file=filename,results=results,lam_list=lam_list,gam_list=gam_list)
        return None


