#!/usr/bin/env python

#import sys
import numpy as np
from math import log


def main(argv):


    '''
        Takes a network as input
        Takes two parameters R0 R1 
        For each node runs the dynamic (ODA model)
        Takes a number of nodes (optional parameter, default, all network)
        Random removal, number of stories (at least 25)
        Returns the average of the asymptotic solution
    '''

    # READING PARAMETERS
    nw         = argv[0]  # Networkx object

    dynp_r0    = argv[1]   #'Value of infection rate lambda: es:0.002
    dynp_r1    = argv[2]   # Field infection rate alpha  es:0.1

    dynp_tmax  = argv[3]  #'Maximum time steps (it stops if the absorbing state is reached): ' es 1000000
    num_sam    = argv[4]  # number of samples es: 10
    # num_nodes  = argv[5]  # number of nodes (defaulf: net_N)



    dynp_mu = 1         # healing rate default: 1

    # import networkx as nx
    net_N         = nw.number_of_nodes()
    my_degrees    = nw.degree()
    degree_values = [v for k, v in my_degrees]
    net_kmax      = max(degree_values)
    
    res_sam = {i: [] for i in nw.nodes()}   # Save results for each node as a list
                                                
    dyn_sig = { i : 0 for i in nw.nodes()}   # 0:operational, 1:affected, 2:disrupted
    

    for init_node in nw.nodes():
        # print(f"\n Node # {init_node}")

        for sam in range(1,num_sam+1):

            # Initialize network
            dyn_sig = dict.fromkeys(dyn_sig, 0) # {1: 0, 136:0, 992: 0...}
            dyn_VI = [None]*net_N # [None, None, None...]
            dyn_VS = [None]*net_N
            dyn_NI  = 0
            dyn_NIk = 0
            dyn_NR  = 0 #new

            # Sort vertices and apply the initial condition
            ver_i_list = []
            ver_i_list.append(init_node)
            dyn_VI[dyn_NI] = init_node
            dyn_NI += 1 # total affected
            dyn_sig[init_node] = 1    # 1: affected, 0: operational
            dyn_NIk += nw.degree(init_node) # somma dei degree degli infetti ? N_IS ?        

            # populate NS and dyn_VS (initial  condition)
            ver_TOT = range(net_N)  
            ver_s_list =  list(set(ver_TOT) - set(ver_i_list)) # set: unordered collection on unique elements
            dyn_NS = 0
            dyn_NSk = 0
            for ver in ver_s_list:
                dyn_VS[dyn_NS] = ver
                dyn_NS += 1 # total susceptible
                dyn_NSk += nw.degree(ver)

            # Running dynamics
            dyn_t  = 0
            dyn_dt = 0.0

            while dyn_t <= dynp_tmax and dyn_NI > 0 and dyn_dt>=0:
                # Modified SIR-OGA ALGORITHM  
            
                # Calculate the total rate
                M = dynp_mu * dyn_NI 
                L = dynp_r0 * dyn_NIk
                A = dynp_r1 * (dyn_NI/net_N)*( dyn_NIk )   #  [TS] Rate due to the field interaction
                dyn_R = M + L + A

    
                # That is total healing rate of the network + total infection rate of the network 
                # that is sum(q=1:Z)ni_q, con Z = NI + NIS spontaneous processes with rates 
                # respectively mu and lambda
            
                # Select the time step
                rnd = max(np.random.uniform(),1e-12) # Avoid u = 0
                dyn_dt = -log(rnd) / dyn_R
                #print('M=' + str(M) + ',L='+ str(L) + ',dt= ' +  str(dyn_dt))
            
                # Update the time
                dyn_t += dyn_dt
            
                # Probability m to heal
                dyn_m = M / dyn_R # m = M/R, probability of a single event of healing
                
                dyn_a  = A / dyn_R            # probablity of a field mediated infection
                dyn_lm = dyn_m + L / dyn_R    # cumulative probability L+M, for rejection method
                rand   = np.random.uniform()  # 

                if rand < dyn_m: # If heal (infected in recovered)
                    # Select a random occupied vertex and heal.
                    pos_inf = np.random.randint(0,dyn_NI) # int random number
                    ver = dyn_VI[pos_inf]

                    # Then, heal it (put it in Recovered)
                    #print('sto guarendo ver= '+ str(ver), ' ,pos_inf=', str(pos_inf))
                    dyn_sig[ver] = 2   # 2: disrupted
                    dyn_NIk -= nw.degree(ver)
                    dyn_NI -= 1
                    dyn_NR += 1 # new
                    dyn_VI[pos_inf] = dyn_VI[dyn_NI]
                    dyn_VI[dyn_NI] = None 

                elif rand >= dyn_m and rand < dyn_lm: # If not, try to infect: w = 1 - m 
                    # Select the infected vertex i with prob. proportional to k_i
                    while True:
                        #print('se scelgo che vertice infetto')

                        pos_inf = np.random.randint(0,dyn_NI)
                        ver = dyn_VI[pos_inf]
                        #print(ver)
                        #print('sto cercando di infettare con ver= '+ str(ver))
                        if np.random.uniform() < 1.0*nw.degree(ver) / (1.0*net_kmax):
                            #print('nw.degree(ver) = ', str(nw.degree(ver)), ',*net_kmax=', str(net_kmax),  sep= ' ')
                            break
                
                    # Select one of its neighbors
                    neighbors_values = [n for n in nw.neighbors(ver)]
                    ver = np.random.choice(neighbors_values)
                    if dyn_sig[ver] == 0: # if not a phantom process, infect
                        dyn_sig[ver] = 1
                        #dyn_NSk -= nw.degree(ver)

                        dyn_NIk += nw.degree(ver)
                        dyn_VI[dyn_NI] = ver    # Add one element to list of I
                        dyn_NI += 1             # Increase by 1 the list
                        dyn_NS -= 1 

                        dyn_VS.remove(ver)
                        
                    # if a absorbing state is reached, exit

                else:    # [TS] else, infect with the field mediating mechanism

                    # Select the infected vertex i with prob. proportional to k_i
                    while True:
                        pos_inf = np.random.randint(0,dyn_NI)
                        ver = dyn_VI[pos_inf]
                        if np.random.uniform() < 1.0*nw.degree(ver) / (1.0*net_kmax):
                            #print('nw.degree(ver) = ', str(nw.degree(ver)), ',*net_kmax=', str(net_kmax),  sep= ' ')
                            break

                    # [TS] In this case, we select uniformly from all susceptible nodes (virtual links) # neighbors_values = [x for x in dyn_VS if x is not None]
                    #  We select all to include phantom processes

                    # CHOICE 1: We sample new infected from all the nodes (only susceptibles, excluding phantom)
                    #           This allows "teleportation" of infection towars other regions
                    neighbors_values = [x for x in ver_TOT if x != ver]

                    # CHOICE 2: We sample only the neighbors of infected (standard SIR, just with an enhanced rate)
                    # neighbors_values = [n for n in nw.neighbors(ver)]

                    ver = np.random.choice(neighbors_values)
                    if dyn_sig[ver] == 0: # if not a phantom process, infect
                        dyn_sig[ver] = 1
                        #dyn_NSk -= nw.degree(ver)

                        dyn_NIk += nw.degree(ver)
                        dyn_VI[dyn_NI] = ver    # Add one element to list of I
                        dyn_NI += 1             # Increase by 1 the list
                        dyn_NS -= 1 

                        dyn_VS.remove(ver)                    
                    # if a absorbing state is reached, exit
            
            # Compute fraction of operational nodes in the absorbing state
            # nf  = nw.subgraph(dyn_VS)
            # gcc = max([len(item) for item in nx.connected_components(nf)])/net_N
            print(dyn_NS)
            res_sam[init_node].append(dyn_NS/net_N)  # Save result here

    
    avg = 0
    for key,value in res_sam.items():
        avg += sum(value)/len(value)
    avg = avg/net_N

    return avg
    





"""
RUN MAIN
"""
if __name__== "__main__":
    main(sys.argv[1:])
  