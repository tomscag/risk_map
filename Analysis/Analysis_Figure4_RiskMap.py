# Analysis risk map

from lib.gillespie import gillespie_optimized
from lib.misc import (load_topology, make_dirs)
from lib.input import Inputs
import multiprocessing as mp


############################################
############################################
############################################

def write_results(foutname,O,node):

    with open(foutname+"/results.dat","a+") as file:
        file.write(f"{node}\t{O:.6f}\n")    
    
    
def analyze(G,node,r0,r1,MAX_TSTEP,NUM_SAMPLES,foutname):
    print(f"Analyzing node {node}")
    O = gillespie_optimized([G,[node],r0,r1,MAX_TSTEP,NUM_SAMPLES])
    write_results(foutname,O,node)

def apply_async_with_callback(r0,r1,foutname):
    """ Parallelize the execution of the function analyze """
    pool = mp.Pool(10)
    for node in range(NUM_NODES):
        pool.apply_async(analyze, args = (G,node,r0,r1,MAX_TSTEP,NUM_SAMPLES,foutname, ))
        # analyze(G,node,r0,r1,MAX_TSTEP,NUM_SAMPLES,foutname)
    pool.close()
    pool.join()



############################################
############################################
############################################



if __name__ == '__main__':

    
    name_topology = "america" # america europe airports random
    inputs_dct = getattr(Inputs(), name_topology)
    NUM_NODES    = inputs_dct['num_nodes']   # 1000
    NUM_SAMPLES  = inputs_dct['num_samples'] # 100
    MAX_TSTEP    = inputs_dct['max_tstep']   # 1000

    NUM_NODES    = inputs_dct['num_nodes']   # 1000
    NUM_SAMPLES  = 10 # 100
    MAX_TSTEP    = inputs_dct['max_tstep']   # 1000
    r0 = inputs_dct["r0"]
    r1 = inputs_dct["r1"]

    G  = load_topology(name_topology,NUM_NODES)

    
    filepath_output = make_dirs(f"./Output_OAD/single_removal_{name_topology}_r0_{r0}_r1_{r1}_samples_{NUM_SAMPLES}_maxtime_{MAX_TSTEP}")


    print(f"Analyze {name_topology} with single node removal")

    apply_async_with_callback(r0,r1,filepath_output)