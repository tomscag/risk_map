import os
import numpy as np


import multiprocessing as mp

########################################
########################################
########################################

def check_arg(node_to_remove):
    for node in node_to_remove:
        if node >= NN or node < 0: #doingg this check inside C script, but just in case...
            raise ValueError("""Manual exiting -- the labels of the nodes I'm removing 
                            should be smaller than %d""" % NN)
        if type(node) != int:
            raise ValueError("""Manual exiting -- the labels of the nodes I'm removing 
                            should be integers, now %f""" % node)



def run_cascades(node_to_remove,NN,alpha):

    script_name_out = "RiskMap_Cascades.c"
    # script_name_out = "RiskMap_Cascades_N%d.c" % NN
    finput_net1 = "NetworkTopologies/network_europe.neighbors"
    finput_net2 = "NetworkTopologies/network_europe.edgelist"
    # alpha = 0.0  # Tolerance

    # os.system("sed '1 i\#define NN %d' %s > %s" % (NN, script_name_in, script_name_out))

    exe_name = './aa_%g_%g.out' % (NN, np.random.randint(0, high = 2147483647, size = None))
    os.system("gcc -g -o %s %s -lm" % (exe_name, script_name_out))


    randseed = np.random.randint(0, high = 2147483647, size = None)
        
    fname='./input_N_'+str(NN)+'_a_'+str(node_to_remove[0])+'.in'

    fout1='Results/N_'+str(NN)+'_a_'+str('{:.3f}'.format(alpha))+'_node_'+str(node_to_remove[0])+'.dat'

    finput=open(fname,'w+')
    finput.write('%d ' % NN)
    finput.write('%f ' % alpha)
    finput.write(finput_net1+' ')
    finput.write(finput_net2+' ')
    finput.write(fout1+' ')
    finput.write('%d\n' % len(node_to_remove))
    for node in node_to_remove:
        finput.write('%d\n' % node)
    finput.close()

    os.system('./'+str(exe_name)+' '+fname+' '+str(randseed)) #run with valgrind when updates are made, to check memory is OK
    # os.system('valgrind --leak-check=yes --track-origins=yes ./'+str(exe_name)+' '+fname+' '+str(randseed))

    os.system("rm %s %s" % (exe_name,fname))



def run_parallel(NN,numcpu,alpha):
    """ Parallelize the execution of the function run_cascades """
    pool = mp.Pool(numcpu)
    for node_to_remove in range(NN):
        node_to_remove = [node_to_remove] #
        pool.apply_async(run_cascades, args = (node_to_remove,NN,alpha  ))
        # run_cascades(node_to_remove,NN)

    pool.close()
    pool.join()



########################################
########################################
########################################

if __name__ == "__main__":

    NN = 13844
    numcpu = 20

    alpha_list = np.arange(0,1,0.1)[1:]
    for alpha in alpha_list:
        run_parallel(NN,numcpu,alpha)