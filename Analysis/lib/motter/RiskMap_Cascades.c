#define NN 13844

#define len(x) (sizeof(x) / sizeof(x[0]))
#define iprint(x) printf("%d\n", x)
#define uep printf("uep\n")
#define hop printf("\n")


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "./PRNG.h"
#include "./RiskMaps.c"


int main(int argc, char *argv[]){
	/*  argv[1] nom de l'arxiu a carregar
	 *  argv[2] numero de la simulacio
	 *  argv[3] seed per al random number  
	 */
	 
    int i, j, k, c, run, flag;
    int intseed;
    int Nl, NN_read, nruns;
    int maxsize;
    int num_initial_removals;
    double alpha;
    double G;

	int *sizes, *failed_nodes, *initial_removals; 
	double *bc, *capacity; 

    int *num_veins1;
	int **veins1 = NULL, **inv_ind = NULL;
	double ave_degree1;

    char finput[100], finput_net1[100], finput_net2[100], fout1[100];

    FILE *fp;

	// PRNG initialization
	intseed = strtol(argv[2],NULL,10);
	srand((unsigned)time(NULL)*intseed);
	double seed = rand();
	init_genrand(seed);
 
	// Malloc arrays
	sizes = malloc((NN) * sizeof(*sizes));
	bc = malloc((NN) * sizeof(*bc));
	capacity = malloc((NN) * sizeof(*capacity));
	failed_nodes = malloc((NN) * sizeof(*failed_nodes));
	initial_removals = malloc((NN) * sizeof(*initial_removals)); //NN is an upper limit

	// Loading params
    sprintf(finput,"%s",argv[1]);
	fp = fopen(finput,"r+");
	rewind(fp);
	
	fscanf(fp, "%d %lf %s %s %s %d", 
			   &Nl, &alpha, finput_net1, finput_net2, fout1, &num_initial_removals);
	who_fails_initially(initial_removals, num_initial_removals, fp);
	fclose(fp);
	printf("%d %f %s %s %s %d\n", Nl, alpha, finput_net1, finput_net2, 
									fout1, num_initial_removals);

	check_inputs(Nl, alpha, num_initial_removals);


	// Simulation loop starts (not needed actually for the Risk Map project)
	nruns = 1;
	for(run=0;run<nruns;run++){
		printf("run %d\n",run);
		
		//Network reading
		ReadEdgeList(&veins1, &num_veins1, &NN_read, finput_net1, finput_net2);

		//inv_ind stands for inverse index. useful for later
		Inverse_Index(veins1, num_veins1, &inv_ind);

		// Computing betweenness and capacity
		BetweennesCentrality_withComponents(veins1, num_veins1, sizes, bc);
		for(i=0;i<NN;i++){
			capacity[i] = (1+alpha)*bc[i];
		}
		
		// Running the cascade
		RemovalAndCascade_Simultaneous_RiskMap(veins1, num_veins1, inv_ind,
					   sizes, bc, capacity, failed_nodes, fout1,
					   initial_removals, num_initial_removals);

		// Freeing
		for(i=0;i<NN;i++){
			free(veins1[i]);
			free(inv_ind[i]);
		}
		free(veins1);
		free(inv_ind);
	}
	
	// Freeing
	free(sizes);
	free(bc);
	free(capacity);
	free(num_veins1);
	free(failed_nodes);
	free(initial_removals);
	
	remove(finput);

	hop;
	return 0;
}

















