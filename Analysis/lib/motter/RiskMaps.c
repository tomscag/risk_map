
void Check_Reprocity(int **veins, int num_veins[]){

	int i, j, kk, flag;
	int multiedges;

	multiedges = 0;	// Multiedges no
	//~ multiedges = 1; // Multiedges sí

	for(i=0; i<NN; i++){
		// Num de veins acceptable
		if((num_veins[i] < -1) || (num_veins[i] > NN)){
			printf("Warning-1 %d %d!\n", i, num_veins[i]);
			printf("\nExiting\n"); exit(EXIT_FAILURE);
		}
		
		for(j=0; j<num_veins[i]; j++){
			flag = 0;
			//Self-loops
			if(veins[i][j] == i){
				printf("Warning0!\n");
				printf("\nExiting\n"); exit(EXIT_FAILURE);
			}
			
			// Check Multiedges
			if(multiedges == 0){
				for(kk=(j+1); kk<num_veins[i]; kk++){
					if(veins[i][kk] == veins[i][j]){ 
						printf("Warning1!\n");
						printf("\nExiting\n"); exit(EXIT_FAILURE);
					}
				}
			}
			//Simetria
			for(kk=0;kk<num_veins[veins[i][j]];kk++){
				if(veins[veins[i][j]][kk] == i){
					flag = 1;
				}
			}
			if(flag != 1){
				printf("Warning2!\n");
				printf("\nExiting\n"); exit(EXIT_FAILURE);
			}
			//Límits
			if ( (veins[i][j] >= NN) || (veins[i][j] < 0) ) {
				printf("Error al construir la network... CR\n");
				printf("\nExiting\n"); exit(EXIT_FAILURE);
			}
		}
	}

	return;
}

void ReadEdgeList(int ***veins, int **num_veins, int *NN_read, 
				  char *filename_neigh, char *filename_network){
	/* Function that inputs the network connections into the matrices
	 * from two files: filename_neigh and filename_network. The former is
	 * a two column files where first column is the node and the second
	 * columns is the number of neighbors of that node. The latter file
	 * is the edgelist: two columns, each column indicating the nodes
	 * at the extreme of the links. 
	 * The other arguments are the neighbor matrix and the number of
	 * neighbors of each node.
	 * 
	 * Assumtpions:
	 * NETWORK IS ASSUMED INDIRECTED
	 * filename_network ONLY CONTAINS THE LINKS ONCE (ie, if node 0 is
	 * connected to node 4, in the file will be written in one row "0 4"
	 * or "4 0", but only one of these in the whole file.
	 * Not self loops allowed. If any, function exits.
	 * 
	 * BE WARE:
	 * Avoid having \n lines at the end of the file. 
	 * Avoid having files with a number of columns different than 2.
	 * The files need to be in the same directory than the main program
	 * 
	 * SUPER IMPORTANT: The nodes need to be integers starting at 0 and
	 * growing continuously. I added a checker, called checker.
	 * 
	 * On how to read files: 
	 * https://stackoverflow.com/questions/23185622/fscanf-reads-the-last-integer-twice
	 * 
	 * To prepare the files in Python, I use:

		def WriteForC(OriginalGrraph):
			H_int = nx.convert_node_labels_to_integers(OriginalGrraph)

			fh = open("network.edgelist",'wb')
			nx.write_edgelist(H_int, fh, data=False)
			fh.close()

			fh = open("network.neighbors",'w+')
			for node in list(H_int.nodes()):
				fh.write("%d %d\n" % (node, len(list(H_int.neighbors(node)))))
			fh.close()
	 * 
	 * If I use the above code, then the filename inputs of this function
	 * should be filename_network = network.edgelist and 
	 * filename_neigh = network.neighbors.
	 * 
	 */
	 
	int i, j, sl, checker, c, c2, lc, lc2;
	int r;
	int node, node1, node2, nveinsnode;
	int *nveins_count;

	FILE *file;

    char *remove_p = malloc(256);
	char remove_string[256];

    char *command_p = malloc(256);
	char command[256];

    char *fname_p = malloc(256);
	char fname[256];

	r = genrand_int32();
	
	strcpy(remove_p, "rm temp_%d.txt");
	sprintf(remove_string, remove_p, r);	
	
	// Trobo el numero total de nodes i el numero de linies dels fitxers//
	strcpy(command_p, "awk 'BEGIN{max=0}{{for(i=1;i<=NF;i++) if($i>max) max=$i;}}; END { print max + 1;}' %s > temp_%d.txt");
	sprintf(command, command_p, filename_neigh, r);
	//~ printf("\n%s\n", command);
	system(command);
	
	strcpy(fname_p, "temp_%d.txt");
	sprintf(fname, fname_p, r);	
	
	file = fopen(fname,"r");
	if(file == NULL){
		perror("Error in opening file");
		printf("\nExiting\n"); exit(EXIT_FAILURE);
	}
	fscanf(file,"%d",&(*NN_read));
	printf("Size of the Network %d\n", (*NN_read));
	fclose(file);

	system(remove_string);
	if((*NN_read)!=NN){
		printf("Global definition of network size needs to be changed ");
		printf("%d %d\nExiting", (*NN_read), NN);
		printf("\nExiting\n"); exit(EXIT_FAILURE);
	}
	
	
	strcpy(command_p, "cat %s | wc -l > temp_%d.txt");
	sprintf(command, command_p, filename_neigh, r);
	//~ printf("\n%s\n", command);
	system(command);
	file = fopen(fname,"r");
	if(file == NULL){
		perror("Error in opening file");
		printf("\nExiting\n"); exit(EXIT_FAILURE);
	}
	fscanf(file,"%d",&lc);
	fclose(file);
	system(remove_string);


	strcpy(command_p, "cat %s | wc -l > temp_%d.txt");
	sprintf(command, command_p, filename_network, r);
//	printf("\n%s\n", command);
	system(command);
	file = fopen(fname,"r");
	if(file == NULL){
		perror("Error in opening file");
		printf("\nExiting\n"); exit(EXIT_FAILURE);
	}
	fscanf(file,"%d",&lc2);
	fclose(file);
	system(remove_string);
	
	

	// Faig les allocations dels veins//
	file = fopen(filename_neigh, "r");
	if(file == NULL){
		perror("Error in opening file");
		printf("\nExiting\n"); exit(EXIT_FAILURE);
	}
	
	(*veins) = malloc((*NN_read) * sizeof(**veins));
	(*num_veins) = malloc((*NN_read) * sizeof(*(*num_veins)));
	nveins_count = malloc((*NN_read) * sizeof(*nveins_count));

	checker = -1;
	c = 0;
	while(!feof(file)){
		if(fscanf(file, "%d %d", &node, &nveinsnode) == 2){ //https://stackoverflow.com/questions/23185622/fscanf-reads-the-last-integer-twice
			//~ printf("%d %d\n", node, nveinsnode);
			if(node == checker + 1){
				(*num_veins)[node] = nveinsnode;
				(*veins)[node] = malloc(((*num_veins)[node] + 1) * sizeof(*(*veins)[node])); 	
				nveins_count[node] = 0;
				checker++;
			}else{
				printf("%d Nodes are not ordered, exiting\n", checker);
				printf("\nExiting\n"); exit(EXIT_FAILURE);
			}
		//	printf("\t%d %d\n", node, (*num_veins)[node]);
		}else{
			if(c<lc){
				printf("\n%d\n", c);
				printf("1 - Wrong data structure in the files, exiting\n");
				printf("\nExiting\n"); exit(EXIT_FAILURE);
			}
		}
		c++;
	}
	fclose(file);

	// Carrego la xarxa //
	file = fopen(filename_network, "r");
	if(file == NULL){
		perror("Error in opening file");
		printf("\nExiting\n"); exit(EXIT_FAILURE);
	}
	
	sl = 0;
	c2 = 0;
	while(!feof(file)){
		if(fscanf(file, "%d %d", &node1, &node2) == 2){
//			printf("%d %d\n", node1, node2);
			(*veins)[node1][nveins_count[node1]] = node2;
			(*veins)[node2][nveins_count[node2]] = node1;
			nveins_count[node1]++;
			nveins_count[node2]++;
		}else{
			if(c2<lc2){
				printf("\n%d\n", c2);
				printf("2 - Wrong data structure in the files, exiting\n");
				printf("\nExiting\n"); exit(EXIT_FAILURE);
			}
		}
		if(node1==node2){
			sl++;
		}
		c2++;
	}
	fclose(file);

	if(sl!=0){
		printf("%d self-loops found, exiting\n", sl);
		printf("\nExiting\n"); exit(EXIT_FAILURE);
	}
	
	
	// Security check...//
	for(i=0;i<(*NN_read);i++){
		if(nveins_count[i] != (*num_veins)[i]){
			printf("\n%d %d %d\n", i, nveins_count[i], (*num_veins)[i]);
			printf("Wrong assigning of nodes, exiting\n");
			printf("\nExiting\n"); exit(EXIT_FAILURE);
		}
	}

	//~ strcpy(command_p, "rm %s %s");
	//~ sprintf(command, command_p, filename_neigh, filename_network);
	//~ system(command);
	
	free(nveins_count);
	free(command_p);
	free(fname_p);
	free(remove_p);

	return;
}

void ReconstructionOfTheNetwork(int **veins, int num_veins[],
								int **veins_out, int num_veins_out[],
								int current_network[]){
	/* Algorithm that gives A single connected component of a global
	 * disconnected network (veins), in the same format (veins_out and
	 * num_veins_out). Since the original network veins can have more 
	 * than disconnected component, the one that we are reconstructing 
	 * is indicated by the array current_network. It is an n-dimensional
	 * array of ints, 0 or 1, indicating the nodes we have to consider
	 * to reconstruct the network (the ones inside the component we are
	 * interested in). That means that before calling this functions we
	 * need to have the network already partioned, with for example an
	 * algorithm of breadth-first search. The aim of this function is 
	 * to help finding the betweenness centrality of a disconnected 
	 * network.
	 * 
	 * The memory allocation for the output network of neighbors needs
	 * to be done outside the function.
	 * 
	 * Created to be used with BetweennesCentrality_InsideComponent() 
	 * and BetweennesCentrality_withComponents().
	 * 
	 */ 
	int i, j, k;
	int *write_index = NULL;
	
	write_index = malloc(NN * sizeof((*write_index)));
	
//	num_veins_out[0]=666;
//	veins_out[0][0]=555;
	for(i=0;i<NN;i++){
		num_veins_out[i]=-1; //Vol dir que node i no participa en la component
		write_index[i]=0;
	}

	for(i=0;i<NN;i++){
		if(current_network[i]==1){
			for(j=0;j<num_veins[i];j++){
				veins_out[i][write_index[i]] = veins[i][j];
				write_index[i]++;
			}
			num_veins_out[i] = write_index[i];
		}
	}

	free(write_index);

	return;
}

void BetweennesCentrality_InsideComponent(int **veins_comp, 
										  int num_veins_comp[], 
										  double bc[], 
										  int current_network[]){

	/* Algorithm that gives the betwennes centrality, that is, the number
	 * of shortest paths that cross each node, when the shortest paths
	 * between all pair of nodes are considered. I have to give to this
	 * function the network of the connected component and the number
	 * of neighbor of its nodes. Therefore, if the whole network forms
	 * a single connected component, this function can be also used to
	 * compute its BC (although I believe that it is a bit slower than
	 * BetweennessCentrality() [needs to be checked].
	 * 
	 * We pass an array bc[], as an argument and update it with this
	 * score. We do not normalize, so the bc is > 1, and == 1 when we
	 * have a connected component of only one node. Another input is
	 * the array current_network[i], which gives me if the node i is 
	 * inside the connected component. This is used to ease the code.
	 * 
	 * This functions has been created to be used together with functions
	 *  ReconstructionOfTheNetwork() and 
	 *  BetweennesCentrality_withComponents().
	 * 
	 * I have modified a bit the algorthm given by Newman in the 
	 * section 10.3.6 of his book. First I compute the shortest path tree
	 * (SPT), which is a directed network of shortest paths from any 
	 * node towards the node indicated as source. The way I assign the
	 * betweenness centrality is the following: I start in the fartest 
	 * node from the source (dist[-1]) and increase its bc one unity. 
	 * The node(s) that it has (they have) 'upward' in the shortest 
	 * path tree (one jump closer to the source) also gains (gain) some 
	 * bc. If this (these) upward node(s) has (have) not been considered
	 * before, it (they) gains (gain) the downward node's bc plus 1 of 
	 * its (their) own. If the node(s) has (have) other downwards 
	 * neighbors, these contibute only by their bc. So the +1 
	 * contribution of an upward node is given by its first appearance 
	 * in the SPT.
	 * 
	 * IMPORTANT REMARKS:
	 * 1 - The algorithm allows for multiple shortest paths of same 
	 * distance between two nodes. In these case, the bc is divided 
	 * proportionally.
	 * 2 - The nodes themselves count in the shortest path score, i.e., 
	 * the path between the neighbours node 2 to node 9, contributes in
	 * a bc of 1 for node 2 and node 9.
	 * 3 - Corollary of the latter point: We should not normalize by
	 * connected components since it is a bit buggy: since we can have
	 * connected components of single nodes, they are not normalizable
	 * (it would diverge using the standard definition). 
	 * 4 - We assume an unweighte network.
	 * 
	 * 
	 */
	
	int i, j, c, source, upward_neigh, dmax, read, write, node;
	int nodes_in_comp;
	int *queue = NULL, *write_n = NULL, *dist = NULL, *weigths = NULL, *selected = NULL;
	double *bc_partial = NULL;

    int *num_veins_SPT = NULL;
	int **veinsSPT = NULL;

	queue = malloc(NN * sizeof((*queue)));
	write_n = malloc(NN * sizeof((*write_n)));
	dist = malloc(NN * sizeof((*dist)));
	weigths = malloc(NN * sizeof((*weigths)));
	selected = malloc(NN * sizeof((*selected)));
	bc_partial = malloc(NN * sizeof((*bc_partial)));
	num_veins_SPT = malloc(NN * sizeof((*num_veins_SPT)));

	nodes_in_comp = 0;
	for(i=0;i<NN;i++){
		num_veins_SPT[i] = num_veins_comp[i];
		 if(current_network[i]==1){
			nodes_in_comp++;
		 }
	}

	veinsSPT = malloc((NN) * sizeof(*veinsSPT));
	for(i=0;i<NN;i++){
		if(num_veins_SPT[i]!=-1){
			veinsSPT[i] = malloc(num_veins_SPT[i]*sizeof(*veinsSPT[i]));
		}
	}

	//Loop sobre tots els nodes
	for(source=0;source<NN;source++){
		if(current_network[source]==1){
			if(num_veins_comp[source]==0){	//If it is a single node as connected component, I directly gives the bc
				bc[source]++;	
			} else {						//If its not, I have to construct the SPT and compute the BC there
				dmax = 0;
				queue[0] = source;
				dist[source] = 0;			//distance from source to source
				weigths[source] = 1;

				for(i=0;i<NN;i++){			//All other distance are unknown (-1)
					bc_partial[i] = 0;
					write_n[i] = 0; 		//it gives me the next empty slot in the array veins[i] to write an i's neighbor
					selected[i] = 0;		//auxilar variable to know the contribution in the bc due to neighbors
					if(i!=source){
						dist[i] = -1; 		//dist from source to node i
						weigths[i] = 0; 	//auxilar variable to count the paths in a given node
					}
				}
				

				read = 0;
				write = 1;

				veinsSPT[source][write_n[source]] = source;
				write_n[source]++;
				
				// Construction of the shortest path tree
				for(read=0;read<NN;read++){
					if(read<nodes_in_comp){
						for(j=0;j<num_veins_comp[queue[read]];j++){							//Miro els veins del node que estic llegint
							node = veins_comp[queue[read]][j];
							if(dist[node] == dist[queue[read]] + 1){						//Aquest if per quan tingui multiples SP
								veinsSPT[node][write_n[node]] = queue[read];
								write_n[node]++;
								weigths[node] += weigths[queue[read]];
							}
							if(dist[node] == -1){											//Si el node seleccionat apareix per primer pic
								dist[node] = dist[queue[read]] + 1;
								weigths[node] = weigths[queue[read]];
								if(dist[node] > dmax){
									dmax = dist[node];
								}
								queue[write] = node;
								veinsSPT[node][write_n[node]] = queue[read];
								write_n[node]++;
								write++;
							}
						}
					}
				}

				//Assignation of the BC values
				for(i=nodes_in_comp-1;i>-1;i--){
					node = queue[i];
					if(bc_partial[node]==0){
						bc_partial[node]++; //initial contribution because it is a leaf
					}

					//I do not go through all the path, but just the above node
					for(j=0;j<write_n[node];j++){
						upward_neigh = veinsSPT[node][j];
						if(node!=source){ 									//No contribution of source because it is the end point
							if(selected[upward_neigh]==0){					//Si ha estat seleccionat per primer pic, li sumo +1 de contribucio propia al SP
								selected[upward_neigh]=1;
								bc_partial[upward_neigh] += (1 + bc_partial[node]*weigths[upward_neigh]/(float)weigths[node]);
							}else{											//Si ja ha estat seleccionat abans, no li sumo el +1
								bc_partial[upward_neigh] += bc_partial[node]*weigths[upward_neigh]/(float)weigths[node];
							}
						}
					}
				}
			
				//Actualitzacio del bc total
				for(i=0;i<NN;i++){
					bc[i] += bc_partial[i];
				}
			}
		}
	}

	for(i=0;i<NN;i++){
		if(num_veins_SPT[i]!=-1){
			free(veinsSPT[i]);
		}
	}
	free(veinsSPT);

	free(queue);
	free(write_n);
	free(dist);
	free(weigths);
	free(selected);
	free(bc_partial);
	free(num_veins_SPT);
	
	return;
}



int BetweennesCentrality_Motter(int **veins, int num_veins[], 
								int sizes[], double bc[]){
	
	/* Algorithm that finds the betweenness centrality in a network
	 * that can be formed of different disconnected components, of 
	 * different size (even of 1 node), 
	 * as well as the number of connected components and their
	 * size in a network, by means of a breadth-first search.
	 * The network can be disconnected in two ways: either two components
	 * do not share links, or I have deleted nodes, breaking the network
	 * in pieces. This will be the case for the Motter model, in which
	 * I delete the node with highest betweenness centrality, 
	 * substituting all its appearances in the network structure by a -1.
	 * What this function does, at odds with 
	 * BetweennesCentrality_withComponents() is to create a new network
	 * where the -1 are not present as neighbors (with the purpose of 
	 * computing the bc, not to modify the original network).
	 * 
	 * Each node can be in three colors: white (non-discovered), 
	 * gray (discovered but still with white neighbors) and 
	 * black (discovered but no white neighbors). 
	 * 
	 * This function has been created to be used together with 
	 * ReconstructionOfTheNetwork() and 
	 * BetweennesCentrality_InsideComponent()
	 * 
	 * With litle work could give the network structure of the 
	 * disconnected components.
	 * 
	 /* Ref: Introduction to Algorithms, T.H. Cormen, C.E. Leiserson, R.L. Rivest & C. Stein, Capítol 22. */ 
 

	#define white 0
	#define gray 1
	#define black 2
		
	int source, ind_bf, nod, in, check_bf, counter, size_index;
	int i, j, k, s, ii, jj;
	int *colors = NULL, *queue = NULL, *current_network = NULL;
	
    int *num_veins_out = NULL, *num_veins_filtrada = NULL;
	int **veins_out = NULL, **veins_filtrada = NULL;

	colors = malloc(NN * sizeof((*colors)));
	queue = malloc(NN * sizeof((*queue)));
	current_network = malloc(NN * sizeof((*current_network)));
	
	num_veins_out = malloc(NN * sizeof((*num_veins_out)));
	num_veins_filtrada = malloc(NN * sizeof((*num_veins_filtrada)));


	
	veins_out = malloc((NN) * sizeof(*veins));
	veins_filtrada = malloc((NN) * sizeof(*veins));
	for(i=0;i<NN;i++){
		veins_out[i] = malloc((num_veins[i]+1)*sizeof(int));	//Dimensio maxima es num_veins[i], en general sera menys
		veins_filtrada[i] = malloc((num_veins[i]+1)*sizeof(int));	//Dimensio maxima es num_veins[i], en general sera menys
	}


	for(i=0;i<NN;i++){
		num_veins_filtrada[i] = 0;
		for(j=0;j<num_veins[i];j++){
			if(veins[i][j] != -1){
				veins_filtrada[i][num_veins_filtrada[i]] = veins[i][j];
				num_veins_filtrada[i]++; 
			}
		}
	}

	Check_Reprocity(veins_filtrada, num_veins_filtrada);

	for(i=0;i<NN;i++){
		if(i<NN){
			colors[i] = white;							// inicialització colors
		}
		sizes[i] = 0;
		bc[i] = 0;
		current_network[i] = 0;
		queue[i] = -1;
	}

	counter = 0;
	size_index = 0;
	ind_bf = 0;
	while(counter<NN){
		for(s=0;s<NN;s++){
			if(colors[s]==white){
				source=s;
				s=NN+1;
			}
		}
		for(i=0;i<ind_bf;i++){ 
			queue[i] = -1;								// inicialització de la queue
		}
		
//		printf("source %d\n",source);
		colors[source] = black;
		sizes[size_index]++;
		current_network[source]=1;
		counter++;

		for(i=0;i<num_veins_filtrada[source];i++){
			queue[i] = veins_filtrada[source][i];
			colors[veins_filtrada[source][i]] = gray;
			sizes[size_index]++;
			counter++;
//			printf("afegit %d\n",veins_filtrada[source][i]);
		}
			
		ind_bf = num_veins_filtrada[source];

		for(i=0;i<NN;i++){	 								//vaig mirant els elements de dins la cua
			if(queue[i] != -1){ 							//que realment sigui un element que he afegit després de la inicialització
				nod = queue[i];
	//			printf("node select %d\n", nod);
				for(j=0;j<num_veins_filtrada[nod];j++){ 				//miro els veïns d'aquest element de la cua
	//				printf("vei, color: %d %d\n", veins_filtrada[nod][j], colors[veins_filtrada[nod][j]]);
					if(colors[veins_filtrada[nod][j]]==white){				//comprovació que el node que trio no estigui dins la cua ja
															// si no és dins la cua, llavors...
						queue[ind_bf] = veins_filtrada[nod][j]; 		// ... l'entro a la cua i ...
						colors[veins_filtrada[nod][j]] = gray;		// ... li canvio el color
						ind_bf++;							// modifico l'índex d'on anirà el següent element de queue[]
						sizes[size_index]++;
						counter++;
	//					printf("afegit %d\n",veins_filtrada[nod][j]);
					}
				}
				colors[nod] = black;						// Quan ja he afegit el veïns de nod, el torno black
				current_network[nod]=1;
			} else {										// Si queue[i] != -1, sent i < NN-1, es que no tinc una global connected component
				size_index++;
				i = NN + 1;
				ReconstructionOfTheNetwork(veins_filtrada,num_veins_filtrada,
										   veins_out,num_veins_out,
										   current_network);
				//Comprova que em treu be les components...
	//			printf("%d \n\n", num_veins_out[0]);
	//			printf("%d \n\n", veins_out[0][0]);
	
				Check_Reprocity(veins_out, num_veins_out);
				BetweennesCentrality_InsideComponent(veins_out, num_veins_out, bc, current_network);

				for(k=0;k<NN;k++){
	//				printf("node, bc: %d %.3f\n", k, bc[k]);
					current_network[k]=0;
				}
			}
		}
	//	printf("counter %d\n", counter);
	}

	for(i=0;i<NN;i++){
		free(veins_out[i]);
		free(veins_filtrada[i]);
	}
	free(veins_out);
	free(veins_filtrada);
	
	for(i=0;i<NN;i++){
		if(sizes[i]!=0){
//			printf("component %d, size %d\n", i, sizes[i]);
		}
	}
	
		
	check_bf = 0;
	for(i=0;i<NN;i++){
		check_bf += colors[i];
	}


	free(colors);
	free(queue);
	free(current_network);
	free(num_veins_out);
	free(num_veins_filtrada);		
		
	if(check_bf == black*NN) {
//		printf("Connected Network!\n"); 
		return 0;
	} else {
		printf("Disconnected Network!\n"); 
	//	printf("\nExiting\n"); exit(EXIT_FAILURE);
		return 1;
	}
}

int BetweennesCentrality_withComponents(int **veins, int num_veins[], 
										int sizes[], double bc[]){
	
	/* Algorithm that finds the betweenness centrality in a network
	 * that can be formed of different disconnected components, of 
	 * different size (even of 1 node), 
	 * as well as the number of connected components and their
	 * size in a network, by means of a breadth-first search.
	 * Each node can be in three colors: white (non-discovered), 
	 * gray (discovered but still with white neighbors) and 
	 * black (discovered but no white neighbors). 
	 * 
	 * This function has been created to be used together with 
	 * ReconstructionOfTheNetwork() and 
	 * BetweennesCentrality_InsideComponent()
	 * 
	 * With litle work could give the network structure of the 
	 * disconnected components.
	 * 
	 /* Ref: Introduction to Algorithms, T.H. Cormen, C.E. Leiserson, R.L. Rivest & C. Stein, Capítol 22. */ 
 


	#define white 0
	#define gray 1
	#define black 2
		
	int source, ind_bf, nod, in, check_bf, counter, size_index;
	int i, j, k, s, ii, jj;
	int *colors = NULL, *queue = NULL, *current_network = NULL;
	
    int *num_veins_out = NULL;
	int **veins_out = NULL;

	colors = malloc(NN * sizeof((*colors)));
	queue = malloc(NN * sizeof((*queue)));
	current_network = malloc(NN * sizeof((*current_network)));
	
	num_veins_out = malloc(NN * sizeof((*num_veins_out)));
	
	veins_out = malloc((NN) * sizeof(*veins));
	for(i=0;i<NN;i++){
		veins_out[i] = malloc((num_veins[i]+1)*sizeof(int));	//Dimensio maxima es num_veins[i], en general sera menys
	}		
	
		
	for(i=0;i<NN;i++){
		if(i<NN){
			colors[i] = white;							// inicialització colors
		}
		sizes[i] = 0;
		bc[i] = 0;
		current_network[i] = 0;
		queue[i] = -1;
	}

	counter = 0;
	size_index = 0;
	ind_bf = 0;
	while(counter<NN){
		for(s=0;s<NN;s++){
			if(colors[s]==white){
				source=s;
				s=NN+1;
			}
		}
		for(i=0;i<ind_bf;i++){ 
			queue[i] = -1;								// inicialització de la queue
		}
		
//		printf("source %d\n",source);
		colors[source] = black;
		sizes[size_index]++;
		current_network[source]=1;
		counter++;

		for(i=0;i<num_veins[source];i++){
			queue[i] = veins[source][i];
			colors[veins[source][i]] = gray;
			sizes[size_index]++;
			counter++;
//			printf("afegit %d\n",veins[source][i]);
		}
			
		ind_bf = num_veins[source];


		for(i=0;i<NN;i++){	 								//vaig mirant els elements de dins la cua
			if(queue[i] != -1){ 							//que realment sigui un element que he afegit després de la inicialització
				nod = queue[i];
	//			printf("node select %d\n", nod);
				for(j=0;j<num_veins[nod];j++){ 				//miro els veïns d'aquest element de la cua
	//				printf("vei, color: %d %d\n", veins[nod][j], colors[veins[nod][j]]);
					if(colors[veins[nod][j]]==white){				//comprovació que el node que trio no estigui dins la cua ja
															// si no és dins la cua, llavors...
						queue[ind_bf] = veins[nod][j]; 		// ... l'entro a la cua i ...
						colors[veins[nod][j]] = gray;		// ... li canvio el color
						ind_bf++;							// modifico l'índex d'on anirà el següent element de queue[]
						sizes[size_index]++;
						counter++;
	//					printf("afegit %d\n",veins[nod][j]);
					}
				}
				colors[nod] = black;						// Quan ja he afegit el veïns de nod, el torno black
				current_network[nod]=1;
			} else {										// Si queue[i] != -1, sent i < NN-1, es que no tinc una global connected component
				size_index++;
				i = NN + 1;
				ReconstructionOfTheNetwork(veins,num_veins,
										   veins_out,num_veins_out,
										   current_network);
				//Comprova que em treu be les components...
	//			printf("%d \n\n", num_veins_out[0]);
	//			printf("%d \n\n", veins_out[0][0]);
				BetweennesCentrality_InsideComponent(veins_out, num_veins_out, bc, current_network);

				for(k=0;k<NN;k++){
	//				printf("node, bc: %d %.3f\n", k, bc[k]);
					current_network[k]=0;
				}
			}
		}
	//	printf("counter %d\n", counter);
	}

	for(i=0;i<NN;i++){
		free(veins_out[i]);
	}
	free(veins_out);
	
	for(i=0;i<NN;i++){
		if(sizes[i]!=0){
	//		printf("component %d, size %d\n", i, sizes[i]);
		}
	}


		
	check_bf = 0;
	for(i=0;i<NN;i++){
		check_bf += colors[i];
	}

	free(colors);
	free(queue);
	free(current_network);
	free(num_veins_out);		
		
	if(check_bf == black*NN) {
//		printf("Connected Network!\n"); 
		return 0;
	} else {
		printf("Disconnected Network!\n"); 
	//	printf("\nExiting\n"); exit(EXIT_FAILURE);
		return 1;
	}
}

void Inverse_Index(int **veins, int num_veins[], int ***inv_ind){
	
	/* 
	 * This function fills a matrix inv_ind that gives indexes, with the
	 * following criteria: let's have two nodes i and j, such that 
	 * veins[i][xx]=j, the goal is then find out the position index xx
	 * in an efficient way. This information is stored in inv_ind. It 
	 * works like:
	 * veins[env][xx] = rec -> inv[env][xx] = l -> veins[rec][l] = env
	 * 
	 * The same function of index_inverses, I just didnt want to delete
	 * the other since I used in the work of temporal correlations. Just
	 * in case... The only difference is that now I malloc the matrix
	 * inside the function.
	 * 
	 */ 

	int i, j, l;
	int env, rec;
	int flag;
	
	(*inv_ind) = malloc((NN) * sizeof(**inv_ind));
	for(i=0;i<NN;i++){
		(*inv_ind)[i] = malloc((num_veins[i]+1)*sizeof(*(*inv_ind)[i]));
		for(j=0;j<num_veins[i];j++){
			(*inv_ind)[i][j] = -1; 
		}
	}
	
	for( i = 0; i < NN; i++ ){
		for( j = 0; j < num_veins[i]; j++ ){
			if( (*inv_ind)[i][j] == -1 ){
				flag = 0;
				env = i;
				rec = veins[i][j];
				for( l = 0; l < num_veins[rec]; l++ ){
					if( veins[rec][l] == env ){
						(*inv_ind)[i][j] = l;
						flag = 1;
					}
				}
				if( flag != 1 ){
					printf("error when finding the inverses\n");
					printf("\nExiting\n"); exit(EXIT_FAILURE);
				}
			}
		}
	}
	
	return;
}

////////////////////////////////////////////////////////////////////////

void check_inputs(int Nl, double alpha, int num_initial_removals){

	if(Nl!=NN){
		printf("Wrong system size; exiting\n");
		exit(EXIT_FAILURE);
	}
	if((alpha<0)||(alpha>1)){
		printf("Wrong alpha; exiting\n");
		exit(EXIT_FAILURE);
	}
	if(num_initial_removals>=Nl){
		printf("Wrong number of initially removed nodes; exiting\n");
		exit(EXIT_FAILURE);
	}

	return;
}

int Check_Capacity(double bc[], double capacity[]){

	int i, flag;

	flag = 0;
	for(i=0;i<NN;i++){
		if(bc[i]>capacity[i]){
			flag = 1;
			i = NN + 1;
		}
	}
	return flag;
}

void LargestCC(int c, int sizes[], int *maxsize){
	
	int i, size_checker;
	
	if(c>NN){
		printf("More deletions than number of nodes... %d\n", c);
		printf("Exiting\n");
		exit(EXIT_FAILURE);
	}
	
	(*maxsize) = 0;
	size_checker = 0;
	for(i=0;i<NN;i++){
		size_checker+=sizes[i];
//				printf("%d %f %f\n", sizes[i], capacity[i], bc[i]);
		if(sizes[i]>(*maxsize)){
			(*maxsize) = sizes[i];
		}
	}
	if(size_checker!=NN){
		printf("\nWrong size checker %d %d; Exiting\n", NN, size_checker);
		exit(EXIT_FAILURE);
	}
	printf("max size %d, total %d\n", (*maxsize), size_checker);

	return;
}

void Removal_Node(int **veins, int num_veins[], 
				  int **inv_ind, double bc[], int deleted_node){
	/* deleted_node is removed.
	 * we put to 0 its bc and substitute it from the matrix of connections
	 * by a -1
	 */
	int i, j, rec;

	bc[deleted_node] = 0;
	for(j=0;j<num_veins[deleted_node];j++){
		rec = veins[deleted_node][j];
		if(rec!=-1){
			veins[deleted_node][j] = -1;
			veins[rec][inv_ind[deleted_node][j]] = -1;
		}
	}
	
	return;
}

void SimultaneousRemoval_RiskMap(int **veins, int num_veins[], 
						 int **inv_ind, double bc[], double capacity[],
						 int failed_nodes[], int *failed_node_index){
	/* Removes in block (simultaneously) all nodes in the overload
	 * condition
	 */
	int i, j, rec, node_to_remove;

	for(i=0;i<NN;i++){
		if(bc[i]>capacity[i]){
			node_to_remove = i;
			failed_nodes[(*failed_node_index)] = node_to_remove;
			(*failed_node_index)++;
			Removal_Node(veins, num_veins, inv_ind, bc, node_to_remove);
		}
	}

	return;
}

void RemovalAndCascade_Simultaneous_RiskMap(
					   int **veins, int *num_veins, int **inv_ind,
					   int sizes[], double bc[], double capacity[],
					   int failed_nodes[], char *fout,
					   int *initial_removals, int num_initial_removals){
	/* In the array initial_removals there are the labels of the nodes 
	 * that are initially removed
	 * This function runs the cascade until all nodes satisfy the 
	 * condition of holding a load (bc[]) below their capacity[].
	 * Time is measured in cascade rounds (variable c). At each round,
	 * all overloaded nodes are removed and load calculated afterwards 
	 * Once the cascade stops, I compute the GC and print it in the file
	 * In sizes[] I keep the sizes of the isolated components, although
	 * this is not used
	 */
	 
	int i, j, c, flag, failed_node_index;
	int maxsize, max_bc_node;
	double G;
	FILE *fp;

	c=0; //counter of removal rounds
	failed_node_index=0;

	//Remove initial nodes
	for(i=0;i<num_initial_removals;i++){
		Removal_Node(veins, num_veins, inv_ind, bc, initial_removals[i]);
		failed_nodes[failed_node_index] = initial_removals[i];
		failed_node_index++;
	}
	
	//Compute betweeness after the initial removal
	BetweennesCentrality_Motter(veins, num_veins, sizes, bc);
	c++;

	//Cascade goes on until all nodes are under their capacity
	while(Check_Capacity(bc, capacity)==1){
		SimultaneousRemoval_RiskMap(veins, num_veins, inv_ind, bc, capacity,
									failed_nodes, &failed_node_index);
		BetweennesCentrality_Motter(veins, num_veins, sizes, bc);
	
		printf("deletion %d; removed nodes %d\n", c, failed_node_index);
		c++;
	}

	//Compute the largest connected component
	LargestCC(c, sizes, &maxsize);
	G = (float)(maxsize)/(float)(NN);

	fp = fopen(fout,"w+");
	fprintf(fp,"%g",G);
	fclose(fp);

	return;
}

void who_fails_initially(int *initial_removals, int num_initial_removals, FILE *fp){
	
	int i, j;
	
    if(!fp){
        perror("fopen()");
        exit(EXIT_FAILURE);
    }
    
    i=0;
    fscanf(fp, "%*[^\n]");  // Read and discard a line
    while(i<num_initial_removals && fscanf(fp, "%d", &initial_removals[i++]) == 1){
		if(initial_removals[i-1]>=NN || initial_removals[i-1]<0){
			printf("\nIntially removed nodes must be between 0 and NN\n");
			exit(EXIT_FAILURE);
		}
        //~ printf("%d\n", initial_removals[i-1]);
    }
	
	printf("\nI'm removing %d nodes: ", num_initial_removals);
	for(j=0;j<num_initial_removals;j++){
		printf("%d ", initial_removals[j]);
	} printf("\n");
	
	return;
}
















