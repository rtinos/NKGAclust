/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * NK Hybrid Genetic Algorithm for Clustering
 * Copyright (C) 2018  Renato Tinos <rtinos@ffclrp.usp.br>
 * See Tinós, Zhao, Chicano & Whitley (2018). "NK Hybrid Genetic Algorithm for Clustering". IEEE Transactions on Evolutionary Computation.
 * 
 * NKGAclust is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * NKGAclust is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "defs.h"
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "nk.h"				// Nk model class

/******************************************************************************\
*				  	Print data							 			 *
\******************************************************************************/
void print_data(population *pop ){
	
	cout <<"Generation:"<<gen<<endl;
	cout <<"Best individual:"<< pop->best_individual << endl;
	cout <<"Fitness of the best individual:"<< pop->min_fitness << endl;
	cout <<"Mean fitness: "<< pop->mean_fitness << endl;
	
/*	for (int i=0;i<pop->pop_size ;i++) {	
		cout <<"("<< pop->ind[i].fitness<<") " ;
		for (int gene=0;gene<lcrom ;gene++) 
			cout << pop->ind[i].chromosome[gene]<<" ";
		cout << endl;
	}*/
}


/*********************************************************************************\
* Mapping the clusters of the red solution to the clusters of the blue solution	  *
\*********************************************************************************/
 void mapSolutions( int *solution_blue, int *solution_red ,  int *solution_red_new){
 	int i, j, k, Nc_blue=0, Nc_red=0, Nc_min, Nc_max, **M_clusters, *lin_aux, *col_aux, lin, col;
	int max_aux, *missing_red, *map_blue_red;
 	
 	// obs: 0 is used for indicating noise; the first label is 1
 	
 	// Number of labels in each solution
 	for (i=0;i<lcrom;i++){
 		if (Nc_blue<solution_blue[i]) 		
 			Nc_blue=solution_blue[i];	 	// remember that label 0 is not used
 	 	if (Nc_red<solution_red[i]) 		
 			Nc_red=solution_red[i];			// remember that label 0 is not used	
 	}
 	if (Nc_red<Nc_blue){
 		Nc_min=Nc_red;
 		Nc_max=Nc_blue;
 	}
 	else { 	
 		Nc_min=Nc_blue;	
 		Nc_max=Nc_red;
 	}
 	
	 // Memmory allocation
 	M_clusters=aloc_matrixi(Nc_blue+1,Nc_red+1);
 	map_blue_red=aloc_vectori(Nc_red+1);
 	lin_aux=aloc_vectori(Nc_blue+1);
 	col_aux=aloc_vectori(Nc_red+1);
 	missing_red=aloc_vectori(Nc_max+1);
 	
 	// Finding the labels with more hits
 	for (i=1;i<=Nc_blue;i++)
 		for (j=1;j<=Nc_red;j++) 		
 			M_clusters[i][j]=0;
 	for (i=0;i<lcrom;i++){
 		if ( solution_blue[i]>0 && solution_red[i]>0 )
 			M_clusters[ solution_blue[i] ][ solution_red[i] ] = M_clusters[ solution_blue[i] ][ solution_red[i] ] + 1;
 	}
	for (i=1;i<=Nc_blue;i++)	
		lin_aux[i]=1;
	for (i=1;i<=Nc_red;i++){
		col_aux[i]=1;
		map_blue_red[i]=0;
	}	
	 	
	// defining the mapping
	for (k=1;k<=Nc_min;k++){
		max_aux=-1;
		for (i=1;i<=Nc_blue;i++){
			if (lin_aux[i]==1){
				for (j=1;j<=Nc_red;j++){
					if (col_aux[j]==1){	
						if (M_clusters[i][j]>max_aux){
							max_aux=M_clusters[i][j];
							lin=i;
							col=j;
						}
					}
				}		
			}
		}
	
		map_blue_red[col]=lin;
		lin_aux[lin]=0;
		col_aux[col]=0;
	
	}	
		
 	// Fixing when the number of labels are different
 	for (i=1;i<=Nc_red;i++)
 		missing_red[i]=1;
 	for (i=1;i<=Nc_red;i++)
 		if (map_blue_red[i] > 0)
 			missing_red[ map_blue_red[i] ]=0;
		 
	i=1;	
	for (k=1;k<=Nc_red;k++){	
		if (map_blue_red[k]==0){	
			while (missing_red[i]==0)
				i++;			 
			map_blue_red[k]=i;
 			missing_red[i]=0; 			
 		}
 	}
 		
 	// Mapping to solution_red_new
 	for (i=0;i<lcrom;i++){
 		if (solution_red[i]==0)	
 			solution_red_new[i]=solution_red[i]; // noise
 		else
 			solution_red_new[i]=map_blue_red[ solution_red[i] ]; 		
 	}

 	
 	delete [] map_blue_red;

 	delete [] lin_aux;
 	delete [] col_aux;
 	delete [] missing_red;
 	
 	desaloc_matrixi(M_clusters, Nc_blue+1);
 	
 }


/******************************************************************************\
*								 Mutation	 						 			*
\******************************************************************************/
double mutation(nk *NK, allele *parent, allele *offspring, double fit_parent){
	double r, fit_offspring;
	
	r=random_dou ();

	if (r<0.6){		
		// reclassify
		fit_offspring=NK->mutationReclassify(parent, offspring, fit_parent);
	}
	else if (r<0.8){
		// merge clusters
		NK->mutationMerge(parent,offspring);
		fit_offspring=NK->comp_fitness(offspring);		
	}
	else{
		// split clusters
		NK->mutationSplit(parent,offspring);
		fit_offspring=NK->comp_fitness(offspring);	
	}
	
	/*double aux=NK->comp_fitness( offspring );
	if ( (aux -fit_offspring)>0.0001 || (fit_offspring -aux)>0.0001 ){
			cout<<"PROBLEM: FITNESS INCORRECT IN MUTATION!"<<endl;	
			cout<<"Fitness: "<<fit_offspring<<", "<<aux<<endl;	
			exit(1);
	}*/
	
	return fit_offspring;
}


/******************************************************************************\
*								 Crossover	 												*
\******************************************************************************/
double crossover( nk *NK, allele *parent1, allele *parent2, allele *offspring )
{
	double fit_offspring;
	int *parent2_mapped;
	
	parent2_mapped=aloc_vectori(lcrom);

	mapSolutions(parent1, parent2, parent2_mapped);			// Mapping the clusters of the red solution to the clusters of the blue solution 
	
	// Partition Crossover
	fit_offspring=NK->px(parent1,parent2_mapped,offspring);		
	NK->fixLabels(offspring);			// fix the labels
	//cout<<"Number of true partitions: "<< NK->n_true_part<<endl;	

	delete [] parent2_mapped;

	return fit_offspring;
	
}


/******************************************************************************\
*		Generation: crossing an individual with all other individuals   		 *
\******************************************************************************/
void generation( nk *NK, int n_run ){
	int gene, j=0 , parent1, parent2;
		
	do {	
		//cout<<"...Initial ind="<<j<<endl;
		// Reproduction
		if ( random_dou () > p_cross ){	
		
			// Selection of one parent	
			parent1=selection( &popold  );
			// cout<<"...mut par="<<popold.ind[parent1].fitness<<endl;	
			// Mutation
			popnew.ind[j].fitness=mutation(NK, popold.ind[parent1].chromosome, popnew.ind[j].chromosome, popold.ind[parent1].fitness);		
			// cout<<"...mut off="<<popnew.ind[j].fitness<<endl;	
		}
		else {
			
			// Selection of two parents	
			parent1=selection( &popold  );
			parent2=selection( &popold  );

			// Crossover
			popnew.ind[j].fitness=crossover( NK, popold.ind[parent1].chromosome , popold.ind[parent2].chromosome,  popnew.ind[j].chromosome );	
            //if (gen<=leng_local){	
    			if ( popnew.ind[j].fitness>popold.ind[parent1].fitness && popnew.ind[j].fitness>popold.ind[parent2].fitness )
    				n_cross_better[n_run]=n_cross_better[n_run]+1;
    			if ( popnew.ind[j].fitness>popold.ind[popold.best_individual].fitness  )
    				n_cross_best[n_run]=n_cross_best[n_run]+1;
            //}
           // cout<<"...cr off="<<popnew.ind[j].fitness<<endl;	
		}
	
		 j = j + 1;	

	} while ( j < popnew.pop_size-1);


	// Elitism
	for (gene=0;gene<lcrom;gene++)
		popnew.ind[popsize-1].chromosome[gene]=popold.ind[popold.best_individual].chromosome[gene];
	popnew.ind[popsize-1].fitness=popold.min_fitness;


}


/******************************************************************************\
*				  	Copy Population							 			 *
\******************************************************************************/
void copy_pop( void ){
	int i, gene;
	
	
	for (i=0;i<popsize ;i++) {	
		popold.ind[i].fitness=popnew.ind[i].fitness;
		for (gene=0;gene<lcrom ;gene++) 
			popold.ind[i].chromosome[gene]=popnew.ind[i].chromosome[gene];
	}
	
}

/******************************************************************************\
*				  	Generate Immigrants					 			 *
\******************************************************************************/
double immigrant(nk *NK, allele *x){
	int gene, k;
	double fit_x;
	
	// generating a random individual
	k=random_int(2, int (sqrt(lcrom)));	// number of clusters is randomly selected between 2 and sqrt(lcrom)
	for (gene=0;gene<lcrom;gene++)
		x[gene]=random_int(1, k);
	fit_x = NK->comp_fitness(x);
	//cout<<"...Immigrant: fitness before ls="<<fit_x;
	// Applying Local Search
	fit_x = NK->LsFi(x,fit_x); 
	//cout<<", fitness after ls ="<<fit_x<<endl;	
		
	return (fit_x);
}


/******************************************************************************\
*				  	Population initialization					 			 *
\******************************************************************************/
void initialization(nk *NK, char *prob_name){

	int num_ind, n_pop;
	
	// Loading individuals optimized by different algorithms to the old population	
	if (load_initsol==1)
		n_pop=load_popold(prob_name);	
	else
		n_pop=0;		
	for (num_ind=0;num_ind<n_pop;num_ind++){	
		popold.ind[num_ind].fitness=NK->comp_fitness( popold.ind[num_ind].chromosome );
		//cout<<"...Initial ind="<<num_ind<<" fitness ="<<popold.ind[num_ind].fitness<<endl;
	}
	
	// Completing with random solutions optimized by local search
	while (num_ind<popold.pop_size){	
		//cout<<"...Initial ind="<<num_ind;
		popold.ind[num_ind].fitness = immigrant(NK, popold.ind[num_ind].chromosome); 
		num_ind++;
	}
	
}

/******************************************************************************\
*				  	Run of the GA 			 *
\******************************************************************************/
void ga(nk *NK, char *prob_name, int n_run){

	int num_ind, gen_init=0, gene;
	clock_t time_start;
	
	// Statistics: crossover
	n_cross_better[n_run]=0;
	n_cross_best[n_run]=0;
	
	// Initialization
	time_start=clock();	
	gen = 0;
	initialization(NK, prob_name);				// initializing the population	
	statistics( &popold);
	//print_data(&popold);
	
	// Generations
	do {
		gen = gen + 1; 				// generation index
		// Applying local search
		if ((gen-gen_init)>leng_local){			
			/*	cout<<"gen: "<<gen<<", initial time: "<<double( clock() - time_start ) / (double)CLOCKS_PER_SEC<<endl;
				cout<<"Applying Immigrants and Local Search...."<<endl;*/
				gen_init=gen;
				for (gene=0;gene<lcrom;gene++)
					popold.ind[0].chromosome[gene]=popold.ind[popold.best_individual].chromosome[gene];
				popold.ind[0].fitness=popold.ind[popold.best_individual].fitness;
				popold.best_individual=0; 
				// Applying local search	
				for (num_ind=0;num_ind<imig_ls_ratio*popold.pop_size;num_ind++){				
					popold.ind[num_ind].fitness = NK->LsFi(popold.ind[num_ind].chromosome,popold.ind[num_ind].fitness); 					
				}
				// Generating random solutions	optimized by local search
				while (num_ind<popold.pop_size){
					popold.ind[num_ind].fitness = immigrant(NK, popold.ind[num_ind].chromosome); 
					num_ind++;
				}												
				statistics( &popold );															
		}
	
		generation(NK, n_run);	
	
		copy_pop();		// popold=popnew	
		statistics( &popold );		
		//print_data(&popold);

		
	//} while ( gen < max_gen );																	// IF STOPPING CRITERIUM IS THE MAXIMUM NUMBER OF GENERATIONS
	} while (double( clock() - time_start ) / (double)CLOCKS_PER_SEC < ((double) (lcrom)/2.0)); 	//IF STOPPING CRITERIUM IS THE MAXIMUM TIME FOR NOT LARGE DATASETS (<=1000)
	//} while (double( clock() - time_start ) / (double)CLOCKS_PER_SEC < ((double) (4.0*lcrom))); 	//IF STOPPING CRITERIUM IS THE MAXIMUM TIME FOR LARGE DATASETS (>1000)



	cout<<"gen: "<<gen<<endl;
	// Data to be saved
	time_run[n_run] = double( clock() - time_start ) / (double)CLOCKS_PER_SEC;
	for (gene=0;gene<lcrom;gene++)
		popwrite.ind[n_run].chromosome[gene]=popold.ind[ popold.best_individual ].chromosome[gene];
	popwrite.ind[n_run].fitness=popold.ind[ popold.best_individual ].fitness;
	
}


/******************************************************************************\
*				  	Main													 *
\******************************************************************************/
int main(int argc , char *argv[])
{
	int n_run, K_input;
	char *prob_name;									// name of the input file 


	// Arguments
	if( argc < 4) {
		cout<<"Insufficient number of arguments!"<<endl;
		cout<<"Call: nkpx_clustering <problem name - without extension>  <K> <load initial solutions from file>"<<endl;
		exit(1);
	}
	else{
		prob_name=argv[1];
		K_input=atoi(argv[2]);			
		load_initsol=atoi(argv[3]);
		if ( K_input<0 || (load_initsol!=0 && load_initsol!=1) ) {
			cout<<"Incorrect arguments!"<<endl;
			cout<<"Call: nkpx_clustering <problem name - without extension> <K> (K>0)  <load initial solutions from file> (0:no; 1: yes)"<<endl;
			exit(1);
		}
	}			
		
	// load dataset
	read_problem(prob_name);	
	
	// Memory Allocation
	time_run=aloc_vectord(n_runs_max);
	n_cross_better=aloc_vectori(n_runs_max);
	n_cross_best=aloc_vectori(n_runs_max);
	// populations
	aloc_pop(&popold, popsize);
	aloc_pop(&popnew, popsize);
	aloc_pop(&popwrite, n_runs_max);	

	// Initializing NK model
	nk *NK = new nk(lcrom, K_input, n_attrib, X_dataset);			// from class nk (nk.h)
	//NK->print ();													// print graph for the NK Landscape
	//NK->save (prob_name);											// save interaction graph for the NK Landscape
	
	cout << "\n ***** Genetic Algorithm ****" << endl;
	cout << "Problem: "<< prob_name <<", K="<<K_input<<", load pop="<<load_initsol<< endl;
	for (n_run=0;n_run<n_runs_max;n_run++) {	
		srand(n_run+1);	// random seed   		
		//cout <<"Run:"<< n_run <<"\n"<< endl;
		ga(NK, prob_name, n_run);											// run of the ga					
	}
		
	file_output(prob_name);													// save data

	
	// Memory Desallocation 
	delete NK;
	desaloc_pop(&popold);
	desaloc_pop(&popnew);
	desaloc_pop(&popwrite);
	desaloc_matrixd (X_dataset,lcrom);
	delete [] d_dataset;
	delete [] n_cross_best;	
	delete [] n_cross_better;
	delete [] time_run;	
	
	return 0;
}
