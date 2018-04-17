/******************************************************************************\
*				  				 Files Manipulation							 *
\******************************************************************************/

#include "defs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string>
#include<cstring>
#include<fstream>

#define CHAR_LEN 1000


/***********************************************************************\
* Read the problem instance (adapted from DTSP generator - Yang et al.) *
\***********************************************************************/
void read_problem(char *prob_name){
	int i, j;
	char line[CHAR_LEN], * keywords,Delimiters[] = " :=\n\t\r\f\v";
	char name[CHAR_LEN];

	//sprintf(name,"datasets/%s.dat", prob_name);
	sprintf(name,"%s.dat", prob_name);
	
	ifstream fin(name);
	if(fin==NULL){
			cout<<"file name error"<<endl;
			exit(0);
	}
	
	while((fin.getline(line, CHAR_LEN-1))){
			if(!(keywords = strtok(line, Delimiters)))
	  			continue;
			if(!strcmp(keywords, "N_ATTRIBUTES")){			
	  			if(!sscanf(strtok(NULL, Delimiters), "%d", &n_attrib)){
					cout<<"N_ATTRIBUTES error"<<endl;
					exit(0);
	  			}
			}
			if(!strcmp(keywords, "N_OBJECTS")){			
	  			if(!sscanf(strtok(NULL, Delimiters), "%d", &lcrom)){
					cout<<"N_OBJECTS"<<endl;
					exit(0);
	  			}
				X_dataset=aloc_matrixd (lcrom,n_attrib);
				d_dataset=aloc_vectori (lcrom);		
			}
			else if(!strcmp(keywords, "DATASET")){
	  			if(lcrom>0){
	  				for(i=0; i<lcrom; i++){	  				
						for(j=0; j<n_attrib; j++)						
							fin>>X_dataset[i][j];						
						fin>>d_dataset[i];
					}
	    		}
			}
	}
	fin.close();
	
	  
	// Test
/*	cout<<"**** Dataset, File: "<<name<<endl;
	cout<<"lcrom= "<<lcrom<<endl;
	cout<<"n_attrib= "<<n_attrib<<endl;
	cout<<"Objects:"<<endl;
	for(i=0; i<lcrom; i++){
		for(j=0; j<n_attrib; j++)
			cout<<X_dataset[i][j]<<" ";
		cout<<endl;
	}
	cout<<"Desired outputs:"<<endl;
	for(i=0; i<lcrom; i++)
		cout<<d_dataset[i]<<" ";
	cout<<endl;
	cout<<endl;*/
	
}



/****************************************************************\
* Read the initial population								     *
\****************************************************************/
int load_popold(char *prob_name){
	int i, j, k=0, n_sol, n_sol_max;
	int n_clust_alg=3, clust_alg;				// Cluster algorithm: 1 - KMeans, 2 - DBSCAN, 3 - Cluster DP			
	char line[CHAR_LEN], * keywords,Delimiters[] = " :=\n\t\r\f\v";
	char name[CHAR_LEN];

	for (clust_alg=0;clust_alg<n_clust_alg;clust_alg++){		
		n_sol=popsize/n_clust_alg;
		//sprintf(name,"datasets/%s_%d.sol", prob_name, clust_alg+1);	
		sprintf(name,"%s_%d.sol", prob_name, clust_alg+1);			
		ifstream fin(name);	
		if(fin==NULL){
				cout<<"file name (solutions) error"<<endl;
				exit(0);
		}
		
		while(fin.getline(line, CHAR_LEN-1)){
				if(!(keywords = strtok(line, Delimiters)))
		  			continue;
		  		if(!strcmp(keywords, "N_SOL")){
			  		if(!sscanf(strtok(NULL, Delimiters), "%d", &n_sol_max)){
						cout<<"N_SOL error"<<endl;
						exit(0);
		  			}		  		
					if (n_sol>n_sol_max)		  			
		  				n_sol=n_sol_max;
					if (popsize<k+n_sol)
						n_sol=popsize-k;							    		
				}
				if(!strcmp(keywords, "SOLUTIONS")){
		  			for(i=0; i<n_sol; i++){		  			
		  				for(j=0; j<lcrom; j++){		  			  								
							fin>>popold.ind[k].chromosome[j];					
							if (popold.ind[k].chromosome[j]==-1)
								popold.ind[k].chromosome[j]=0; 	
						}
						k++;
					}
				}
				
		}
		fin.close();
	}

	// Test
/*	cout<<"**** Inital Population, File: "<<name<<endl;
	cout<<"popsize= "<<k<<endl;
		cout<<"Solutions:"<<endl;	
	for(i=0; i<k; i++){	
		for(j=0; j<lcrom; j++)	  
			cout<<popold.ind[i].chromosome[j]<<" ";
		cout<<endl;
	}
	cout<<endl;*/
	
	return k;	  
}



/******************************************************************************\
* 										Save data : end of the simulation						 *
\******************************************************************************/
void file_output(char *prob_name)
{
	int i, j;
	FILE *Sol_file;
	FILE *Bestfit_file;
	//FILE *Time_file;
	//FILE *Imp_file;
	//FILE *Impb_file;
	char *name_p;
	char name[CHAR_LEN];

   name_p = name;

	// Best solution
	//sprintf(name,"data/sol_%s.dat",prob_name);
	sprintf(name,"sol_%s.dat",prob_name);
	if ((Sol_file = fopen(name_p,"w"))==NULL) {
		puts("The file sol to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<popwrite.pop_size;i++){	
		for (j=0;j<lcrom;j++) 
			fprintf(Sol_file,"%d ",popwrite.ind[i].chromosome[j]);
		fprintf(Sol_file,"\n");
	}
	fclose(Sol_file);

   // Best fitness 
	//sprintf(name,"data/bfi_%s.dat",prob_name);
	sprintf(name,"bfi_%s.dat",prob_name);
	if ((Bestfit_file = fopen(name_p,"w"))==NULL) {
		puts("The file bfi to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<popwrite.pop_size;i++) 
		fprintf(Bestfit_file,"%1.5f ",popwrite.ind[i].fitness);
	fclose(Bestfit_file);


  	/*
	// Time for each run
	sprintf(name,"data/time_%s.dat",prob_name);
	if ((Time_file = fopen(name_p,"w"))==NULL) {
		puts("The file time to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<n_runs_max;i++) 
		fprintf(Time_file,"%1.3f ",time_run[i]);
	fclose(Time_file);

	// Number of improvements caused by crossover 
	sprintf(name,"data/imp_%s.dat",prob_name);
	if ((Imp_file = fopen(name_p,"w"))==NULL) {
		puts("The file imp to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<n_runs_max;i++) 
		fprintf(Imp_file,"%d ",n_cross_better[i]);
	fclose(Imp_file);

	// Number of improvements of the best fitness caused by crossover 
	sprintf(name,"data/imb_%s.dat",prob_name);
	if ((Impb_file = fopen(name_p,"w"))==NULL) {
		puts("The file impb to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<n_runs_max;i++) 
		fprintf(Impb_file,"%d ",n_cross_best[i]);
	fclose(Impb_file);
	*/

}



