/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * Class: NK clustering (used in NK Hybrid Genetic Algorithm for Clustering)
 * Copyright (C) 2018  Renato Tinos <rtinos@ffclrp.usp.br>
 * See Tinós, Zhao, Chicano & Whitley (2018). "NK Hybrid Genetic Algorithm for Clustering". IEEE Transactions on Evolutionary Computation.
 * 
 * nk.h is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * nk.h is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.*/

/******************************************************************************\
*								 Class: NK clustering						 *
\******************************************************************************/
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include "Grafo.h"							// Graph class
#include "BuscaEmProfundidade.h"			// Deep First Search  class		(it is used inn compGraph) 	

using cap7_listaadj_autoreferencia::Grafo; 
#define CHAR_LEN 1000

class nk {
private:	
	Grafo *G_ep;											// Epistasis graph 
	Grafo *G_ep_transp;										// Transpose of the Epistasis graph 
	double **D;												// Distance matrix (distance between objects of the dataset) 
	int N;													// number of objects
	int K;													// epistasis degreee
	int n_at;												// dimension: number of atributtes
	double *rho;											// density vector 	
	double dt0;												// distance parameter 0: used to compute NK fitness function	
	double dt1;												// distance parameter 1: used to compute NK fitness function
	double dt2;												// distance parameter 2: used to compute NK fitness function
	double d_rho;											// density parameter: used to compute NK fitness function	
	typedef struct{
		int size;	   										// number of nodes in the candidate partition
		int test;   										// 1 if it is a true recombining partition and 0 otherwise
		double cost_blue;									// cost of the candidate partition (if it is a true recombining partition) - blue solution
		double cost_red;									// cost of the candidate partition (if it is a true recombining partition) - red solution
		int select;											// 0 if selected from blue solutionn and 1 from red solution
		int label;											// label for the cluster (used in fixLabels)
	} candidate;											// candidate partitions: used in the partiton crossover
	double comp_fi(int *x, int i);							// Compute fitness subfunction fi	
	void threeSigmaRule(double **X_dataset);				// 3-Sigma Rule: used to comptute dt1 and dt2	
	double dfElement (int *x, int i, int xi_old);			// Difference in the fitness caused by changing element i	
public:
	int n_true_part;															// number of true candidate partitions (for the last application of crossover)
	nk(int N_input, int K_input, int n_attrib_input, double **X_dataset);
    ~nk(void);
	double comp_fitness (int *x);												// Compute NK custering internal evaluation criterion (fitness)
	double px(int *solution_blue, int *solution_red, int *offspring );			// Partition crossover for clustering
	void fixLabels(int *offspring);												// Fix labels (create new label for isolated clusters with same label)
	double LsFi (int *x, double f);												// Local search
	double mutationReclassify(int *parent, int *offspring, double fit_parent );	// Mutation 1: Reclassify Points	
	void mutationMerge(int *parent, int *offspring );							// Mutation 2: Merge Clusters
	void mutationSplit(int *parent, int *offspring );							// Mutation 3: Split Clusters
	void print(void);															// Print NK information	
	void save(char *prob_name);													// Save  NK information (transpose of the interaction graph)
};


/******************************************************************************\
*								Constructor									   *
\******************************************************************************/
nk::nk(int N_input, int K_input, int n_attrib_input, double **X_dataset){
	int i, j, k, *deg_in, ind_min, percent, position, **M_ep;
	double max_D1=0.0, **D_aux, *vD, delta, dc, rho_max, m_rho, s_rho;										
										
	N = N_input;
	K = K_input;
	n_at = n_attrib_input;
		
	// Graphs (allocation)
	G_ep = new Grafo (N); 					
	G_ep_transp = new Grafo (N); 
	
	//	Memmory allocation (obs.: some matrices and vector are allocated in other places)
	M_ep=aloc_matrixi (N,N);	// Epistasis matrix (obs.: we keep both representations - graph and matrix) 
	D=aloc_matrixd (N,N);	
	deg_in=aloc_vectori(N);
	vD=aloc_vectord( ( (N*(N-1))/2 ) );	// vector with all distances
	rho=aloc_vectord(N);

	// Distance matrix
	k=0;
	for (i=0;i<N;i++){
		D[i][i]=0.0;
		M_ep[i][i]=0;
		for (j=i+1;j<N;j++){		
			D[i][j]=sqrEucDist(X_dataset,i,j,n_at);					//squared Euclidean distance
			D[j][i]=D[i][j];
			vD[k]=D[i][j];
			k++;
			M_ep[i][j]=0;
			M_ep[j][i]=0;
			if (D[i][j]>max_D1)
				max_D1=D[i][j];									// maximum element of the matrix
		}
	}
	max_D1=max_D1+1.0;									// maximum element of the matrix + 1

	// Computing the density rho
	// Defining the cutoff distance dc (see ref. Rodriguez & Laio (2014))
	percent=2;
	position=(( (N*(N-1))/2 )*percent)/100;
	quicksort(vD, ( (N*(N-1))/2 ));
	dc=vD[position];		// cutoff distance: used to compute the density
	//cout<<"dc="<<dc<<endl;
	// Defining the density rho of each point( see ref. Rodriguez & Laio (2014))
	for (i=0;i<N;i++)
		rho[i]=0.0;	
	// Gaussian Kernel	
	for (i=0;i<N-1;i++){	
		for (j=i+1;j<N;j++){
			rho[i]=rho[i]+exp(-(D[i][j]/dc)*(D[i][j]/dc));
     		rho[j]=rho[j]+exp(-(D[i][j]/dc)*(D[i][j]/dc));
		}
	}
	// Scaling rho (between 0 and 1) and computing mean and std of rho
	rho_max=rho[0];				// maximum value of rho
	for (i=1;i<N;i++)
		if (rho[i]>rho_max)
			rho_max=rho[i];
	rho[0]=rho[0]/rho_max;
	m_rho=rho[0]; 			// mean of rho
	for (i=1;i<N;i++){
		rho[i]=rho[i]/rho_max;	
		m_rho=m_rho+rho[i];
	}
	m_rho=m_rho/N;
	s_rho=0.0;					// std of rho
	for (i=0;i<N;i++)
		s_rho=s_rho+(rho[i]-m_rho)*(rho[i]-m_rho);
	s_rho=sqrt(s_rho/N);
	d_rho=m_rho-s_rho; 			//68–95–99.7 rule or three-sigma rule of thumb:
	if (d_rho<=0.0)
		d_rho=m_rho/2.0;
	
	// Defining the Epistasis Graph	
	// Step 1: Connecting each object to the closest object with higher density
	for (i=0;i<N;i++){
		deg_in[i]=0;
		delta=max_D1;
		for (j=0;j<N;j++){
			if (rho[i]<rho[j]){			
				if (D[i][j]<delta){				
					delta=D[i][j];
					ind_min=j;
				}
			}		
		}
		if (delta==max_D1){	
			// Object with highest density (
			if (i>0)
				ind_min=0;	
		    else
		    	ind_min=1;
			for (j=0;j<N;j++){			
				//	if (D[i][ind_min]<D[i][j]) // definition from Rodriguez & Laio (2014): we use another definition here
				if (j!=i && D[i][ind_min]>D[i][j])
					ind_min=j;
			}
		}
		G_ep->insereAresta (ind_min, i, 1);				// insert edge: ind_min influences fi
		G_ep_transp->insereAresta (i, ind_min, 1);		// insert edge: list i: all nodes that influence fi
		M_ep[ind_min][i]=1;		
		deg_in[i]=deg_in[i]+1;		
	}	
		
	// Step 2: Adding edges to each vertex (according to the distance) until the indegree is equal to K
    delete [] vD;	
	D_aux=aloc_matrixd (N,N);
	for (i=0;i<N;i++){
		for (j=0;j<N;j++){		
			D_aux[i][j]=D[i][j];
			if (i==j)
				D_aux[i][j]=D_aux[i][j]+max_D1;
		}
	}
	for (i=0;i<N;i++){		
		while ( deg_in[i]<K ){				
			ind_min=0;
			for (j=1;j<N;j++){
				if (D_aux[i][j]<D_aux[i][ind_min])
					ind_min=j;
			}
			
			D_aux[i][ind_min]=max_D1;
			if ( M_ep[ind_min][i]==0 ){
				G_ep->insereAresta (ind_min, i, 1);				// insert edge: ind_min influences fi
				G_ep_transp->insereAresta (i, ind_min, 1);		// insert edge: list i: all nodes that influence fi
				M_ep[ind_min][i]=1;
				deg_in[i]=deg_in[i]+1;
			}		
		}	
	}
	// Compute parameters dt1 and dt2
	threeSigmaRule(X_dataset);
	
	desaloc_matrixi (M_ep,N);
    desaloc_matrixd (D_aux,N);
	delete [] deg_in;
    
}


/******************************************************************************\
*								 Destructor													   *
\******************************************************************************/
nk::~nk(void){

	desaloc_matrixd (D,N);
	delete G_ep;
	delete G_ep_transp;
	delete [] rho;
               
}  


/******************************************************************************\
*				3-Sigma Rule: used to comptute dt1 and dt2					   *
\******************************************************************************/
void nk::threeSigmaRule(double **X_dataset){
	int i, j, k, flag_null;
	double *xg, m_adj, s_adj;
	//Grafo::Aresta *adj;
	
	xg=aloc_vectord(N*K);

	k=0;
	for (i=0;i<N;i++){
		flag_null=0;		
		Grafo::Aresta *adj = G_ep_transp->primeiroListaAdj(i);	  // first edge to node i
		if (adj == NULL)
			flag_null=1;
		else			
			j=adj->_v2 ();		
		while (flag_null==0) {
			xg[k]=D[i][j];
			k++;
			delete adj;
			adj = G_ep_transp->proxAdj(i);	  // next edge to node i
			if (adj == NULL)
				flag_null=1;
			else			
				j=adj->_v2 ();	
		}			
		
	}
	m_adj=0.0;
	for (i=0;i<N*K;i++)
		m_adj=m_adj+xg[i];
	m_adj=m_adj/(N*K);
	s_adj=0.0;
	for (i=0;i<N*K;i++)
		s_adj=s_adj+(xg[i]-m_adj)*(xg[i]-m_adj);
	s_adj=sqrt(s_adj/(N*K));
	
	//dt0=m_adj+s_adj; 	//rho: 68.27%, 95.45% and 99.73% of the values lie within one, two and three standard deviations of the mean, respectively
	//dt1=m_adj+2*s_adj; 
    //dt2=m_adj+3*s_adj;
   	dt0=m_adj; 	//rho: 68.27%, 95.45% and 99.73% of the values lie within one, two and three standard deviations of the mean, respectively
	dt1=m_adj+2*s_adj; 
   	dt2=m_adj+s_adj;
       		
	delete [] xg;
}


/******************************************************************************\
*								Compute contribution of subfunction fi		   *
\******************************************************************************/
double nk::comp_fi(int *x, int i){
	int j, flag_null;
	double fi=0.0, fmax, fmin, fmax0=1.0, fmin0=0.0;

	// obs.: the first label is 1  (label 0 indicates noise)

	flag_null=0;		
	Grafo::Aresta *adj = G_ep_transp->primeiroListaAdj(i);	  // first edge to node i
	if (adj == NULL)
		flag_null=1;
	else			
		j=adj->_v2 ();		
	while (flag_null==0) {
		fmax=rho[j]*fmax0;
		fmin=rho[j]*fmin0;
		if (x[i]==0){
			if (rho[i]<=d_rho && D[i][j]>dt2)
				fi=fi+fmin;	
			else
				fi=fi+fmax;					
		}
		else if (x[i]==x[j]){
			if	(D[i][j]>dt0 && D[i][j]<=dt1)	
				fi=fi+fmin+(fmax-fmin)*(D[i][j]-dt0)/(dt1-dt0);
			else if	(D[i][j]>dt1)
				fi=fi+fmax;	
			else
				fi=fi+fmin;			
		}
		else{
			if	(D[i][j]<=dt0)	
				fi=fi+fmax;
			else if	(D[i][j]>dt0 && D[i][j]<=dt1)
				fi=fi+fmax-(fmax-fmin)*(D[i][j]-dt0)/(dt1-dt0);		
			else
				fi=fi+fmin;	
		}
	
		delete adj;
		adj= G_ep_transp->proxAdj(i);	  // next edge o node i
		if (adj == NULL)
			flag_null=1;
		else			
			j=adj->_v2 ();	
	}	
	
	// Scaling	
	fi=fi/(N*K);
		
	
	return (fi);
}


/******************************************************************************\
*	Compute NK custering internal evaluation criterion (fitness)			   *
\******************************************************************************/
double nk::comp_fitness (int *x){
	int i;
	double f=0.0, fi;
	
	for (i=0;i<N;i++){	
	
		fi=comp_fi(x,i);		// Computing fi				
		f=f+fi;					// Sum of fi	
	//	cout<<x[i]<<"("<<fi<<") ";
	}
				
//	cout<<endl;
	return f;
}


/******************************************************************************\
*					Partition Crossover for clustering						   *
\******************************************************************************/
double nk::px(int *solution_blue, int *solution_red, int *offspring ){
	int i, k, flag_null, n_cand, *comp_id;	
	double fit_offspring, common_cost, f_aux;
	//Grafo::Aresta *adj;

	n_true_part=0;			// number of true candidate partitions 

	// Identifying the recombining partitions (Phase 1)
	// Step 1.1.: Creating the symmetric union graph without the edges of the common vertices  
	Grafo *G_p1 = new Grafo (N);			// from class graph (grafo.h)
	for (i = 0; i <N; i++) {		
		if ( solution_blue[i] != solution_red[i] ){
			flag_null=0;
			Grafo::Aresta *adj = G_ep->primeiroListaAdj(i);	  // first edge of the vertex	
			if (adj == NULL)
				flag_null=1;
			else 
				k=adj->_v2 ();
			while (flag_null==0) {				
				G_p1->insereAresta (i, k, 1);				// insert edge: i influences k				
				G_p1->insereAresta (k, i, 1);				// insert edge: k influences i (symmetric)
				delete adj;
				adj = G_ep->proxAdj(i);	  // next edge of the vertex	
				if (adj == NULL)
					flag_null=1;					
				else
					k=adj->_v2 ();			
			}
		}
	}
	
	
	// Step 1.2: Finding the connected components of the union graph G_p1
	comp_id=new int[N];							//candidate component for each node 	
	cap7::BuscaEmProfundidade cfc(G_p1);					// object connected componets (from BuscaEmProfundidade.h)
	cfc.compCon (comp_id);							// find the connected components for the graph
	cfc.~BuscaEmProfundidade ();	
	//cap7::Cfc cfc(G_p1);								// from class cfc (Cfc.h)
    //cfc.obterCfc (comp_id);					// find the connected components for the graph
	//cfc.~Cfc ();
	delete G_p1;
	// Step 1.3: Size of the candidates
	n_cand=0;				// number of candidates
	for (i=0;i<N;i++)
			if ( comp_id[i]>n_cand )
				n_cand=comp_id[i];	
	n_cand++;	   // remember that the first component has label 0
	candidate *candidates=new candidate[n_cand];							// candidates		
	for (i=0;i<n_cand;i++){
			candidates[i].size=0;
			candidates[i].test=0;
			candidates[i].cost_blue=0.0;
			candidates[i].cost_red=0.0;
	}
	for (i=0;i<N;i++)
			candidates[ comp_id[i] ].size=candidates[ comp_id[i] ].size + 1;	
	// Step 1.4: Testing the candidates		
	for (i=0;i<N;i++)
			if ( (candidates[ comp_id[i] ].size >1) || (solution_blue[i] != solution_red[i]) ){
				candidates[ comp_id[i] ].test=1;
			}
	for (i=0;i<n_cand;i++)
		if ( candidates[i].test == 1)
			n_true_part++;

	// Generating the best offspring (Phase 2)
	// Step 2.1: computing the cost of each partition for the blue solution		
	common_cost=0.0;
	for (i=0;i<N;i++){
		f_aux=comp_fi(solution_blue,i);		
		if (candidates[comp_id[i]].test==1)
			candidates[comp_id[i]].cost_blue = candidates[comp_id[i]].cost_blue + f_aux;
		else
			common_cost = common_cost + f_aux;
	}
	// Step 2.2: computing the cost of each partition for the red solution		
	for (i=0;i<N;i++){
		f_aux=comp_fi(solution_red,i);		
		if (candidates[comp_id[i]].test==1)
			candidates[comp_id[i]].cost_red = candidates[comp_id[i]].cost_red + f_aux;	
	}
	// Step 2.3: generating the offspring
	fit_offspring=common_cost;		// first offspring (best individual)
	for (i=0;i<n_cand;i++){
		if ( candidates[i].test == 1){
			if ( candidates[i].cost_blue < candidates[i].cost_red ){
				candidates[i].select = 0;
				fit_offspring=fit_offspring+candidates[i].cost_blue;
			}
			else {
				candidates[i].select = 1;
				fit_offspring=fit_offspring+candidates[i].cost_red;
			}
		} 	
		else
			candidates[i].select = 0;	  // does not matter the solution: the cost for the common partition is the same
	}
	for (i=0;i<N;i++){
		if (candidates[comp_id[i]].select==0){		
			offspring[i]=solution_blue[i];
		}
		else{
			offspring[i]=solution_red[i];
		}
	}
				
	delete [] candidates;	
	delete [] comp_id;

	return fit_offspring;
}


/******************************************************************************\
* Fix labels (create new label for isolated clusters with same label)		   *
\******************************************************************************/
void nk::fixLabels(int *offspring){
	int i, j, k, *comp_id, *n_labels, n_comp, label_max, flag_null, new_label;
			
	// Step 1: Creating a symmetric union graph without the edges between vertices with different labels 
	Grafo *G_u = new Grafo (N);			// from class graph (grafo.h)
	for (i = 0; i <N; i++) {				
		flag_null=0;
		Grafo::Aresta *adj = G_ep->primeiroListaAdj(i);	  	// first edge of the vertex	
		if (adj == NULL)
			flag_null=1;
		else 
			k=adj->_v2 ();
		while (flag_null==0) {		
			if (offspring[i]==offspring[k]){			
				G_u->insereAresta (i, k, 1);				// insert edge: i influences k				
				G_u->insereAresta (k, i, 1);				// insert edge: k influences i (symmetric)
			}
			delete adj;
			adj = G_ep->proxAdj(i);	  // next edge of the vertex	
			if (adj == NULL)
				flag_null=1;					
			else
				k=adj->_v2 ();						
		}		
	}
	// Step 2: Finding the connected components of the union graph G_u
	comp_id=aloc_vectori(N);			// candidate component for each node 
	cap7::BuscaEmProfundidade cfc(G_u);					// object connected componets (from BuscaEmProfundidade.h)
	cfc.compCon (comp_id);							// find the connected components for the graph
	cfc.~BuscaEmProfundidade ();	
	//cap7::Cfc cfc(G_u);										// from class cfc (Cfc.h)
    //cfc.obterCfc (comp_id);									// find the connected components for the graph
	//cfc.~Cfc ();
	// Step 3: Computing the number of times that labels appear in different clusters
	// Finding the number of clusters
	label_max=0;
	for (i=0;i<N;i++)
		if (offspring[i]>label_max)
			label_max=offspring[i];
	if (label_max<=1){	
		delete G_u;
		delete [] comp_id;
		return;			
	}
	// Defining the labels of the partitions
	n_comp=0;				// number of connected components
	for (i=0;i<N;i++)
		if ( comp_id[i]>n_comp )
			n_comp=comp_id[i];	
	n_comp++;	   											// remember that the first component has label 0
	candidate *components=new candidate[n_comp];			//  connected components		
	for (i=0;i<N;i++)		
		components[comp_id[i]].label=offspring[i];
	// Computing the number of times for each label
	n_labels=aloc_vectori(label_max+1);			// numbel of times that each label appear in isolated clusters
	for (i=1;i<=label_max;i++)
		n_labels[i]=0;
	for (i=0;i<n_comp;i++)	
		if (components[i].label>0)
			n_labels[components[i].label]=n_labels[components[i].label]+1;
	// Step 4: Defining new labels for clusters with repeated labels
	new_label=label_max+1;
	for (i=1;i<=label_max;i++){
		j=0;
		while (n_labels[i]>1){
			// Finding first occurence
			while (j<n_comp && components[j].label!=i)
				j++;
			// Finding next ocurrences
			while (j<n_comp && components[j].label!=i)
				j++;
			if (j<n_comp){			
				components[j].label=new_label;
				new_label=new_label+1;
			}
			n_labels[i]=n_labels[i]-1;				
		}
	}
	for (i=0;i<N;i++)
		offspring[i]=components[comp_id[i]].label;

	
	delete [] components;	
	delete G_u;
	delete [] comp_id;
	delete [] n_labels;
}


/******************************************************************************\
*					Difference in the fitness caused by changing element i	   *
\******************************************************************************/
double nk::dfElement (int *x, int i, int xi_old){
	int j, flag_null=0, xi_new;
	double df=0.0;
	//Grafo::Aresta *adj;
	
	xi_new=x[i];
	
	// Adding the i-th contribution of the new element i
	df=df+comp_fi(x,i);
	// Removing the i-th contribution of the old element i
	x[i]=xi_old;
	df=df-comp_fi(x,i);
	x[i]=xi_new;
	
	// Other elements influenced by bit i
	Grafo::Aresta *adj = G_ep->primeiroListaAdj(i);	  // first edge of the vertex	(component influenced by bit i)
	if (adj == NULL)
		flag_null=1;
	else 
		j=adj->_v2 ();
	while (flag_null==0) {
		// Adding the j-th contribution of the new element i 
		df=df+comp_fi(x,j);		
		// Removing the j-th contribution of the old element i
		x[i]=xi_old;
		df=df-comp_fi(x,j);
		x[i]=xi_new;

		delete adj;
		adj = G_ep->proxAdj(i);	  // next edge of the vertex	
		if (adj == NULL)
			flag_null=1;	
		else 
			j=adj->_v2 ();		
	}

	return df;
}


/******************************************************************************\
*					Local Search: First Improvement						   	   *
\******************************************************************************/
double nk::LsFi (int *x, double f){
	int i, j, k, count_it=0, count_plato=0, flag_cont=1, *v_in, *v_out, xi_old, flag_null, l, flag_pointer_adj;
	double df;
	// Grafo::Aresta *adj;
	
	// Randomly defining the order of the search
	v_in=aloc_vectori(N);
	v_out=aloc_vectori(N);
	for (i=0;i<N;i++)
		v_in[i]=i;
	rand_perm(v_in, v_out, N);

	// Loop
	while (flag_cont==1 && count_it<N/2){
		
		// finding a best (or with same fitness) neigbour of the current best solution
		df = 0.0;
		j=0;		
		while(j<N && df>=0.0){
			i=v_out[j];
			xi_old=x[i];		// value of x[i] before the changing			
		
			// Changing element x[i] (according to its neighbours)
			flag_null=0;		
			Grafo::Aresta *adj = G_ep_transp->primeiroListaAdj(i);	  // first edge to node i
			flag_pointer_adj=0;
			if (adj == NULL)
				flag_null=1;
			else			
				k=adj->_v2 ();		
			while (flag_null==0) {
				// changing x[i] to the label of x[k]
				x[i]=x[k];
				//cout<<k<<","<<x[i]<<","<<flag_pointer_adj<<endl	;		

				// Computing the difference in the fitness
				delete adj;
				df=dfElement(x,i,xi_old);							
				/*double aux=comp_fitness( x );
				for (int ii=0;ii<N;ii++)
					v_in[ii]=x[ii];
				v_in[i]=xi_old;	
				double aux2=aux-comp_fitness( v_in );
			
				if ( (df - aux2)>0.0001 || (aux2 -df)>0.0001 ){
					cout<<"PROBLEM: FITNESS INCORRECT IN LS AAXXXXX!"<<endl;	
					cout<<"Fitness: "<<aux2<<", "<<df<<endl;	
					df=dfElement2(x,i,xi_old);	
					exit(1);
				}	*/		
				if (df>=0.0){	
					if (df>0.0)							
						x[i]=xi_old;		// Reversing the changing		
					// next neighboir					
					l=0;
					while(l<=flag_pointer_adj){
						if (l==0)
							adj = G_ep_transp->primeiroListaAdj(i);	  // first edge to node i
						else{
							delete adj;						
							adj = G_ep_transp->proxAdj(i);	  // next edge o node i
						}
						l++;
					}					
					delete adj;	
					adj = G_ep_transp->proxAdj(i);	  // next edge o node i
					flag_pointer_adj++;
					if (adj == NULL)
						flag_null=1;
					else			
						k=adj->_v2 ();	
						
				}
				else
					flag_null=1;
			}	
							
			j++;
		}

		if (df<=0.0){	
				
			if (df<0.0){
				count_plato=0;
				f=f+df;		// new fitness value
			}
			else{
				// checking platos
				count_plato=count_plato+1;
				 if (count_plato>3)				 
			 		flag_cont=0;
			}		
		}	
		else{	
			flag_cont=0;			
		}
		
		count_it++;
		
	}	

	delete [] v_out;
	delete [] v_in;
	
	return f;

}


/******************************************************************************\
*					Mutation 1: Reclassify Points						   	   *
\******************************************************************************/
double nk::mutationReclassify(int *parent, int *offspring, double fit_parent){
	int i, ri, k, flag_count=0, flag_null, offspringi_old, best_mut, l, flag_pointer_adj;
	double df, df_min=0.0;
	//Grafo::Aresta *adj;
	
	for (i=0;i<N;i++)
		offspring[i]=parent[i];
	
	ri=random_int(0, N-1);			// selecting randomly one element of the vector to be mutated
	offspringi_old=offspring[ri];
	
	// Changing element offspring[i] (according to its neighbours)
	flag_null=0;		
	Grafo::Aresta *adj = G_ep_transp->primeiroListaAdj(ri);	  // first edge to node ri
	flag_pointer_adj=0;
	if (adj == NULL)
		flag_null=1;
	else			
		k=adj->_v2 ();		
	while (flag_null==0) {
		// changing offspring[i] to offspring[k]
		offspring[ri]=offspring[k];					
		// Computing the difference in the fitness
		df=dfElement(offspring,ri,offspringi_old);						
		if (flag_count==0){
			df_min=df;
			best_mut=offspring[ri];
			flag_count=1;	
		}
		else{
			if (df<df_min){
				df_min=df;
				best_mut=offspring[ri];	
			}			
		}				
		
		delete adj;				
		// next neighboir
		l=0;
		while(l<=flag_pointer_adj){
			if (l==0)
				adj = G_ep_transp->primeiroListaAdj(ri);	  // first edge to node i
			else{	
				delete adj;		
				adj = G_ep_transp->proxAdj(ri);	  // next edge o node i
			}
			l++;
		}
		delete adj;	
		adj = G_ep_transp->proxAdj(ri);	  // next edge o node i
		flag_pointer_adj++;
		if (adj == NULL)
			flag_null=1;
		else			
			k=adj->_v2 ();	
	}
	// Testing mutation to '0' (noise)		
	offspring[ri]=0;
	df=dfElement(offspring,ri,offspringi_old);						
	if (flag_count==0 || df<df_min){		
		best_mut=offspring[ri];	
		df_min=df;	
	}	
	
	// Changing offspring[ri] to the best mutation
	offspring[ri]=best_mut;	
	
	return (fit_parent+df_min);
}


/******************************************************************************\
*					Mutation 2: Merge Clusters						   		   *
\******************************************************************************/
void nk::mutationMerge(int *parent, int *offspring){
	int i, l, label_max, *prototype, pcluster1=-1, pcluster2=-1;
	double aux, r, *p_chosen;
	
	// Finding the number of clusters
	label_max=0;
	for (i=0;i<N;i++){
		offspring[i]=parent[i];
		if (parent[i]>label_max)
			label_max=parent[i];
	}
	if (label_max<=1)
		return;	

	// Finding the prototypes of each cluster (object with highest density)
	prototype=aloc_vectori(label_max+1);
	for (l=0;l<=label_max;l++)
		prototype[l]=-1;
	for (i=0;i<N;i++){
		if (parent[i]>0){
			if (prototype[parent[i]]==-1)	
				prototype[parent[i]]=i;
			else{
				if ( rho[i]>rho[prototype[parent[i]]] )
					prototype[parent[i]]=i;				
			}							
		}		
	}
	
	// Finding the two clusters with closest prototypes
	/*
	int j;
	for (l=1;l<label_max;l++){	
		for (j=l+1;j<=label_max;j++){
			if (prototype[l]>=0 && prototype[j]>=0){
				if (pcluster1<0){
					pcluster1=l;
					pcluster2=j;
				}
				else {
					if ( D[prototype[l]][prototype[j]]<D[prototype[pcluster1]][prototype[pcluster2]] ){
						pcluster1=l;
						pcluster2=j;
					}					
				}								
			}	
		}
	}*/
	// Finding two clusters: the first one is random and the second has probability of being select inversely proportional to its distance to first cluster	
	p_chosen=aloc_vectord(label_max+1);
	aux=-1.0;
	for (l=1;l<=label_max;l++){
		p_chosen[l]=0.0;
		if (prototype[l]>=0){
			r=random_dou();
			if (r>aux){
				aux=r;
				pcluster1=l;
			}
		}
	}
	aux=0.0;
	for (l=1;l<=label_max;l++){
		if ( prototype[l]>=0 && l!=pcluster1 && D[prototype[l]][prototype[pcluster1]]>0.0){
			p_chosen[l]=1.0/D[prototype[l]][prototype[pcluster1]];
			aux=aux+p_chosen[l];
		}
	}
	for (l=1;l<=label_max;l++)
		p_chosen[l]=p_chosen[l]/aux;
	r=random_dou();
	aux=0.0;
	l=1;
	while(l<=label_max && aux<=r){
		aux=aux+p_chosen[l];
		l++;
	}
	if (l<=label_max)
		pcluster2=l-1;
	else
		pcluster2=label_max;
	
	// Merging the two clusters
	if (pcluster1>0 && pcluster2>0)
		for (i=0;i<N;i++)
			if (parent[i]==pcluster2)
				offspring[i]=pcluster1;
			
	delete [] prototype;
	delete [] p_chosen;
	
}


/******************************************************************************\
*					Mutation 3: Split Clusters						   	       *
\******************************************************************************/
void nk::mutationSplit(int *parent, int *offspring){
	int i, l, label_max, *size_clust, ri, partial_sum=0, size_total=0, chosen=0;
	int new_label, prototype1, prototype2;
	double auxd, r1, r2, *p_chosen;
	
	// Finding the number of clusters
	label_max=0;
	for (i=0;i<N;i++){
		offspring[i]=parent[i];
		if (parent[i]>label_max)
			label_max=parent[i];
	}
	if (label_max<1)
		return;

	// Computing the probability of spliting each cluster (according to its size)
	size_clust=aloc_vectori(label_max+1);
	for (i=0;i<=label_max;i++)
		size_clust[i]=0;
	for (i=0;i<N;i++){	
		if (parent[i]>0){		
			size_clust[parent[i]]=size_clust[parent[i]]+1;	// size of each cluster
			size_total=size_total+1;						// number of all points minus number of points with label ´0'
		}
	}
	
	// Selecting the cluster to be splitted according to its probability
	ri=random_int(1, size_total);
	l=1;
	while (l<=label_max && chosen==0){	
		partial_sum=partial_sum+size_clust[l];
		if (ri<=partial_sum)
			chosen=l;
		l++;		
	}

	if (size_clust[chosen]>1){	
		// Finding two prototypes in the chosen cluster according the densities of the points		
		/*
		int aux;
		prototype1=-1;
		prototype2=-1;
		for (i=0;i<N;i++){
			if (parent[i]==chosen){
				if (prototype1==-1)
					prototype1=i;
				else if (prototype2==-1){			
					prototype2=i;
					if (rho[prototype2]>rho[prototype1]){
						aux=prototype2;
						prototype2=prototype1;
						prototype1=aux;					
					}
				}
				else{
					if (rho[i]>=rho[prototype1]){
						prototype2=prototype1;
						prototype1=i;
					}
					else if (rho[i]>rho[prototype2])
						prototype2=i;
				}
				
			}
		}*/
		// Randomly finding two prototypes in the chosen cluster according to a probability proportional to the density
	
		p_chosen=aloc_vectord(N);
		auxd=0.0;
		for (i=0;i<N;i++){
			p_chosen[i]=0.0;
			if (parent[i]==chosen){
				p_chosen[i]=rho[i];
				auxd=auxd+p_chosen[i];
			}
		}
		if (auxd==0.0){
			delete [] p_chosen;	
			delete [] size_clust;
			return;
		}
		
		for (i=0;i<N;i++)
			p_chosen[i]=p_chosen[i]/auxd;
		r1=random_dou();
		r2=random_dou();
		if (r1>r2){
			auxd=r1;
			r1=r2;
			r2=auxd;
		}
		auxd=0.0;
		prototype1=0;
		prototype2=0;

		for (i=0;i<N;i++){
			auxd=auxd+p_chosen[i];
			if (auxd>=r1){
				prototype1=i;
				r1=2.0;				
			}
			if (auxd>=r2){
				prototype2=i;
				r2=2.0;				
			}	
		}

		// Finding label for the new cluster
		new_label=0;
		l=1;
		while (l<=label_max && new_label==0){	
			if (size_clust[l]==0)
				new_label=l;
			l++;		
		}
		if (new_label==0)
			new_label=label_max+1;

		// Labeling points to the new cluster
		for (i=0;i<N;i++)
			if (parent[i]==chosen)
				if (D[i][prototype2]<D[i][prototype1])
					offspring[i]=new_label;
		
		delete [] p_chosen;			
	}

	delete [] size_clust;
}

/******************************************************************************\
*								Print NK information						   *
\******************************************************************************/
void nk::print(void){
	
	cout<< "NK Landscape: "<<endl<<"  N="<<N<<endl<<"  K: "<< K<<endl;			
	cout<<"Epistasis Graph: "<<endl;
	G_ep->imprime ();
	cout<<"Transpose of the Epistasis Graph: "<<endl;
	G_ep_transp->imprime();
	
}


/******************************************************************************\
*			Save  NK information (transpose of the interaction graph)		   *
\******************************************************************************/
void nk::save(char *prob_name){
	
	int i, j, flag_null;
	FILE *NK_file;
	char *name_p;
	char name[CHAR_LEN];

   name_p = name;

	// Best solution
	//sprintf(name,"data/NK_%s.dat",prob_name);
	sprintf(name,"NK_%s.dat",prob_name);
	if ((NK_file = fopen(name_p,"w"))==NULL) {
		cout<<"The file NK to be saved cannot be open"<<endl;
		exit(1);
	}
	for (i=0;i<N;i++){	
		fprintf(NK_file,"%d ",i);
		flag_null=0;		
		Grafo::Aresta *adj = G_ep_transp->primeiroListaAdj(i);	  // first edge to node i
		if (adj == NULL)
			flag_null=1;
		else			
			j=adj->_v2 ();		
		while (flag_null==0) {
			fprintf(NK_file,"%d ",j);
			delete adj;
			adj= G_ep_transp->proxAdj(i);	  // next edge o node i
			if (adj == NULL)
				flag_null=1;
			else			
				j=adj->_v2 ();	
		}	
		fprintf(NK_file,"\n");
	}
	
	
	fclose(NK_file);
	
}

