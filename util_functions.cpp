/******************************************************************************\
*								Diverse Functions						 *
\******************************************************************************/
#include "defs.h"
#include <cstdlib>
#include <cmath>


/******************************************************************************\
*								 Dynamic Allocation: Matrix of Integers					 *
\******************************************************************************/
int **aloc_matrixi(int lines , int collums)
{
	int i, **Matrix;
	
	Matrix = new int*[lines];
	for (i=0;i<lines;i++) {
		Matrix[i] = new int[collums];
	}
	if (!Matrix) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}

	return Matrix;
}

/******************************************************************************\
*								 Dynamic Allocation: Matrix of Doubles					 *
\******************************************************************************/
double **aloc_matrixd(int lines , int collums)
{
	int i;
	double **Matrix;
	
	Matrix = new double*[lines];
	for (i=0;i<lines;i++) {
		Matrix[i] = new double[collums];
	}
	if (!Matrix) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}

	return Matrix;
}


/******************************************************************************\
*								Dynamic Allocation: Vector of Integers						 *
\******************************************************************************/
int *aloc_vectori(int lines)
{
	int *vector;

	vector = new int[lines];
	if (!vector) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}
	return vector;
}
/******************************************************************************\
*								Dynamic Allocation: Vector of Doubles						 *
\******************************************************************************/
double *aloc_vectord(int lines)
{
	double *vector;

	vector = new double[lines];
	if (!vector) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}
	return vector;
}


/******************************************************************************\
*								Dynamic Allocation: Vector of individuals						 *
\******************************************************************************/
individual *aloc_vectorind(int lines)
{
	individual *vector;

	vector = new individual[lines];
	if (!vector) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}
	return vector;
}

/******************************************************************************\
*								Dynamic Allocation: Population						 *
\******************************************************************************/
void aloc_pop(population *pop, int n_ind)
{
	int i;
	
	pop->pop_size=n_ind;
	pop->ind = aloc_vectorind(pop->pop_size);
	for (i=0;i<pop->pop_size;i++)
		pop->ind[i].chromosome = aloc_vectori(lcrom);
}

/******************************************************************************\
*								 Dynamic Desallocation: Matrix of Integers					 *
\******************************************************************************/
void desaloc_matrixi(int **Matrix , int lines)
{
	int i;

	for(i=0;i<lines;i++) {
		delete [] Matrix[i];
	}
	delete [] Matrix;
}

/******************************************************************************\
*								 Dynamic Desallocation: Matrix of Doubles				 *
\******************************************************************************/
void desaloc_matrixd(double **Matrix , int lines)
{
	int i;

	for(i=0;i<lines;i++) {
		delete [] Matrix[i];
	}
	delete [] Matrix;
}

/******************************************************************************\
*								 Dynamic Desallocation: Population        	   *
\******************************************************************************/
void desaloc_pop( population *pop )
{
	int i;
	
	for (i=0;i<pop->pop_size;i++)
		delete [] pop->ind[i].chromosome;
		
	delete [] pop->ind;
}

/******************************************************************************\
*								 Random Integer between L_range and H_range	   *
\******************************************************************************/
int random_int(int L_range, int H_range)
{
	return(  (int) ( (rand()/(RAND_MAX+1.0))*(H_range-L_range+1)+L_range ) );  // random integer beteween [L_range and H_range]
}

/******************************************************************************\
*								 Random double in [0.0,1.0]			 		   *
\******************************************************************************/
double random_dou(void)
{
	return(  rand() / double(RAND_MAX) );  //  random double in [0.0, 1.0]:
}

/******************************************************************************\
*		Random Permutation of the vector of integers v				 		   *
\******************************************************************************/
void rand_perm(int *inp, int *out, int size)
{
	int i, j;
	
	out[0]=inp[0];
	for(i=1;i<size;i++) {
		j= random_int(0,i);  
		if (i != j)
			out[i]=out[j];
		out[j]=inp[i];
	}
}

/******************************************************************************\
*								 Quicksort: qsort							 *
\******************************************************************************/
void qsort( double *a , int L, int R)
{
	int i, j;
	double x, w;

	i = L;
	j = R;
	x = a [  (L + R ) / 2 ] ;
	do {
		while ( a[ i ] < x ){
    			i = i + 1;
		}
		while ( x < a[ j ] ){
   			 j = j - 1;
		}
		if ( i <= j ){
  			w = a[ i ]; 
    		a[ i ] = a[ j ];
    		a[ j ] = w; 
    		i = i + 1;
    		j = j - 1;
		}
	} while ( i <= j );
	if ( L < j ){
		qsort( a , L, j );
	}
	if ( i < R ){
		qsort( a , i, R );
	}
	return;
}

/******************************************************************************\
*								 Quicksort: qsort							 *
\******************************************************************************/
void quicksort(double *v, int N){

	qsort( v , 0, N-1 );	
	
}


/******************************************************************************\
*				Squared Euclidean Distance									   *
\******************************************************************************/
double sqrEucDist(double **X, int i, int j, int l){
	int k;
	double sum=0.0;

	for (k=0;k<l;k++)
		sum=sum+pow( (X[i][k]-X[j][k]) , 2 );

	return sum;
}
