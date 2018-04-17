/******************************************************************************\
*								 Selection									 *
\******************************************************************************/
#include "defs.h"
#include <cstdlib>


/******************************************************************************\
*								Tournament Selection							 *
\*******************************************************************************/
int selection_tournament(population *pop)
{
	int i, individual_rand, individual_chosen;
	
	individual_chosen=random_int(0, (pop->pop_size-1) );
	for (i=1;i<tournament_size;i++){
		individual_rand=random_int(0, (pop->pop_size-1) );
		if ( pop->ind[individual_rand].fitness < pop->ind[individual_chosen].fitness  )
			individual_chosen=individual_rand;
	}
		
	return individual_chosen;
}


/******************************************************************************\
*								Selection										 *
\*******************************************************************************/
int selection( population *pop ) 
{
	int individual_chosen;
	
	individual_chosen = selection_tournament( pop );

 	return individual_chosen;
}



