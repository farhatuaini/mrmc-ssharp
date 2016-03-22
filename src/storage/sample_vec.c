/**
*	WARNING: Do Not Remove This Section
*
*       $LastChangedRevision: 415 $
*       $LastChangedDate: 2010-12-18 17:21:05 +0100 (Sa, 18. Dez 2010) $
*       $LastChangedBy: davidjansen $
*
*	MRMC is a model checker for discrete-time and continuous-time Markov
*	reward models. It supports reward extensions of PCTL and CSL (PRCTL
*	and CSRL), and allows for the automated verification of properties
*	concerning long-run and instantaneous rewards as well as cumulative
*	rewards.
*
*	Copyright (C) The University of Twente, 2004-2008.
*	Copyright (C) RWTH Aachen, 2008-2009.
*	Authors: Ivan Zapreev, Christina Jansen
*
*	This program is free software; you can redistribute it and/or
*	modify it under the terms of the GNU General Public License
*	as published by the Free Software Foundation; either version 2
*	of the License, or (at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program; if not, write to the Free Software
*	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*	MA  02110-1301, USA.
*
*	Main contact:
*		Lehrstuhl für Informatik 2, RWTH Aachen University
*		Ahornstrasse 55, 52074 Aachen, Germany
*		E-mail: info@mrmc-tool.org
*
*       Old contact:
*		Formal Methods and Tools Group, University of Twente,
*		P.O. Box 217, 7500 AE Enschede, The Netherlands,
*		Phone: +31 53 4893767, Fax: +31 53 4893247,
*		E-mail: mrmc@cs.utwente.nl
*
*	Source description: This is a source file for the sample-vector
*	structure and operations on it. It is used in the simulation engine.
*/

#include "sample_vec.h"
#include "sparse.h"

#include <string.h>

/****************************************************************************/
/************THE ALLOCATION AND DEALLOCATION OF THE SAMPLE VECTOR************/
/****************************************************************************/

/**
* This method is used to initialize the base of every sample vector
* @param pSampleVecBase the structure to initialize
* @param sample_size the initial sample size
* @param initial_state the initial observation
*/
static inline void initializeSampleVectorBase(PTSampleVec pSampleVecBase, const int sample_size, const int initial_state ){
	IF_SAFETY( pSampleVecBase != NULL )
		/* Assign the initial state */
		pSampleVecBase->initial_state = initial_state;

		/* Assign the sample size right away */
		pSampleVecBase->curr_sample_size = sample_size;

		/* Assign the initial number of visited states */
		pSampleVecBase->num_visited_states = sample_size;
	ELSE_SAFETY
		printf("ERROR: Trying to initialize an unallocated base sample.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* This method is used to allocate a new sample-vector structure for the
* unbounded-until simulations CTMC/DTMC:
* We set "curr_sample_size" to "sample_size".
* We set "curr_simulation_depth" to zero.
* We set "sum_good" to zero.
* We set "sum_trans" to "sample_size".
* We allocate "pTransStateInd" and "pObservationsVec" of size "sample_size".
* We set "pTransStateInd" to all ones.
* We fill "pObservationsVec" with "initial_state".
* @param sample_size the initial sample size
* @param initial_state the initial observation.
* @return the newly-created sample vector.
*/
inline PTSampleVecUntil allocateSampleVectorUntil( const int sample_size, const int initial_state ){
        PTSampleVecUntil pSampleVecUntil = (PTSampleVecUntil) calloc((size_t) 1,
                        sizeof(TSampleVecUntil));
	int i = 0;

	IF_SAFETY( pSampleVecUntil != NULL )
		/* Initialize the sample base */
		initializeSampleVectorBase( (PTSampleVec) pSampleVecUntil, sample_size, initial_state );

                /* It is always zero in the beginning */
		pSampleVecUntil->curr_simulation_depth = 0;

		/* Asuming the initial state is transient we */
		/* have that all states in the sample are */
		/* identical and transient */
		pSampleVecUntil->sum_good = 0;
		pSampleVecUntil->sum_trans = sample_size;

                pSampleVecUntil->pObservationsVec = (int*) calloc(
                                (size_t) sample_size, sizeof(int));
		/* Can't use memset here since "initial_state" is an integer */
		for( i = 0; i < sample_size; i++ ){
			pSampleVecUntil->pObservationsVec[i] = initial_state;
		}

		pSampleVecUntil->pTransStateInd = get_new_bitset( sample_size );
		/* Sets all bits, including unused bits of the last block, to one! */
		fill_bitset_one( pSampleVecUntil->pTransStateInd );
	ELSE_SAFETY
		printf("ERROR: We've run out of memory.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY

	return pSampleVecUntil;
}

/**
* This method is used to allocate a new sample-vector structure for the
* interval-until simulations CTMC/DTMC:
* We set "curr_sample_size" to "sample_size".
* We set "sum_good" to zero.
* @param sample_size the initial sample size
* @param initial_state the initial observation.
* @return the newly-created sample vector.
*/
inline PTSampleVecIntUntil allocateSampleVectorIntUntil( const int sample_size, const int initial_state ){
        PTSampleVecIntUntil pSampleVecIntUntil = (PTSampleVecIntUntil) calloc(
                        (size_t) 1, sizeof(TSampleVecIntUntil));

	IF_SAFETY( pSampleVecIntUntil != NULL )
		/* Initialize the sample base */
		initializeSampleVectorBase( (PTSampleVec) pSampleVecIntUntil, sample_size, initial_state );

		/* The good states are originally set to zero */
		pSampleVecIntUntil->sum_good = 0;
	ELSE_SAFETY
		printf("ERROR: We've run out of memory.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY

	return pSampleVecIntUntil;
}

/**
* This method is used to allocate a new sample-vector structure for the
* steady-state simulations CTMC:
* We set "initial_state" to "initial_state".
* We set "curr_reg_cycles" to reg_cycles.
* We set "pRegenerationPoints" to all zeros.
* We set "pTimeObservationsVec" to all zeros.
* We set "pGoodObservationsVec" to all zeros.
* @param reg_cycles the initial number of regeneration cycles.
* @param initial_state the regeneration point (the state index).
* @return the newly-created sample vector.
*/
inline PTSampleVecSteady allocateSampleVectorSteady(const int reg_cycles, const int initial_state){
        PTSampleVecSteady pSampleVecSteady=(PTSampleVecSteady) calloc((size_t)1,
                        sizeof(TSampleVecSteady));

	IF_SAFETY( pSampleVecSteady != NULL )
		/* Initialize the sample base */
		initializeSampleVectorBase( (PTSampleVec) pSampleVecSteady, reg_cycles, initial_state );

		/* As no regeneration cycle was calculated yet, we allocate an array filled with 0 */
                pSampleVecSteady->pTimeObservationsVec = (double *) calloc(
                                (size_t) reg_cycles, sizeof(double));

		/* As no regeneration cycle was calculated yet, we allocate an array filled with 0 */
                pSampleVecSteady->pGoodObservationsVec = (double *) calloc(
                                (size_t) reg_cycles, sizeof(double));
	ELSE_SAFETY
		printf("ERROR: We've run out of memory.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY

	return pSampleVecSteady;
}

/**
* This method is used to extend a sample base vector
* We can only increase the sample size, i.e. we expect
*	"new_sample_size >= pSampleVecBase->curr_sample_size"
* @param pSampleVecBase the sample-vector base of the sample scheduled for reallocation.
* @param sample_size the new sample size
* @param pVectorName the name of the extended vector (the original)
* @return the old sample size
*/
static int extendSampleVectorBase(PTSampleVec pSampleVecBase,
                const int new_sample_size, const char * pVectorName)
{
	int old_sample_size = 0;
	IF_SAFETY( pSampleVecBase != NULL )

		/* Save the old sample size*/
		old_sample_size = pSampleVecBase->curr_sample_size;

		IF_SAFETY( old_sample_size <= new_sample_size )
			/* Reset the sample size*/
			pSampleVecBase->curr_sample_size = new_sample_size;

			/* The newly added initial states shall be counted as visited */
			pSampleVecBase->num_visited_states += (new_sample_size - old_sample_size);
		ELSE_SAFETY
			printf("ERROR: Trying to decrease the sample size of %s.\n", pVectorName);
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: Trying to reallocate a NULL pointer %s.\n", pVectorName);
                exit(EXIT_FAILURE);
	ENDIF_SAFETY

	return old_sample_size;
}

/**
* This method is used to extend a sample-vector structure for the
* unbounded-until simulations CTMC/DTMC:
* We can only increase the sample size, i.e. we expect
*	"new_sample_size >= pSampleVecUntil->curr_sample_size"
* @param pSampleVecUntil the sample-vector scheduled for reallocation.
* @param sample_size the new sample size
* @return the newly-created sample vector.
*/
inline void extendSampleVectorUntil( PTSampleVecUntil pSampleVecUntil, const int new_sample_size ){
	int i = 0;
	bitset * temp_pTransStateInd;

	/* Cast to the parent structure which contains initial_state and curr_sample_size */
	PTSampleVec pSampleVecUntilBase = (PTSampleVec) pSampleVecUntil;

	/* NOTE: Safety assertions are done inside extendSampleVectorBase  */
	const int old_sample_size = extendSampleVectorBase( pSampleVecUntilBase, new_sample_size, "PTSampleVecUntil" );
	const int initial_state = pSampleVecUntilBase->initial_state;

	/* Add extra observations */
	pSampleVecUntil->pObservationsVec = ( int* ) realloc( pSampleVecUntil->pObservationsVec,
								new_sample_size * sizeof( int ) );
	for( i = old_sample_size; i < new_sample_size ; i++){
		pSampleVecUntil->pObservationsVec[i] = initial_state;
	}

	/* Add the new transient states */
	pSampleVecUntil->sum_trans += new_sample_size - old_sample_size;
	temp_pTransStateInd = bs_extend(pSampleVecUntil->pTransStateInd,
			new_sample_size, TRUE);
	if ( (bitset *) NULL == temp_pTransStateInd ) {
                exit(err_macro_2(err_MEMORY, "extendSampleVectorUntil(%p,%d)",
				(void *) pSampleVecUntil, new_sample_size,
                                EXIT_FAILURE));
	}
	pSampleVecUntil->pTransStateInd = temp_pTransStateInd;
}

/**
* This method is used to extend a sample-vector structure for the
* interval-until simulations CTMC/DTMC:
* We can only increase the sample size, i.e. we expect
*	"new_sample_size >= pSampleVecUntil->curr_sample_size"
* @param pSampleVecIntUntil the sample-vector scheduled for reallocation.
* @param sample_size the new sample size
* @return the newly-created sample vector.
*/
inline void extendSampleVectorIntUntil( PTSampleVecIntUntil pSampleVecIntUntil, const int new_sample_size ){
	extendSampleVectorBase( (PTSampleVec) pSampleVecIntUntil, new_sample_size, "PTSampleVecIntUntil");
}

/**
* This method is used to extend a sample-vector structure for the
* steady-state simulations CTMC:
* We can only increase the regeneration cycles, i.e. we expect
*	"new_reg_cycles >= pSampleVecSteady->curr_reg_cycles"
* @param pSampleVecSteady the sample-vector scheduled for reallocation.
* @param new_reg_cycles the new number of regeneration cycles.
* @return the newly-created sample vector.
*/
inline void extendSampleVectorSteady( PTSampleVecSteady pSampleVecSteady, const int new_reg_cycles ){
	/* Cast to the parent structure which contains initial_state and curr_sample_size */
	PTSampleVec pSampleVecSteadyBase = (PTSampleVec) pSampleVecSteady;

	/* NOTE: Safety assertions are done inside extendSampleVectorBase  */
	const int old_reg_cycles = extendSampleVectorBase( pSampleVecSteadyBase, new_reg_cycles, "PTSampleVecSteady" );

	/* Add extra regeneration cycles */
	pSampleVecSteady->pTimeObservationsVec = ( double* ) realloc( pSampleVecSteady->pTimeObservationsVec,
								new_reg_cycles * sizeof( double ) );
	/* Set the newly added array elements to zero */
	memset( &pSampleVecSteady->pTimeObservationsVec[old_reg_cycles], 0, (new_reg_cycles - old_reg_cycles) * sizeof( double ) );

	pSampleVecSteady->pGoodObservationsVec = ( double* ) realloc( pSampleVecSteady->pGoodObservationsVec,
								new_reg_cycles * sizeof( double ) );
	/* Set the newly added array elements to zero */
	memset( &pSampleVecSteady->pGoodObservationsVec[old_reg_cycles], 0, (new_reg_cycles - old_reg_cycles) * sizeof( double ) );
}

/****************************************************************************/
/***************FREE THE MEMORY ALLOCATED FOR THE SAMPLE VECTOR**************/
/****************************************************************************/

/**
* This method is used to free the sample-vector structure for the
* unbounded-until simulations CTMC/DTMC:
* @param pSampleVecUntil the sample vector that has to be freed.
*/
inline void freeSampleVectorUntil( PTSampleVecUntil pSampleVecUntil ){
	if( pSampleVecUntil != NULL ){
		if( pSampleVecUntil->pObservationsVec != NULL ){
			free( pSampleVecUntil->pObservationsVec );
		}
		if( pSampleVecUntil->pTransStateInd != NULL ){
			free_bitset( pSampleVecUntil->pTransStateInd );
		}
		free( pSampleVecUntil );
	}
}

/**
* This method is used to free the sample-vector structure for the
* interval-until simulations CTMC/DTMC:
* @param pSampleVecIntUntil the sample vector that has to be freed.
*/
inline void freeSampleVectorIntUntil( PTSampleVecIntUntil pSampleVecIntUntil ){
	if( pSampleVecIntUntil != NULL ){
		free( pSampleVecIntUntil );
	}
}

/**
* This method is used to free the sample-vector structure for the steady-state
* simulations CTMC:
* @param pSampleVecSteady the sample vector that has to be freed.
*/
inline void freeSampleVectorSteady ( PTSampleVecSteady pSampleVecSteady ){
	if (pSampleVecSteady != NULL ){
		if( pSampleVecSteady->pTimeObservationsVec != NULL ){
			free( pSampleVecSteady->pTimeObservationsVec );
		}
		if( pSampleVecSteady->pGoodObservationsVec != NULL ){
			free( pSampleVecSteady->pGoodObservationsVec );
		}
		free( pSampleVecSteady);
	}
}

/****************************************************************************/
/***************************PRINT THE SAMPLE VECTOR**************************/
/****************************************************************************/

/**
* This method prints the sample-vector base structure for all other samples
* @param pSampleVecBase the sample-vector base
* WARNING: There is no check for pSampleVecBase == NULL
*/
static void printSampleVector(PTSampleVec pSampleVecBase,
                const char * pVectorName)
{
	printf( "------- Printing %s -------\n", pVectorName );
	printf( "\tinitial_state = %d\n", pSampleVecBase->initial_state + 1 );
	printf( "\tcurr_sample_size = %d\n", pSampleVecBase->curr_sample_size );
	printf( "\tnum_visited_states = %u\n", pSampleVecBase->num_visited_states );
}

/**
* This method prints the sample-vector structure for unbounded-until simulation
* @param pSampleVecUntil the sample-vector to be printed
*/
inline void printSampleVectorUntil( PTSampleVecUntil pSampleVecUntil ){
	PTSampleVec pSampleVecUntilBase = (PTSampleVec) pSampleVecUntil;
	if( pSampleVecUntil != NULL ){
		printSampleVector( pSampleVecUntilBase, "PTSampleVecUntil");

		printf( "\tcurr_simulation_depth = %d\n", pSampleVecUntil->curr_simulation_depth );
		printf( "\tsum_good = %d\n", pSampleVecUntil->sum_good );
		printf( "\tsum_trans = %d\n", pSampleVecUntil->sum_trans );
		printf( "\tpObservationsVec :\n\t\t" );
                if ( err_state_iserror(print_vec_int(
                                        pSampleVecUntilBase->curr_sample_size,
                                        pSampleVecUntil->pObservationsVec)) )
                {
                        exit(err_macro_1(err_CALLBY, "printSampleVectorUntil("
                                "%p)", (void *) pSampleVecUntil, EXIT_FAILURE));
                }
		printf( "\tWARNING: For getting an external (a user-level) state index do +1.\n" );
		printf( "\tpTransStateInd :\n\t\t" );
		print_bitset_states( pSampleVecUntil->pTransStateInd );
		printf( "\n" );
	} else {
		printf( "\tWARNING: Trying to print the NULL pointer of type PTSampleVecUntil.\n" );
	}
}

/**
* This method prints the sample-vector structure for interval-until simulation CTMC/DTMC
* @param pSampleVecIntUntil the sample-vector to be printed
*/
inline void printSampleVectorIntUntil( PTSampleVecIntUntil pSampleVecIntUntil ){
	if( pSampleVecIntUntil != NULL ){
		printSampleVector( (PTSampleVec) pSampleVecIntUntil, "PTSampleVecIntUntil");

		printf( "\tsum_good = %d\n", pSampleVecIntUntil->sum_good );
		printf( "\n" );
	} else {
		printf( "\tWARNING: Trying to print the NULL pointer of type PTSampleVecIntUntil.\n" );
	}
}

/**
* This method prints the sample-vector structure for steady-state simulation CTMC
* @param pSampleVecSteady the sample-vector to be printed
*/
inline void printSampleVectorSteady( PTSampleVecSteady pSampleVecSteady ){
	PTSampleVec pSampleVecSteadyBase = (PTSampleVec) pSampleVecSteady;
	if( pSampleVecSteady != NULL ){
		printSampleVector( (PTSampleVec) pSampleVecSteady, "PTSampleVecSteady");

                if ( err_state_iserror(print_vec_double(
                                pSampleVecSteadyBase->curr_sample_size,
                                pSampleVecSteady->pTimeObservationsVec)) )
                {
                        exit(err_macro_1(err_CALLBY, "printSampleVectorSteady("
                                "%p)", (void*) pSampleVecSteady, EXIT_FAILURE));
                }
		printf( "\n" );
                if ( err_state_iserror(print_vec_double(
                                pSampleVecSteadyBase->curr_sample_size,
                                pSampleVecSteady->pGoodObservationsVec)) )
                {
                        exit(err_macro_1(err_CALLBY, "printSampleVectorSteady("
                                "%p)", (void*) pSampleVecSteady, EXIT_FAILURE));
                }
		printf( "\n\n" );
	} else {
		printf( "\tWARNING: Trying to print the NULL pointer of type PTSampleVecSteady.\n" );
	}
}
