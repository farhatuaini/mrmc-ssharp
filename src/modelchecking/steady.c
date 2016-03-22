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
*	Authors: Ivan Zapreev
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
*	Source description: This file contains sources for the steady
*		state analysis of the PCTL formulas for CTMCs.
*/

#include "steady.h"

#include "bscc.h"
#include "transient.h"
#include "iterative_solvers.h"

#include "runtime.h"

/* This part is required in order to prevent reuse of probabilities precomputed cached */
/* for the case MRMC runtime settings are re-set. */
#define UNDEFINED -1  /*Undefined*/
/**
* This method is used in order to check whether the runtime settings
* that might influence the computation of the steady-sate probabilities have changed.
* It also stores the new values in case of such a change.
* @return This method returns FALSE if the error bound or the
*         maximum number of iterations or the iteration method have NOT changed.
*         Otherwise, i.e. if at least one of them has changed, it returns TRUE,
*         unless the last_* values are UNDEFINED.
*         Then this is the first run and nothing has been computed yet, so
*         the method returns FALSE.
*/
static BOOL isSettingsChange(void) {
        /* Store the previously used error bound */
        static double last_err = UNDEFINED;
        /* Store the previously used maximum number of iterations */
        static int last_max_iterations = UNDEFINED;
        /* Store the previously used iteration method */
        static int last_method = UNDEFINED;

	double err = get_error_bound();
	int max_iterations = get_max_iterations();
	int method = get_method_steady();
	BOOL isDefined = ( last_method != UNDEFINED ) || ( last_max_iterations != UNDEFINED ) || ( last_err != UNDEFINED );
	BOOL hasChanged = ( last_method != method ) || ( last_max_iterations != max_iterations ) || ( last_err != err );

        if ( hasChanged ) {
		last_method = method;
		last_max_iterations = max_iterations;
		last_err = err;
	}

	return isDefined && hasChanged;
}

/*This variable stores all steady state probabilities for the resolved bsccs.*/
static double * pSteadyStateProbs = NULL;
/*The First time call marker*/
static BOOL isFirstTime = TRUE;
/*This variable contains the state space i.e. the initial Rate matrix*/
static sparse * pStateSpace = NULL;
/*This pointer stores the Q matrix pointer*/
static sparse * pQ = NULL;
/*This stores the number of states in the state space*/
static int N_STATES = 0;
/*Stores E values for the ergodic matrix*/
static double * pErgodicE = NULL;
/*The true bitset*/
static bitset * pbTRUEBitSet = NULL;
/*The marker of ergodic CTMC*/
static BOOL isErgodicCTMC = FALSE;
/*This variable contains the BSCC search data for steady-state operator and initial matrix*/
static TBSCCs * pBSCCsHolder = NULL;
/**
* This method returns the Q matrix for the BSCC (Q=R-diag(E))
* Where R is the rate matrix of the BSCC
* @param pBSCC the BSCC nodes of the original matrix
*/
static void computeQMatrix(sparse * pM, const int * pValidStates)
{
	if( pValidStates )
	{
		/*Get the E vector*/
		double * pE = get_mtx_cer_row_sums(pM, pValidStates);
                if ( NULL == pE
                /*Get the Q matrix ( Q=R-diag(E) )*/
                                || err_state_iserror(sub_mtx_cer_diagonal(pM,
                                                pValidStates, pE)) )
                {
                        exit(err_macro_4(err_CALLBY, "computeQMatrix(%p[%dx%d],"
                                "%p)", (void *) pM, mtx_rows(pM), mtx_cols(pM),
                                (const void *) pValidStates, EXIT_FAILURE));
                }
		/*Clean the Row Sum array*/
		free(pE);
	}
	else
	{
		/*If we did not compute E vector for this matrix then get the E vector*/
		pErgodicE = get_mtx_row_sums(pM);
                if ( NULL == pErgodicE
                /*Get the Q matrix ( Q=R-diag(E) )*/
                                || err_state_iserror(sub_mtx_diagonal(pM,
                                                                pErgodicE)) )
                {
                        exit(err_macro_4(err_CALLBY, "computeQMatrix(%p[%dx%d],"
                                "%p)", (void *) pM, mtx_rows(pM), mtx_cols(pM),
                                (const void *) pValidStates, EXIT_FAILURE));
                }
	}
}

/**
* This method solves BSCC in order to get steady state probabilities
* @param pM the matrix we need to solve
* @param pValidStates the valid states, pValidStates[0] is the number of
*                     states can be NULL if all states are valid
*/
static void obtainSteadyStateProbabilities(sparse * pM, int * pValidStates)
{
	double * pB = NULL;
	double * pResult = NULL;
	int i, length;

	/* NOTE: Here we store the current error_bound, max_iterations and method settings */
	/* In order to avoid reuse of cached probabilities computed with different parameters later. */

	double err = get_error_bound(); /* Retrieve the Error bound from the runtime.c */
	int max_iterations = get_max_iterations(); /* Retrieve the Max number of iterations */
											   /* from the runtime.c */
	int method;
	if( get_method_steady() == GJ) /* Retrieve the l.e. solution method from the runtime.c */
		method = GAUSS_JACOBI_INV;
	else /* This should be "GS" otherwise */
		method = GAUSS_SEIDEL_INV;

	pResult = solve_initial(method, pM, NULL, pB, err, max_iterations, pValidStates);

	/*Normalize the solution*/
	normalizeSolution( N_STATES, pResult, pValidStates, TRUE );
	/* print_vec_double(N_STATES, pResult); */

	/*Store results*/
	if( pValidStates ){
		/*If the CTMC is not ergodic*/
		length = pValidStates[0];
		/*Copy colution to the solutions vector*/
		for( i = 1; i <= length; i++)
		{
			pSteadyStateProbs[pValidStates[i]] = pResult[pValidStates[i]];
		}
		/*Free result*/
		free( pResult );
	}
	else
	{
		/*Free the old results*/
		free( pSteadyStateProbs );
		/*For ergodic CTMC the steady state probabilities computed
		Copy colution to the solutions vector*/
		pSteadyStateProbs = pResult;
	}
}

/**
* This method returns the states that belong to a BSCC
* @param pBSCCInfo the BSCC info a 2 elements array, first element if a BSCC id
*		   the second one is the number of states in the BSCC
* @return the array of states belonging to the given BSCC, the first value
*         of the array is the number of states
*/
static int * initValidStates(const int * pBSCCInfo)
{
	/*Allocate memory for Valid states storage*/
        int * pValidStates = (int *) calloc((size_t) (pBSCCInfo[1] + 1),
                        sizeof(int));
	int theBSCCId = pBSCCInfo[0], i, counter = 1;
        const
	int *pStatesBSCCMapping = getStatesBSCCsMapping(pBSCCsHolder);
	/*Store the number of valid states (the states opf BSCC)*/
	pValidStates[0] = pBSCCInfo[1];

	/*Get valid states;*/
	for( i=0 ; i < N_STATES ; i++)
	{
		if( pStatesBSCCMapping[i] == theBSCCId )
		{
			pValidStates[counter] = i;
			counter++;
		}
	}

	return pValidStates;
}

/**
* Restores the pStateSpace by adding pErgodicE
* to diagonal, and frees the pErgodicE vector.
*/
static void restoreStateSpace(sparse * pStateSpace_local) {
        if ( err_state_iserror(add_mtx_diagonal(pStateSpace_local, pErgodicE)) )
        {
                exit(err_macro_3(err_CALLBY, "restoreStateSpace(%p[%dx%d])",
                        (void *) pStateSpace_local, mtx_rows(pStateSpace_local),
                        mtx_cols(pStateSpace_local), EXIT_FAILURE));
        }
	free( pErgodicE );
	pErgodicE = NULL;
}

/**
* This method computes the probability to reach the newly found bscc
* Compute the probabilities P(true U BSCC)
* @return the probabilities
*/
static double * getReachProbability(const int bscc_id)
{
	bitset * pBSCC = get_new_bitset( N_STATES );
	int i;
        const
	int * bscc_mapping = getStatesBSCCsMapping(pBSCCsHolder);
	double *pUntilResults=NULL;
	/* Get BSCC */
	for( i = 0; i < N_STATES ; i++ ) {
		if( bscc_mapping[ i ] == bscc_id ) {
			set_bit_val( pBSCC, i, BIT_ON );
		}
	}

	/* Compute the probabilities P(true U BSCC) */
	if( isRunMode(CTMC_MODE) || isRunMode(DTMC_MODE) || isRunMode(DMRM_MODE) ){
		pUntilResults = until( TIME_UNBOUNDED_FORM, pbTRUEBitSet, pBSCC, 0.0, 0.0, TRUE );
                if ( NULL == pUntilResults ) {
                        err_msg_1(err_CALLBY, "getReachProbability(%d)",bscc_id,
                                (free_bitset(pBSCC),free(pUntilResults), NULL));
                }
	}

	/*Free memory*/
	free_bitset(pBSCC);

	return pUntilResults;
}

/**
* This method is used to obtain the steady state probabilities of the found BSCCs
* @param ppNewBSCCs the list of newly found BSCCs and its lengths
*/
static void solveBSCCs(const int * const * ppNewBSCCs)
{
	int i;
	int num = getBSCCCount(ppNewBSCCs);
        const
	int * pBSCCInfo = NULL;
	int * pValidStates = NULL;
	/*Iterate through the BSCCs*/
	for( i=1; i <= num; i++ )
	{
		pBSCCInfo = ppNewBSCCs[i];
		/* Check that this is not a trivial - 1 node BSCC */
		if( pBSCCInfo[1] != 1 )
		{
			/*This code works for CSL, PCTL and PRCTL because of the
			  computeQMatrix implementation, i.e. for PCTL and PRCTL
			  it computes Q = P-I*/
			if( pBSCCInfo[1] != N_STATES )
			{
				pValidStates = initValidStates(pBSCCInfo);
                                if ( err_state_iserror(initMatrix(pStateSpace,
                                                        pQ, pValidStates))
                                        || (computeQMatrix(pQ, pValidStates),
                                                obtainSteadyStateProbabilities(
                                                        pQ, pValidStates),
                                                err_state_iserror(cleanMatrix(
                                                        pQ, pValidStates))) )
                                {
                                        exit(err_macro_1(err_CALLBY,
                                                "solveBSCCs(%p)",
                                                (const void *) ppNewBSCCs,
                                                EXIT_FAILURE));
                                }
				free(pValidStates);
			}
			else
			{
				/*The ergodic CTMC*/
				printf("The given Process is ergodic.\n");
				isErgodicCTMC = TRUE;
				computeQMatrix( pStateSpace, NULL );
				obtainSteadyStateProbabilities( pStateSpace, NULL );
				restoreStateSpace(pStateSpace);
			}
		}
		else
		{
			/*Store 1.0 probability for a 1 node BSCC*/
			pSteadyStateProbs[ pBSCCInfo[2] ] = 1.0;
		}
	}
        /* printf("Steady State Probabilities :\n"); */
	/* print_vec_double( N_STATES, pSteadyStateProbs ); */
}

/**
* This method comutes the S(F) probabilities taking into account the steady state probabilities
* ToDo: When we know all BSCCs then we know the number of transient states and thus can store
* the probabilities to reach BSCC for these transient states in memory.
* @param the F states
* @return the array of probabilities
*/
static double * getS_F_Solution(const bitset * pStates)
{
	/*Allocate final results storage*/
        double * pResult = (double *) calloc((size_t) N_STATES, sizeof(double));
	/*Allocate storage for the sum of steady state probabilities */
	int BSCC_NUM = getBSCCCounter(pBSCCsHolder)+1; /*This +1 allows to remove bscc_id-1 expressions*/
        double * pPiBSCCProbs = (double *) calloc((size_t) BSCC_NUM,
                        sizeof(double));
	/*Define other local variables*/
	double * pProbToReach = NULL;
	int i, j, bscc_id = 0;
        const
	int * bscc_mapping = getStatesBSCCsMapping(pBSCCsHolder);
	bitset *pUsedBSCCs = NULL;

	/*Compute steady state probabilities sums for BSCC intersect Set(F)*/
        /* get_idx_next_non_zero() is more efficient than checking every bit in
           pStates individually. David N. Jansen. */
        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(pStates, i)) != state_index_NONE ) {
		/*If BSCC intersect Set(F) is not 0 then*/
                if ( (bscc_id = bscc_mapping[i]) != 0 ) {
			pPiBSCCProbs[ bscc_id ] += pSteadyStateProbs[ i ];
		}
	}

	/* if (bscc_id != 0) printf("SS-Probablitity for bscc  = %lf\n", pPiBSCCProbs[ bscc_id ]); */

	/*Look through all known bsccs*/
	if( ! isErgodicCTMC ) {
		pUsedBSCCs = get_new_bitset( BSCC_NUM );

		printf("Getting the probability of reaching BSCCs ...\n");

                i = state_index_NONE;
                while ( (i = get_idx_next_non_zero(pStates, i))
                                        != state_index_NONE )
                {
			/*If BSCC intersect Set(F) is not 0 then*/
                        if ( (bscc_id = bscc_mapping[i]) != 0 ) {
				if( ! get_bit_val( pUsedBSCCs, bscc_id ) ) {
					/*Get probability to reach the given BSCC*/
					pProbToReach = getReachProbability( bscc_id );
                                        if ( NULL == pProbToReach ) {
                                                err_msg_2(err_CALLBY,
                                                        "getS_F_Solution("
                                                        "%p[%d])",
                                                        (const void *) pStates,
                                                        bitset_size(pStates),
                                                        NULL);
                                        }

                                        /*printf("Probabilities to reach BSCC "
                                                        "%d:\n", bscc_id);
					print_vec_double( N_STATES, pProbToReach );*/

					/*Compute final values*/
					for( j = 0; j < N_STATES; j++) {
						pResult[j] += pPiBSCCProbs[ bscc_id ] * pProbToReach[j];
					}
					/*Free memory, later we will store this data for reuse*/
					free( pProbToReach );
					/*Mark BSCC as used in summation*/
					set_bit_val( pUsedBSCCs, bscc_id, BIT_ON);
				}
			}
		}
		free_bitset( pUsedBSCCs );
	} else {
		/*If the CTMC is ergodic (there is only one BSCC) then
		the probability is equal to the steady state probability
		to in the set(F) states*/
		for( j = 0; j < N_STATES; j++) {
			pResult[j] = pPiBSCCProbs[ 1 ]; /*Here 1 is because the BSCC is only one*/
		}
	}

	free(pPiBSCCProbs);
	return pResult;
}

/**
* This method gets the probabilities for the S(F) formula for every state.
* @param pStates the pointer to the bit set that indicates the set(F) states.
* @return the array of probabilities. For each state there is a steady state
*               probability value.
*/
double* steady(const bitset * pStates)
{
	double * pResult = NULL;
	int **ppNewBSCCs = NULL;

	/* If the MRMC runtime settings have changed we */
	/* have to recompute the probabilities for all BSCCs. */
	if( isSettingsChange() ){
		/* After this method call we should run as if it */
		/* is the first time the steady(...) is called. */
		freeSteady();
	}

	/*Do first time initializations*/
	if(isFirstTime)
	{
		/*Store the state space*/
		pStateSpace = get_state_space();
		/*Get the number of states*/
                N_STATES = bitset_size(pStates);
		/*Allocate the theSteadyStateProbs*/
                pSteadyStateProbs = (double *) calloc((size_t) N_STATES,
                                sizeof(double));
		/*Init matrix*/
		pQ = allocate_sparse_matrix(N_STATES, N_STATES);
                if ( NULL == pQ ) {
                        err_msg_2(err_CALLBY, "steady(%p[%d])",
                                (const void *) pStates, bitset_size(pStates),
                                NULL);
                }
		/*Create True bitset*/
		pbTRUEBitSet = get_new_bitset( N_STATES );
		fill_bitset_one( pbTRUEBitSet );

		isFirstTime = FALSE;
		isErgodicCTMC = FALSE;

		/*Create the initial structure for storing BSCC search data*/
		pBSCCsHolder = allocateTBSCCs(pStateSpace);
	}
	/*Find the BSCCs*/
	printf("Find new BSCCs....\n");
	ppNewBSCCs = getNewBSCCs(pBSCCsHolder , pStates );

	/*Solve the system of linear equations for all newly
	found BSCC to find the steady state probabilities*/
	printf("Solve new BSCCs....\n");
        solveBSCCs((const int **) ppNewBSCCs);

	freeBSCCs(ppNewBSCCs);

	/*Compute the S(F) probabilities using steady state probabilities*/
	printf("Compute S(F)....\n");
	pResult = getS_F_Solution(pStates);
        if ( NULL == pResult ) {
                err_msg_2(err_CALLBY, "steady(%p[%d])", (const void *) pStates,
                        bitset_size(pStates), NULL);
        }

        /* printf("Result :\n"); */
	/* print_vec_double( N_STATES, pResult ); */

	return pResult;
}

/**
* This method is used to reset the steady state analysis before
* model checking the new matrix.
*/
void freeSteady(void)
{
	freeTBSCC(pBSCCsHolder);
	pBSCCsHolder = NULL;

	if( pSteadyStateProbs ){
		free( pSteadyStateProbs );
		pSteadyStateProbs = NULL;
	}

	isErgodicCTMC = FALSE;
	pStateSpace = NULL;
	N_STATES = 0;

	/*Clean up the Q matrix*/
        if ( NULL != pQ ) {
                if ( err_state_iserror(free_mtx_wrap(pQ)) ) {
                        exit(err_macro_0(err_CALLBY, "freeSteady()",
                                        EXIT_FAILURE));
                }
        }
	pQ = NULL;

	if( pbTRUEBitSet ) free_bitset( pbTRUEBitSet );
	pbTRUEBitSet = NULL;

	isFirstTime = TRUE;
}
