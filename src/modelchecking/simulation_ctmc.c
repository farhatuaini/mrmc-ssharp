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
*	Source description: This is a header file for the simulation engine
*	here we intend to define the CTMC-related functions and data
*	structures.
*/

#include "simulation_ctmc.h"

#include "simulation_common.h"
#include "simulation.h"
#include "transient_common.h"

#include <gsl/gsl_cdf.h>
#include <math.h>

/*This macro computes the N'th root of X */
#define ROOT(X, N) (pow((X), 1.0 / (N)))

/**********************************************************************/
/********THE FUNCTIONS FOR SIMULATING UNBOUNDED UNTIL CTMC*************/
/**********************************************************************/

/**
* This method updates the confidence interval borders for states with
* trivially computed probability 1.0. It also allocates the conf.
* int. arrays, in case it was not done beforehand, i.e. the complete
* initial state set consists only of states, whose probabilities can
* be computed trivially.
* NOTE: pAlwaysBitSet is allowed to be NULL!
* @param isSimOneInitState TRUE if we simulate one initial state only
* @param ppProbCILeftBorder for storing the left conf. int. border
* @param ppProbCIRightBorder for storing the right conf. int. border
* @param pAlwaysBitSet the states that satisfy the formula with probability >= 1.0
* @param initial_state the initial state index if isSimOneInitState == TRUE
*/
static void updateConfIntAlwaysStates(const BOOL isSimOneInitState_local,
                double ** ppProbCILeftBorder,
		double ** ppProbCIRightBorder, int * pResultSize, const bitset * pAlwaysBitSet, const int initial_state)
{
	int curr_index, state_space_size;
        if ( isSimOneInitState_local ) {
		/* Is the conf. int. array allocated yet? */
		if( ( *ppProbCILeftBorder ) == NULL ){
			*pResultSize = 1;
                        *ppProbCILeftBorder = (double *) calloc((size_t) 1,
                                        sizeof(double));
                        *ppProbCIRightBorder = (double *) calloc((size_t) 1,
                                        sizeof(double));
		}
		/* If we simulate just one initial state, test if it belongs to pAlwaysBitSet */
		if( ( pAlwaysBitSet != NULL ) && get_bit_val( pAlwaysBitSet, initial_state ) ){
			(*ppProbCILeftBorder)[0] = (*ppProbCIRightBorder)[0] = 1.0;
		}
	} else {
		/* Is the conf. int. array allocated yet? */
		if( ( *ppProbCILeftBorder ) == NULL ){
                        state_space_size = bitset_size(pAlwaysBitSet);
			*pResultSize = state_space_size;
                        *ppProbCILeftBorder = (double *) calloc(
                                (size_t) state_space_size, sizeof(double));
                        *ppProbCIRightBorder = (double *) calloc(
                                (size_t) state_space_size, sizeof(double));
		}
		/* All pAlwaysBitSet states satisfy the formula with probability >= 1.0 */
		if( pAlwaysBitSet != NULL ){
			/* Set the conf. int. borders to 1.0 for every state, which trivially */
			/* fulfills the formula. */
			curr_index = -1;
			while( ( curr_index =  get_idx_next_non_zero( pAlwaysBitSet, curr_index ) ) != -1 ){
				(*ppProbCILeftBorder)[curr_index]  = (*ppProbCIRightBorder)[curr_index] = 1.0;
			}
		}
	}
}

/**
* This function is an implementation of the unboundedUntil function from the PhD Thesis
* of Ivan S. Zapreev. The function does the model checking of unbounded until operator
* on the provided CTMC and uses the runtime-simulation settings from the "simulation.h"
* @param pStateSpace the sparse matrix of the embedded DTMC with good and bad
*			states made absorbing
* @param confidence this is is confidence we should use for simulations.
*			Note that since the formulas can be nested one can not use the
*			overall confidence set in simulation.h
* @param pPhiBitSet the bitsets containing all the good absorbing states
* @param pPsiBitSet the bitsets containing all the transient states
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param isSimOneInitState TRUE if we simulate one initial state only
* @param initial_state the initial state index if isSimOneInitState == TRUE
*
* NOTE: Below we have pointers to the pointers, these are going to be the return values
*
* @param ppYesBitSet this bitset is going to be filled with the states which
*			satisfy the formula
* @param ppNoBitSet this bitset is going to be filled with the states which
*			do not satisfy the formula
* @param ppCiLeftBorders will contain the left conf. int. borders
* @param ppCiRightBorders will contain the right conf. int. borders
* @param pResultSize the pointer to the size of arrays *ppCiLeftBorders and *ppCiRightBorders
* @param pMaxNumUsedObserv the return-value variable that will store the maximum number
* (over all initial states) of used observations
*/
void modelCheckUnboundedUntilCTMC(sparse * pStateSpace,
                                                const double * pCTMCRowSums,
						const double confidence, const bitset *pPhiBitSet,
						const bitset * pPsiBitSet, bitset ** ppYesBitsetResult,
						bitset ** ppNoBitsetResult, double ** ppProbCILeftBorder,
						double ** ppProbCIRightBorder, int * pResultSize,
						const int comparator, const double prob_bound,
                                                const int initial_state,
                                                const BOOL
                                                isSimOneInitState_local,
						unsigned int * pMaxNumUsedObserv ){
	int *pTransientStateIdx=NULL;
	/* This vector will store the exit rates for the elements of pTransientBitSet*/
	double * pExitRatesOfTransStates = NULL;
        bitset * pAUBitSet = NULL, *pBadBitSet = NULL, *pTransientBitSet = NULL;

        /* 0: Sort out the 1.0 and 0.0 prob reachable states */
	bitset * pEUBitSet = get_exist_until( pStateSpace, pPhiBitSet, pPsiBitSet );
        if ( NULL == pEUBitSet
                        || (pAUBitSet = get_always_until(pStateSpace,pPhiBitSet,
                                        pPsiBitSet, pEUBitSet)) == NULL

	/* NOTE: Using pEUBitSet there is NO NEED to compute BSCCs! */
	/* Because "pEUBitSet = get_exist_until( pStateSpace, pPhiBitSet, pPsiBitSet )" */
	/* returns the states from which it is possible to go into Psi states via Phi states */

	/* 1: Compute the states that has to be absorbing */
	/* 1.1: At this point pBadBitSet contains states from which */
	/* psi states are not reachable via phi states */
                        || (pBadBitSet = not(pEUBitSet)) == NULL
	/* 1.2: At this point pTransientBitSet contains all absorbing states, i.e. */
	/* Good and Bad absorbing states (including Phi BSCC states) */
                        || (pTransientBitSet = or(pAUBitSet, pBadBitSet))
                                        == NULL
	/* 1.3: At this point pTransientBitSet contains only allowed transient states */
	/* from which there is a way to go to Psi states via Phi states with prob < 1.0 */
                        || err_state_iserror(not_result(pTransientBitSet))

	/* print_mtx_sparse(pStateSpace); */

	/* 3: Make states absorbing */
                        || (pTransientStateIdx = count_set(pTransientBitSet))
                                        == NULL

	/* 4: Obtain the embedded DTMC (WITHOUT THE SELF LOOPS!!!) */
                        || (pExitRatesOfTransStates =
                                        make_embedded_dtmc_out_of_rate_mtx_vs(
                                        pStateSpace, pCTMCRowSums,
                                        pTransientStateIdx)) == NULL )
        {
                exit(err_macro_19(err_CALLBY, "modelCheckUnboundedUntilCTMC(%p"
                        "[%dx%d],%p,%g,%p[%d],%p[%d],%p,%p,%p,%p,%p,%d,%g,%d,"
                        "%d,%p)", (void *) pStateSpace, mtx_rows(pStateSpace),
                        mtx_cols(pStateSpace), (const void *) pCTMCRowSums,
                        confidence, (const void *) pPhiBitSet,
                        bitset_size(pPhiBitSet), (const void *) pPsiBitSet,
                        bitset_size(pPsiBitSet), (void *) ppYesBitsetResult,
                        (void *) ppNoBitsetResult, (void *) ppProbCILeftBorder,
                        (void *) ppProbCIRightBorder, (void *) pResultSize,
                        comparator, prob_bound, initial_state,
                        isSimOneInitState_local, (void *) pMaxNumUsedObserv,
                        ((void) (NULL == pEUBitSet || (free_bitset(pEUBitSet),
                            NULL == pAUBitSet || (free_bitset(pAUBitSet),
                                NULL == pBadBitSet || (free_bitset(pBadBitSet),
                                    NULL == pTransientBitSet
                                            || (free_bitset(pTransientBitSet),
                                        free(pTransientStateIdx), FALSE))))),
                        EXIT_FAILURE)));
        }

	/* print_mtx_sparse(pStateSpace); */

	/* 5: Construct the pGoodStates and pTransientStates sets */
	/* NOTE: They are already constructed, namely: pAUBitSet and pTransientBitSet */

	/* 6: Do simulations */
	modelCheckStatesCommon( pStateSpace, confidence, ppYesBitsetResult, ppNoBitsetResult, ppProbCILeftBorder,
				ppProbCIRightBorder, pResultSize, comparator, prob_bound, initial_state,
                                isSimOneInitState_local, pTransientBitSet,
                                pMaxNumUsedObserv, modelCheckOneStateUUCommon,
				pTransientBitSet, pAUBitSet );

	/* 7.1: Update ppYesBitsetResult and ppNoBitsetResult with pAUBitSet and pBadBitSet states. */
        considerAlwaysAndNeverStates(isSimOneInitState_local,
                                        * ppYesBitsetResult, * ppNoBitsetResult,
					pAUBitSet, pBadBitSet, initial_state, comparator, prob_bound );
	/* 7.2: Update ppProbCILeftBorder and ppProbCIRightBorder with probability 1.0 for pAUBitSet states. */
        updateConfIntAlwaysStates(isSimOneInitState_local, ppProbCILeftBorder,
                                         ppProbCIRightBorder, pResultSize,
					 pAUBitSet, initial_state );

	/* 8: Restore the state space */
	/* NOTE: not required for diagonal elements */
        if ( err_state_iserror(restore_rate_mtx_out_of_embedded_dtmc_vs(
                                pStateSpace, pTransientStateIdx,
                                pExitRatesOfTransStates)) )
        {
                exit(err_macro_19(err_CALLBY, "modelCheckUnboundedUntilCTMC(%p"
                        "[%dx%d],%p,%g,%p[%d],%p[%d],%p,%p,%p,%p,%p,%d,%g,%d,"
                        "%d,%p)", (void *) pStateSpace, mtx_rows(pStateSpace),
                        mtx_cols(pStateSpace), (const void *) pCTMCRowSums,
                        confidence, (const void *) pPhiBitSet,
                        bitset_size(pPhiBitSet), (const void *) pPsiBitSet,
                        bitset_size(pPsiBitSet), (void *) ppYesBitsetResult,
                        (void *) ppNoBitsetResult, (void *) ppProbCILeftBorder,
                        (void *) ppProbCIRightBorder, (void *) pResultSize,
                        comparator, prob_bound, initial_state,
                        isSimOneInitState_local, (void *) pMaxNumUsedObserv,
                        EXIT_FAILURE));
        }

	/* print_mtx_sparse(pStateSpace); */

	/* 9: Free the resources */
	free( pExitRatesOfTransStates );
	free( pTransientStateIdx );
	free_bitset( pEUBitSet );
	free_bitset( pAUBitSet );
	free_bitset( pBadBitSet );
	free_bitset( pTransientBitSet );
}

/*********************************************************************/
/********THE FUNCTIONS FOR SIMULATING INTERVAL UNTIL CTMC*************/
/*********************************************************************/

/**
 * This function computes when we will leave the state.
 * In addition to that for an absorbing current state the time of leaving the state is
 * set to be right_time_bound + 1.0 in order to avoid further simulations
 * @param current_exit_rate the exit rate of the current state
 * @param right_time_bound the right time bound of the until
 * @return the next exit time
 */
static inline double computeExitTime( const double current_exit_rate, const double right_time_bound ){
	/* If we are not in an absorbing state */
	if( current_exit_rate != 0.0 ){
		/* Compute when we leave the current state */
		return generateRandNumberExp( current_exit_rate );
	} else {
		/* If we are absorbed in a Phi state which might be also a Psi state */
		/* then we just go beyond the time interval to stop further iterations */
		return right_time_bound + 1.0;
	}
}

/**
* This function is designed to simulate the new sample observations for interval
* until CTMC: "Phi U[tl, tr] Psi".
* WARNING: We simulate the states ASSUMING there are no self loops!
* @param pStateSpace the sparse matrix of the embedded DTMC
* @param pSampleVecIntUntil the sample vector obtained on the previous iteration
* @param old_sample_size the old sample size
* @param pNotPhiAndNotPsiBitSet the bitsets containing all
*	not Phi and not Psi states
* @param pPhiBitSet the bitsets containing all the Phi
* @param pPsiBitSet the bitsets containing all the Psi
* @param left_time_bound the left time bound "tl"
* @param right_time_bound the right time bound "tr"
* @param pExitRatesOfAllowedStates the array of exit rates for the "not pNotPhiAndNotPsiBitSet"
*	states, not that its size corresponds to the dimensions of pStateSpace
*/
static void simulateIntUntilSampleCTMC( const sparse * pStateSpace, PTSampleVecIntUntil pSampleVecIntUntil,
										const int old_sample_size, const bitset * pNotPhiAndNotPsiBitSet,
										const bitset * pPhiBitSet, const bitset * pPsiBitSet,
										const double left_time_bound, const double right_time_bound,
										const double * pExitRatesOfAllowedStates ){
	int i, current_obs_state;
	double current_obs_enter_time, current_obs_exit_time;
	/* Will store the number of newly simulated observations */
	unsigned int newlyVisitedStates = 0;

	/* Cast to the parent structure which contains initial_state and curr_sample_size */
	PTSampleVec pSampleVecIntUntilBase = (PTSampleVec) pSampleVecIntUntil;

	/* For all the newly added observation in the sample */
	for( i = 0; i < ( pSampleVecIntUntilBase->curr_sample_size - old_sample_size ) ; i++ ){
		/* Take the initial state index */
		current_obs_state = pSampleVecIntUntilBase->initial_state;
		/* We start in the initial state at time 0.0 */
		current_obs_enter_time = 0.0;
		/* Compute when we leave the initial state */
		current_obs_exit_time = computeExitTime( pExitRatesOfAllowedStates[current_obs_state], right_time_bound );

		/* printf("\n----------------------------------------------------\n"); */
		/* printf("I: (%s) state %d exit at %e \n", (isInitialPartFailed? "-": "+"),
				current_obs_state+1, current_obs_exit_time ); */

		/* Iterate while we go through Phi states until left_time_bound */
		while( ( current_obs_exit_time <= left_time_bound ) && get_bit_val( pPhiBitSet, current_obs_state ) ){
			/* Go through the states until */
                        current_obs_state = computeNextState(pStateSpace,
                                                current_obs_state);
			current_obs_enter_time = current_obs_exit_time;
			current_obs_exit_time += computeExitTime( pExitRatesOfAllowedStates[current_obs_state],
									right_time_bound );
			/* The next state was simulated */
			newlyVisitedStates++;
		}

		/* printf("M: (%s) state %d exit at %e \n", (isInitialPartFailed? "-": "+"),
				current_obs_state+1, current_obs_exit_time ); */

		if( ( current_obs_enter_time == left_time_bound ) || get_bit_val( pPhiBitSet, current_obs_state ) ){
			/* Iterate while we go through Phi states */
			while( ! get_bit_val( pNotPhiAndNotPsiBitSet, current_obs_state ) &&
				! get_bit_val( pPsiBitSet, current_obs_state) &&
				( current_obs_exit_time <= right_time_bound ) ){
				/* While we are not in a bad state not in the target state and before right_time_bound */
                                current_obs_state=computeNextState(pStateSpace,
									current_obs_state );
				current_obs_exit_time += computeExitTime( pExitRatesOfAllowedStates[current_obs_state],
									right_time_bound );
				/* The next state was simulated */
				newlyVisitedStates++;
			}

			/* Update the sum_good field if we ended up in a Psi state in the right time */
			if( get_bit_val( pPsiBitSet, current_obs_state) ){
				pSampleVecIntUntil->sum_good += 1;
			}

			/* printf("F: (%s) state %d exit at %e \n",
				( get_bit_val( pPsiBitSet, current_obs_state) ? "+": "-"),
				current_obs_state+1, current_obs_exit_time ); */
		}
	}

	/* printSampleVectorIntUntil( pSampleVecIntUntil ); */

	/* Update the number of states visited while preparing this sample */
	pSampleVecIntUntilBase->num_visited_states += newlyVisitedStates;
}

/**
* This function is designed to extended and simulate the sample for interval
* until CTMC: "Phi U[tl, tr] Psi". This function is an implementation of
* the "extendSample" function from the PhD Thesis of Ivan S. Zapreev. It adds
* new observations to the sample and simulates.
* WARNING: We simulate the states ASSUMING there are no self loops!
* @param pStateSpace the sparse matrix of the embedded DTMC
* @param pSampleVecIntUntil the sample vector obtained on the previous iteration
* @param new_sample_size the new sample size
* @param pNotPhiAndNotPsiBitSet the bitsets containing all
*	not Phi and not Psi states
* @param pPhiBitSet the bitsets containing all the Phi
* @param pPsiBitSet the bitsets containing all the Psi
* @param left_time_bound the left time bound "tl"
* @param right_time_bound the right time bound "tr"
* @param pExitRatesOfAllowedStates the array of exit rates for the "not pNotPhiAndNotPsiBitSet"
*	states, not that its size corresponds to the dimensions of pStateSpace
*/
static void simulateSampleVectorIntUntilCTMC( const sparse * pStateSpace, PTSampleVecIntUntil pSampleVecIntUntil,
					const int new_sample_size, const bitset * pNotPhiAndNotPsiBitSet,
					const bitset * pPhiBitSet, const bitset * pPsiBitSet, const double left_time_bound,
					const double right_time_bound, const double * pExitRatesOfAllowedStates ){
	IF_SAFETY( pSampleVecIntUntil != NULL )
		/* Cast to the parent structure which contains initial_state and curr_sample_size */
		PTSampleVec pSampleVecIntUntilBase = (PTSampleVec) pSampleVecIntUntil;

		IF_SAFETY( ( pNotPhiAndNotPsiBitSet != NULL ) && ( pPhiBitSet != NULL ) && ( pPsiBitSet != NULL ) )
			IF_SAFETY( pStateSpace != NULL )
				IF_SAFETY( pExitRatesOfAllowedStates != NULL )
					IF_SAFETY( pSampleVecIntUntilBase->curr_sample_size < new_sample_size )
						const int old_sample_size = pSampleVecIntUntilBase->curr_sample_size;

						/* Extend the sample size */
						extendSampleVectorIntUntil( pSampleVecIntUntil, new_sample_size );

						/*Do simulations for the newly added observations */
						simulateIntUntilSampleCTMC( pStateSpace, pSampleVecIntUntil, old_sample_size,
									pNotPhiAndNotPsiBitSet, pPhiBitSet, pPsiBitSet,
									left_time_bound, right_time_bound, pExitRatesOfAllowedStates );
					ELSE_SAFETY
						printf("ERROR: The sample size is not going to increase, there is nothing to simulate.\n");
                                                exit(EXIT_FAILURE);
					ENDIF_SAFETY
				ELSE_SAFETY
					printf("ERROR: The exit-rates array is NULL.\n");
                                        exit(EXIT_FAILURE);
				ENDIF_SAFETY
			ELSE_SAFETY
				printf("ERROR: The sparse matrix is NULL.\n");
                                exit(EXIT_FAILURE);
			ENDIF_SAFETY
		ELSE_SAFETY
			printf("ERROR: The Phi, Psi or not Phi and not Psi states are undefined.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: The sample that has to be simulated further is NULL.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* Model checks the interval-until operator "Phi U[tl, tr] Psi" for one initial state only.
* Note that this is a univesral procedure that assumes that we work with the embedded DTMC.
* WARNING: We simulate the states ASSUMING there are no self loops!
* @param pStateSpace the sparse matrix of the embedded DTMC with good and bad
*			states made absorbing
* @param initial_state the initial state for simulations
* @param indiff_width the width of the indifference region
* @param confidence the desired confidence of the answer
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param arr_index indicates the index in arrays pCiLeftBorders and pCiRightBorders
*			that corresponds to the initial_state
*
* NOTE: The parameters below will contain the return values
*
* @param pCiLeftBorders for storing the left conf. int. border
* @param pCiRightBorders for storing the right conf. int. border
* @param pYesBitSet contains curr_state if it satisfies the formula
* @param pNoBitSet contains curr_state if it does not satisfy the formula
* @param pNumUsedObserv the pointer to the return-value variable that will store the
*		number of observations used for simulating this one state.
* @param pArguments contains the following (ordered!):
*		1. pNotPhiAndNotPsiBitSet the bitsets containing all
*			not Phi and not Psi states
*		2. pPhiBitSet the bitsets containing all the Phi
*		3. pPsiBitSet the bitsets containing all the Psi
*		4. left_time_bound the left time bound "tl"
*		5. right_time_bound the right time bound "tr"
*		6. pExitRatesOfAllowedStates the array of exit rates for the
*			"not pNotPhiAndNotPsiBitSet" states, not that its size corresponds
*			to the dimensions of pStateSpace
*/
static void modelCheckOneStateIUCTMC( const sparse* pStateSpace, const int initial_state, const double indiff_width,
					const double confidence, const int comparator, const double prob_bound,
					const int arr_index, double *pCiLeftBorders, double *pCiRightBorders,
					bitset * pYesBitSet, bitset * pNoBitSet, unsigned int *pNumUsedObserv,
					va_list pArguments ){
	/* Contain the pNotPhiAndNotPsiBitSet, pPhiBitSet and pPsiBitSet states respectively */
	const bitset * pNotPhiAndNotPsiBitSet = va_arg( pArguments, bitset * );
	const bitset * pPhiBitSet = va_arg( pArguments, bitset * );
	const bitset * pPsiBitSet = va_arg( pArguments, bitset * );
	/* Contain the left and the right time bounds respectively */
	const double left_time_bound = va_arg( pArguments, double );
	const double right_time_bound = va_arg( pArguments, double );
	/* Contains the exit rates for the "not pNotPhiAndNotPsiBitSet" states */
	const double * pExitRatesOfAllowedStates = va_arg( pArguments, double * );

	/* The model checking result*/
	TV_LOGIC mc_result = TVL_NN;
	/* The conf. int. left and right borders */
	double leftBorder = 0.0, rightBorder = 1.0;

	/* Compute the general confidence zeta value */
	/* We compute it for the 1 - beta/2 with confidence = 1 - beta */
	const double gen_conf_zeta = gsl_cdf_ugaussian_Pinv( 0.5 * ( 1 + confidence ) );
	/* NOTE: The confidence we have is used to derive z for the left and the */
	/* right border BUT the TRUE/FALSE answer is based only on one border of */
	/* the conf. int., i.e. we could use gsl_cdf_ugaussian_Pinv( confidence ) */
	/* for model checking, but then again we have a different confidence for the */
	/* produced conf. int. this is bad because then the user does not know that. */

	/* Initialize the sample size parameters */
	int sample_size = getSimMinSampleSize();
	const int sample_size_step = getSimSampleSizeStep();
	const int max_sample_size = getSimMaxSampleSize();

	/* Allocate an empty sample vector */
	PTSampleVecIntUntil pSampleVecIntUntil = allocateSampleVectorIntUntil( 0, initial_state );
	/* Cast to the parent structure which contains initial_state and curr_sample_size */
	PTSampleVec pSampleVecIntUntilBase = (PTSampleVec) pSampleVecIntUntil;

	/* Extend the initial vector with sample_size observations*/
	simulateSampleVectorIntUntilCTMC( pStateSpace, pSampleVecIntUntil, sample_size, pNotPhiAndNotPsiBitSet,
					pPhiBitSet, pPsiBitSet, left_time_bound, right_time_bound, pExitRatesOfAllowedStates );

	/* printSampleVectorIntUntil( pSampleVecIntUntil ); */
	/* Perform the check of the initial conf. int. */
	computeBordersConfInt( gen_conf_zeta, &leftBorder, &rightBorder, pSampleVecIntUntilBase->curr_sample_size,
				pSampleVecIntUntil->sum_good, pSampleVecIntUntil->sum_good, AGRESTI_COULL_CONF_INT);

	/* Check the conf. int. against the probability constraint */
	mc_result = checkBoundVSConfInt( comparator, prob_bound, leftBorder, rightBorder, indiff_width );

	/* The main simulation cycle */
	while( ( mc_result == TVL_NN ) && ( sample_size < max_sample_size ) ){
		/* Increase the sample size */
		increment( &sample_size, sample_size_step, max_sample_size );

		/* Extend the initial vector with sample_size observations*/
		simulateSampleVectorIntUntilCTMC( pStateSpace, pSampleVecIntUntil, sample_size, pNotPhiAndNotPsiBitSet,
						pPhiBitSet, pPsiBitSet, left_time_bound, right_time_bound, pExitRatesOfAllowedStates );

		/* printSampleVectorIntUntil( pSampleVecIntUntil ); */
		/* Perform the check of the initial conf. int. */
		computeBordersConfInt( gen_conf_zeta, &leftBorder, &rightBorder, pSampleVecIntUntilBase->curr_sample_size,
					pSampleVecIntUntil->sum_good, pSampleVecIntUntil->sum_good, AGRESTI_COULL_CONF_INT);

		/* Check the conf. int. against the probability constraint */
		mc_result = checkBoundVSConfInt( comparator, prob_bound, leftBorder, rightBorder, indiff_width );
	}
	/* Update the pCiLeftBorders and pCiRightBorders */
	pCiLeftBorders[arr_index] = leftBorder;
	pCiRightBorders[arr_index] = rightBorder;

	/* Update the pYesBitSet or pNoBitSet if we have a definite answer */
	markYesNoSetEntree(mc_result, initial_state, pYesBitSet, pNoBitSet);

	*pNumUsedObserv	= pSampleVecIntUntilBase->num_visited_states;

	/* Free the sample vectors */
	freeSampleVectorIntUntil( pSampleVecIntUntil );
}

/**
* This function is an implementation of the intervalUntil function from the PhD Thesis
* of Ivan S. Zapreev. The function does the model checking of time-interval until operator
* on the provided CTMC and uses the runtime-simulation settings from the "simulation.h"
* @param pStateSpace the sparse matrix of the embedded DTMC with good and bad
*			states made absorbing
* @param confidence this is is confidence we should use for simulations.
*			Note that since the formulas can be nested one can not use the
*			overall confidence set in simulation.h
* @param pPhiBitSet the bitsets containing all the good absorbing states
* @param pPsiBitSet the bitsets containing all the transient states
* @param left_time_bound the left time bound
* @param right_time_bound the right time bound
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param isSimOneInitState TRUE if we simulate one initial state only
* @param initial_state the initial state index if isSimOneInitState == TRUE
*
* NOTE: Below we have pointers to the pointers, these are going to be the return values
*
* @param ppYesBitSet this bitset is going to be filled with the states which
*			satisfy the formula
* @param ppNoBitSet this bitset is going to be filled with the states which
*			do not satisfy the formula
* @param ppCiLeftBorders will contain the left conf. int. borders
* @param ppCiRightBorders will contain the right conf. int. borders
* @param pResultSize the pointer to the size of arrays *ppCiLeftBorders and *ppCiRightBorders
* @param pMaxNumUsedObserv the return-value variable that will store the maximum number
* (over all initial states) of used observations
*/
void modelCheckTimeIntervalUntilCTMC(sparse * pStateSpace,
                                     const double * pCTMCRowSums,
                                     const double confidence, const bitset *pPhiBitSet,
                                     const bitset * pPsiBitSet, const double left_time_bound,
                                     const double right_time_bound, bitset ** ppYesBitsetResult,
                                     bitset ** ppNoBitsetResult, double ** ppProbCILeftBorder,
                                     double ** ppProbCIRightBorder, int * pResultSize,
                                     const int comparator, const double prob_bound,
                                     const int initial_state, const BOOL
                                     isSimOneInitState_local, unsigned int * pMaxNumUsedObserv ){
        bitset *pNotPhiAndNotPsiBitSet = NULL, * pTmpBitSet = NULL;
        /* This vector will store the exit rates for the elements of pNotPhiAndNotPsiBitSet*/
        double * pExitRatesOfAllowedStates = NULL;
        /* This array will hold indices of the possible initial states */
        int *pNAIStatesIdx = NULL;

        /* To compute pNotPhiAndNotPsiBitset directly, one needed three
           negations, and a fourth to compute its complement. Changing the order
           and calculating the complement first saves a number of bitset
           operations. David N. Jansen. */
        /* Compute possible initial-state indexes */
        pTmpBitSet = or(pPsiBitSet, pPhiBitSet);
        /* 1: Compute the set of "Always illegal" states, i.e. not Psi and not Phi states */
        pNotPhiAndNotPsiBitSet = not(pTmpBitSet);
        pNAIStatesIdx = count_set( pTmpBitSet );

        /* 2: Obtain the embedded DTMC and the vector of exit rates */
        /* (the diagonal elements stay intact) */
        pExitRatesOfAllowedStates = make_embedded_dtmc_out_of_rate_mtx_vs( pStateSpace, pCTMCRowSums, pNAIStatesIdx );
        if ( NULL == pExitRatesOfAllowedStates ) {
                exit(err_macro_21(err_CALLBY, "modelCheckTimeIntervalUntilCTMC("
                        "%p[%dx%d],%p,%g,%p[%d],%p[%d],%g,%g,%p,%p,%p,%p,%p,%d,"
                        "%g,%d,%d,%p)", (void *) pStateSpace,
                        mtx_rows(pStateSpace), mtx_cols(pStateSpace),
                        (const void *) pCTMCRowSums, confidence,
                        (const void *) pPhiBitSet, bitset_size(pPhiBitSet),
                        (const void *) pPsiBitSet, bitset_size(pPsiBitSet),
                        left_time_bound, right_time_bound,
                        (void *) ppYesBitsetResult, (void *) ppNoBitsetResult,
                        (void*) ppProbCILeftBorder, (void*) ppProbCIRightBorder,
                        (void *) pResultSize, comparator, prob_bound,
                        initial_state, isSimOneInitState_local,
                        (void *) pMaxNumUsedObserv, EXIT_FAILURE));
        }

        /* print_mtx_sparse(pStateSpace); */
        /* printf("The valid-state exit rates are listed here: "); */
        /* print_vec_double(mtx_rows(pStateSpace), pExitRatesOfAllowedStates);*/
        /* printf("\n"); */

        /* 3: Do model checking via simulations */
        modelCheckStatesCommon( pStateSpace, confidence, ppYesBitsetResult, ppNoBitsetResult,
                                ppProbCILeftBorder, ppProbCIRightBorder, pResultSize, comparator,
                                prob_bound, initial_state, isSimOneInitState_local, pTmpBitSet,
                                pMaxNumUsedObserv, modelCheckOneStateIUCTMC, pNotPhiAndNotPsiBitSet,
                                pPhiBitSet, pPsiBitSet, left_time_bound, right_time_bound,
                                pExitRatesOfAllowedStates );

        /* 4: Restore the state space */
        /* NOTE: not required for diagonal elements */
        if ( err_state_iserror(restore_rate_mtx_out_of_embedded_dtmc_vs(
                                pStateSpace, pNAIStatesIdx,
                                pExitRatesOfAllowedStates)) )
        {
                exit(err_macro_21(err_CALLBY, "modelCheckTimeIntervalUntilCTMC("
                        "%p[%dx%d],%p,%g,%p[%d],%p[%d],%g,%g,%p,%p,%p,%p,%p,%d,"
                        "%g,%d,%d,%p)", (void *) pStateSpace,
                        mtx_rows(pStateSpace), mtx_cols(pStateSpace),
                        (const void *) pCTMCRowSums, confidence,
                        (const void *) pPhiBitSet, bitset_size(pPhiBitSet),
                        (const void *) pPsiBitSet, bitset_size(pPsiBitSet),
                        left_time_bound, right_time_bound,
                        (void *) ppYesBitsetResult, (void *) ppNoBitsetResult,
                        (void*) ppProbCILeftBorder, (void*) ppProbCIRightBorder,
                        (void *) pResultSize, comparator, prob_bound,
                        initial_state, isSimOneInitState_local,
                        (void *) pMaxNumUsedObserv, EXIT_FAILURE));
        }

        /* 5: Update the ppYesBitsetResult and ppNoBitsetResult with respect to pNotPhiAndNotPsiBitSet */
        considerAlwaysAndNeverStates( isSimOneInitState_local, * ppYesBitsetResult, * ppNoBitsetResult,
                                      NULL, pNotPhiAndNotPsiBitSet, initial_state, comparator, prob_bound );

        /* 6: Free the allocated memory */
        free_bitset( pTmpBitSet );
        free_bitset( pNotPhiAndNotPsiBitSet );
        free( pExitRatesOfAllowedStates );
        free( pNAIStatesIdx );
}

/*******************************************************************************/
/********THE FUNCTIONS FOR SIMULATING STEADY-STATE OPERATOR ON A CTMC***********/
/*******************************************************************************/

/**
* This function is designed to simulate the new regeneration cycles for steady
* state CTMC: "S<>p(Psi)".
*
* WARNING: We simulate the states ASSUMING there are no self loops!
*	The absorbing states have to be found based on the absence of transitions to other states
*
* @param pStateSpace the sparse matrix of the embedded DTMC
* @param pSampleVecSteady the sample vector obtained on the previous iteration
* @param old_sample_size the old number of regeneration cycles (sample size)
* @param pPsiBitSet the bitsets containing all Psi states
* @param pExitRatesOfAllowedStates the array of exit rates for the Psi states,
*	note that its size corresponds to the dimensions of pStateSpace
*/
static void simulateSSSampleCTMC( const sparse * pStateSpace, PTSampleVecSteady pSampleVecSteady,
					const int old_samle_size, const bitset * pPsiBitSet,
					const double * pExitRatesOfAllowedStates ){
	int i, current_obs_state;
	/* Calculate the expected time spent in the i'th regeneration cycle in here */
	double accumulated_time;
	/* Calculate the accumulated number of "good" states divided by exit rates */
	/* in the i'th regeneration cycle in here */
	double accumulated_good;
	/* Will store the number of newly visited states */
	unsigned int newlyVisitedStates = 0;

	/* Cast to the parent structure which contains initial_state and curr_sample_size */
	PTSampleVec pSampleVecSteadyBase = (PTSampleVec) pSampleVecSteady;

	/* For every regeneration cycle we start in the same initial state */
	current_obs_state = pSampleVecSteadyBase->initial_state;

	/* Simulate the newly added regeneration cycles */
	for( i = old_samle_size; i < pSampleVecSteadyBase->curr_sample_size; i++ ){
		/* Set the initial counter values to zero */
		accumulated_time = 0.0;
		accumulated_good = 0.0;

		/* Simulate one regeneration cycle, i.e. until the initial state is reached again */
		do{
			accumulated_time += 1.0 / pExitRatesOfAllowedStates[current_obs_state];
			/* Increment only if the current state is a Psi state */
			if( get_bit_val( pPsiBitSet, current_obs_state ) ){
				accumulated_good += 1.0 / pExitRatesOfAllowedStates[current_obs_state];
			}

			/* Go to the next state */
                        current_obs_state = computeNextState(pStateSpace,
                                        current_obs_state);

			/* The next state was simulated */
			newlyVisitedStates++;
		} while( current_obs_state != pSampleVecSteadyBase->initial_state );

		/* printf("REGENERATION SYCLE IS DONE\n"); */

		/* Store the results for this regeneration cycle in the sample structure */
		pSampleVecSteady->pTimeObservationsVec[i] = accumulated_time;
		pSampleVecSteady->pGoodObservationsVec[i] = accumulated_good;

	}

	/* printSampleVectorSteady( pSampleVecSteady ); */

	/* Update the number of simulations used up to this point */
	pSampleVecSteadyBase->num_visited_states += newlyVisitedStates;
}

/**
* This function is designed to extended and simulate the sample for steady state
* CTMC: "S<>p(Psi)". This function is an implementation of the regenerative method
* mentioned in PhD Thesis of Ivan S. Zapreev. It adds new regeneration cycles to
* the sample and simulates.
* WARNING: We simulate the states ASSUMING there are no self loops!
* @param pStateSpace the sparse matrix of the embedded DTMC
* @param pSampleVecSteady the sample vector obtained on the previous iteration
* @param new_sample_size the new sample size
* @param pPsiBitSet the bitsets containing all Psi states
* @param pExitRatesOfAllowedStates the array of exit rates for the Psi states,
*	note that its size corresponds to the dimensions of pStateSpace
*/
static void simulateSampleVectorSSCTMC( const sparse * pStateSpace, PTSampleVecSteady pSampleVecSteady,
					const int new_sample_size, const bitset * pPsiBitSet,
					const double * pExitRatesOfAllowedStates ){
	IF_SAFETY( pSampleVecSteady != NULL )
		/* Cast to the parent structure which contains initial_state and curr_sample_size */
		PTSampleVec pSampleVecSteadyBase = (PTSampleVec) pSampleVecSteady;

		IF_SAFETY( pPsiBitSet != NULL )
			IF_SAFETY( pStateSpace != NULL )
				if( pSampleVecSteadyBase->curr_sample_size < new_sample_size ){
					const int old_sample_size = pSampleVecSteadyBase->curr_sample_size;

					/* Add new regeneration cycles */
					extendSampleVectorSteady( pSampleVecSteady, new_sample_size );

					/*Do simulations for the newly added observations */
					simulateSSSampleCTMC( pStateSpace, pSampleVecSteady, old_sample_size,
								pPsiBitSet, pExitRatesOfAllowedStates );
				}
			ELSE_SAFETY
				printf("ERROR: The sparse matrix is NULL.\n");
                                exit(EXIT_FAILURE);
			ENDIF_SAFETY
		ELSE_SAFETY
			printf("ERROR: The Psi states are undefined.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: The sample that has to be simulated further is NULL.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* This function is an implementation of the checkCIHybrid function from
* the PhD Thesis of Ivan S. Zapreev. The function does the recomputation of
* the confidence interval, afterwards check it against the probability
* constraint p given by the formula "S<>p(Psi)".
* @param pReachableTrivialBSCCs the bitset containing the indices of reachable trivial BSCCs.
* @param pReachableNonTrivBSCCs the bitset containing the indices of reachable non-trivial BSCCs.
* @param pBSCCReachProbability the array containing the reachability probabilities of every BSCC.
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param error_bound the error bounf for the numerically computed reachability
*			probabilities (Hybrid mode).
* @param pLeftCiSSBorder the array containing the left conf. int. border for every BSCC.
* @param pRightCiSSBorder the array containing the right conf. int. border for every BSCC.
* @param indiff_width the width of the indifference region
*
* NOTE: The parameters below will contain the return values

* @param pLeftBorder the overall left conf. int. border for the actual initial state.
* @param pRightBorder the overall right conf. int. border for the actual initial state.
* @return returns one of: TVL_TT, TVL_FF, TVL_NN.
*/
static TV_LOGIC checkConfIntSSHybrid(const bitset * pReachableTrivialBSCCs,
                const bitset * pReachableNonTrivBSCCs,
										const double * pBSCCReachProbability, const int comparator,
										const double prob_bound, const double error_bound,
                const double * pLeftCiSSBorder, const double * pRightCiSSBorder,
										const double indiff_width, double * pLeftBorder, double * pRightBorder){
	int i, curr_index;
	/* The confidence interval borders, which hold the end result finally */
	double LeftCiBorder = 0;
	double RightCiBorder = 0;

	int numberOfTrivialBSCCs = pReachableTrivialBSCCs->n;

	i = 0, curr_index = -1;
	/* Compute the combined confidence interval for the trivial BSCCs */
	/* First we should iterate through the trivial BSCCs since the reachability probabilities for */
	/* them are located in the pBSCCReachProbability in the very beginning */
	/* this part of the array is "indexed by" the number of BSCCs with Psi states whereas pLeftCiSSBorder */
	/* and pRightCiSSBorder are "indexed by" the number of reachable Psi BSCCs */
	while( ( curr_index =  get_idx_next_non_zero( pReachableTrivialBSCCs, curr_index ) ) != -1 ){
		LeftCiBorder += ( pBSCCReachProbability[curr_index] - error_bound ) * pLeftCiSSBorder[i];
		RightCiBorder += ( pBSCCReachProbability[curr_index] + error_bound ) * pRightCiSSBorder[i];
		i++;
	}

	/* WARNING: At this point we do not change the value of i because we */
	/* continue iterating through the pLeftCiSSBorder and pRightCiSSBorder arrays */
	curr_index = -1;
	/* Compute the combined confidence interval for the non-trivial BSCCs */
	/* Note that to get the indexes of pBSCCReachProbability from pReachableNonTrivBSCCs */
	/* related to the reachable non-trivial Psi BSCCs we have to add numberOfTrivialBSCCs to the index */
	while( ( curr_index =  get_idx_next_non_zero( pReachableNonTrivBSCCs, curr_index ) ) != -1 ){
		LeftCiBorder += ( pBSCCReachProbability[curr_index + numberOfTrivialBSCCs] - error_bound ) * pLeftCiSSBorder[i];
		RightCiBorder += ( pBSCCReachProbability[curr_index + numberOfTrivialBSCCs] + error_bound ) * pRightCiSSBorder[i];
		i++;
	}

	/* Check if we've gotten values < 0.0 and > 1.0 */
	if( LeftCiBorder < 0.0 ){
		LeftCiBorder = 0.0;
	}
	if( RightCiBorder > 1.0 ){
		RightCiBorder = 1.0;
	}

	/* Assign the return values for the combined confidence interval borders */
	( * pLeftBorder ) = LeftCiBorder;
	( * pRightBorder ) = RightCiBorder;

	/* Check the probability constraint */
	return checkBoundVSConfInt( comparator, prob_bound, LeftCiBorder, RightCiBorder, indiff_width );
}

/**
 * This method initializes the confidence intervals for our simulations and
 * also set up some parameters such as the sample size (initial value), the
 * sample size step and the maximum sample size.
 * @param numberOfReachableBSCCs the total number of reachable BSCCs
 * @param numberOfNonTrivBSCCStates the average number of states in non-trivial BSCCs
 *                                  it is used to set up the sample step size in case of
 *                                  automatic choice for that value
 * @param numberOfReachableTrivialBSCCs The number of reachable but trivial BSCCs
 * @param ppLeftBorder the pointer to the pLeftBorder parameter which will be initialized here
 * @param ppRightBorder the pointer to the pRightBorder parameter which will be initialized here
 * @param pSampleSize the pointer to the sample_size parameter which will be initialized here
 * @param pMaxSampleSize the pointer to the max__sample_size parameter which will be initialized here
 * @param pSampleSizeStep the pointer to the sample_size_step parameter which will be initialized here
 */
static void initializeSHybridSimulationConfIntRelated( const int numberOfReachableBSCCs,
														const int numberOfNonTrivBSCCStates,
														const int numberOfReachableTrivialBSCCs,
														double **ppLeftBorder, double **ppRightBorder,
														int *pSampleSize, int *pMaxSampleSize,
														int * pSampleSizeStep ) {
	/* Some internal counter */
	int i;
	/* NOTE: The confidence we have is used to derive z for the left and the */
	/* right border BUT the TRUE/FALSE answer is based only on one border of */
	/* the conf. int., i.e. we could use a smaller confidence for model checking, */
	/* but then again we have a different confidence for the produced conf. int. */
	/* this is bad because then the user does not know that. */

	/* The conf. int. left and right borders for the simulation part */
        double *pLeftBorder = (double *) calloc((size_t) numberOfReachableBSCCs,
                        sizeof(double));
        double *pRightBorder = (double*) calloc((size_t) numberOfReachableBSCCs,
                        sizeof(double));

	/* ...and their initialisation, First for the trivial BSCCs */
	for( i = 0; i < numberOfReachableTrivialBSCCs; i++ ){
		/* If the BSCC contains only one state, the long run probability is exactly 1 */
		pLeftBorder[i] = pRightBorder[i] = 1.0;
	}
	/* Next for the non-trivial BSCCs*/
	/* NOTE: There is no need to do pLeftBorder[i] = 0.0 because we use calloc */
	for(i = numberOfReachableTrivialBSCCs; i < numberOfReachableBSCCs; i++){
		pRightBorder[i] = 1.0;
	}

	/* Set the return argument values */
	(* ppLeftBorder) = pLeftBorder;
	(* ppRightBorder) = pRightBorder;
	(* pSampleSize) = getSimMinSampleSize();
	(* pMaxSampleSize) = getSimMaxSampleSize();
	/* Is the sample size step set manually or automatically? */
	if( isSimSampleSizeStepAuto() ){
		/* For automatic sample size step calculate the initial sample size step */
                *pSampleSizeStep = sqrt((double) numberOfNonTrivBSCCStates);
	} else {
		/* For manual sample size step get the value set */
		(* pSampleSizeStep) = getSimSampleSizeStep();
	}
}

/**
 * This is the body of the steady-state simulation procedure for the reachable \Psi BSCCs.
 * @param pStateSpace the state space
 * @param pAvgCycleLength the average number of the regeneration cycle, is updated here
 * @param pPsiBitSet the set of all \Psi states in the model
 * @param pExitRatesOfAllowedStates the array of the exit rates for the allowed states
 * @param pReachableNonTrivBSCCBitSet the bitset that stores the indices of
 *                                    non-trivial \Psi BSCCs
 *                                    that are also reachable from the actual initial state
 * @param numberOfReachableTrivialBSCCs the number of trivial \Psi BSCCs that are reachable
 *                                      from the actual initial state.
 * @param numberOfNonTrivBSCCStates the number of states in the reachable \Psi BSCCs
 * @param pSampleVecSteady the vector that contains simulation samples for each non-trivial BSCC
 * @param isDynamicSampleStep true if we want dynamic change in the sample_size_step
 * @param pSampleSize the pointer to the sample size, the pointed values is updated here
 * @param sample_size_step the default sample size step, as set by the user if
 *                         isDynamicSampleStep == true then it is the actual value
 *                         that is recomputed here, but initially taken from the user setting
 * @param max_sample_size the maximum allowed sample size
 * @param zeta_ss the zeta of the confidence for computing the conf int of each BSCC
 * @param pLeftBorder the array of left conf. int. borders for all reachable BSCCs
 * @param pRightBorder the array of right conf. int. borders for all reachable BSCCs
 */
static void simulateAllReachablePsiBSCCs( const sparse* const pStateSpace, int * pAvgCycleLength,
									const bitset * const pPsiBitSet, const double * const pExitRatesOfAllowedStates,
									const bitset * const pReachableNonTrivBSCCBitSet,
									const int numberOfReachableTrivialBSCCs, const int numberOfNonTrivBSCCStates,
									PTSampleVecSteady * pSampleVecSteady, BOOL isDynamicSampleStep,
									int *pSampleSize, int * pSampleSizeStep, const int max_sample_size,
									const double zeta_ss, double * pLeftBorder, double * pRightBorder ) {
	/* The local indice holders for the iterations */
	int curr_index, i;
	/* This will store the local value of the incremented sample size */
	int sample_size;
	/* This will store the local copy of the average regeneration cycle length */
	double avg_cycle_length = (* pAvgCycleLength);
	/* This will store the local value of the sample size step */
	int sample_size_step = (* pSampleSizeStep);

    /* Increment the sample size, i.e. the number of regeneration cycles */
    increment( pSampleSize, sample_size_step, max_sample_size);
    /* Get the local copy of the sample size value */
    sample_size = *pSampleSize;

    /* Simulate for every BSCC its sample vector and compute the new confidence borders*/
    curr_index = -1;
    i = numberOfReachableTrivialBSCCs;
    while((curr_index = get_idx_next_non_zero( pReachableNonTrivBSCCBitSet, curr_index) ) != -1){
        /* Do more simulations for this BSCC */
        simulateSampleVectorSSCTMC( pStateSpace, pSampleVecSteady[curr_index], sample_size, pPsiBitSet, pExitRatesOfAllowedStates );
        /* Re-compute the confidence interval borders */
        computeBordersConfIntSS( sample_size, zeta_ss, &pLeftBorder[i], &pRightBorder[i], pSampleVecSteady[curr_index] );
        /* Is the sample size step set manually or automatically? */
        if(isDynamicSampleStep){
            avg_cycle_length += ((PTSampleVec) pSampleVecSteady[curr_index] )->num_visited_states;
        }
        /* We move on to the conf. int . corresponding to the following BSCC */
        i++;
    }

    /* Is the sample size step set manually or automatically? */
    if(isDynamicSampleStep) {
        /* Calculate the average cycle length per BSCC and completed cycle */
        avg_cycle_length = avg_cycle_length / ((i - numberOfReachableTrivialBSCCs) * sample_size);
        /* For automatic sample size step calculate the new sample size step */
        sample_size_step = (numberOfNonTrivBSCCStates / avg_cycle_length) + sample_size;
    }

    /* Update the return values */
    (* pAvgCycleLength) = avg_cycle_length;
    (* pSampleSizeStep) = sample_size_step;
}

/**
* Model checks the steady-state operator "S<>p(Psi)" for one initial state only.
* Here we use hybrid approach, i.e. the reachability probabilities of BSCCs are pre computed
* and we only do regeneration simulations for the BSCCs themselves. Also note that the BSCC
* samples are reused from one call of this method to another.
* Note that this is a univesral procedure that assumes that we work with the embedded DTMC.
*
* WARNING: We simulate the states ASSUMING there are no self loops! The self loops of absorbing
* states have to be detected by the absence of other outgoing transitions.
*
* @param pStateSpace the sparse matrix of the embedded DTMC
* @param initial_state the initial state for simulations
* @param indiff_width the width of the indifference region
* @param confidence the desired confidence of the answer
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param arr_index indicates the index in arrays pCiLeftBorders and pCiRightBorders
*			that corresponds to the initial_state
*
* NOTE: The parameters below will contain the return values
*
* @param pCiLeftBorders for storing the left conf. int. border
* @param pCiRightBorders for storing the right conf. int. border
* @param pYesBitSet contains curr_state if it satisfies the formula
* @param pNoBitSet contains curr_state if it does not satisfy the formula
* @param pNumUsedObserv the pointer to the return-value variable that will store the
*		number of observations used for simulating this one state.
* @param pArguments contains the following (ordered!):
*		1. pPsiBitSet the bitset containing the Psi states
*		2. ppReachProbability The reachability probabilities from all
			initial states to BSCCs with Psi states.
*		3. error_bound the error bound for the numerically computed reachability
*			probabilities
*		4. pExitRatesOfAllowedStates the array of exit rates for the
*			"pPsiBitSet" states, note that its size corresponds
*			to the dimensions of pStateSpace
*		5. pSampleVecSteady the array of sample structures for the non-trivial BSCCs
*		6. numberOfNonTrivBSCCs the number of non-trivial BSCCs containing Psi states
*		7. numberOfTrivialBSCCs the number of trivial BSCCs consisting of good states
*		8. numberOfNonTrivBSCCStates the average number of states in the non-trivial \Psi BSCCs
*/
static void modelCheckOneStateSSHybridCTMC( const sparse* pStateSpace, const int initial_state, const double indiff_width,
											const double confidence, const int comparator, const double prob_bound,
											const int arr_index, double *pCiLeftBorders, double *pCiRightBorders,
											bitset * pYesBitSet, bitset * pNoBitSet, unsigned int *pNumUsedObserv,
											va_list pArguments ){
	/* Contains the pPsiBitSet */
	const bitset * pPsiBitSet = va_arg( pArguments, bitset * );
	/* The reachability probabilities from all initial states to BSCCs with Psi states */
	double ** ppReachProbability = va_arg( pArguments, double ** );
	/* Contains the error bound for numerical computation */
	const double error_bound = va_arg( pArguments, double );
	/* Contains the exit rates for the "pPsiBitSet" states */
	const double * pExitRatesOfAllowedStates = va_arg( pArguments, double * );
	/* Contains a sample structure for every non-trivial BSCC */
	PTSampleVecSteady * pSampleVecSteady = va_arg( pArguments, PTSampleVecSteady * );
	/* The number of non-trivial BSCCs containing Psi states */
	const int numberOfNonTrivBSCCs = va_arg( pArguments, int );
	/* The number of trivial BSCCs containing Psi states */
	const int numberOfTrivialBSCCs = va_arg( pArguments, int );
	/* The average number of states in non-trivial BSCCs */
	const int numberOfNonTrivBSCCStates = va_arg( pArguments, int );

	/* The model checking result*/
	TV_LOGIC mc_result = TVL_NN;
	/* The combined conf. int. left and right borders */
	double leftBorder = 0.0, rightBorder = 1.0;
	int curr_index;
	/* The bitset indicating the trivial BSCCs containing Psi states, which are also reachable */
	bitset * pReachableTrivialBSCCBitSet;
	/* The bitset indicating the non-trivial BSCCs containing Psi states, which are also reachable */
	bitset * pReachableNonTrivBSCCBitSet;
	/* The number of non-trivial BSCCs, that contain Psi states and are reachable from the */
	/* initial state */
	int numberOfReachableNonTrivBSCCs = 0;
	/* The number of trivial BSCCs, that contain Psi states and are reachable from the */
	/* initial state */
	int numberOfReachableTrivialBSCCs = 0;
	/* The number of overall BSCCs containing Psi states, that are reachable from the initial state */
	int numberOfReachableBSCCs = 0;

	/* The array holding the reachability probabilities from the initial state */
	/* to the reachable BSCCs. The first numberOfTrivialBSCCs are devoted to the */
	/* reachability probs. of trivial BSCCs and the rest numberOfNonTrivBSCCs */
	/* elements for the reachability probs. of non-trivial BSCCs*/
	double *pBSCCReachProbability = NULL;

	/* Extract the reachability probabilities */
	numberOfReachableBSCCs = extractReachabilityProbs( initial_state, numberOfTrivialBSCCs, numberOfNonTrivBSCCs,
														ppReachProbability, &pBSCCReachProbability,
														&pReachableTrivialBSCCBitSet, &pReachableNonTrivBSCCBitSet,
														&numberOfReachableTrivialBSCCs, &numberOfReachableNonTrivBSCCs );

	/* Do computations only if there is at least one reachable BSCC, that contains Psi states */
	if ( numberOfReachableBSCCs > 0 ){
		/* Compute the confidence zeta value */
		const double zeta_ss = gsl_cdf_ugaussian_Pinv( 0.5 * ( 1 + ROOT(confidence, numberOfReachableNonTrivBSCCs) ) );
		/* Check if we need a dynamic change for the sample size step */
		const BOOL isDynamicSampleStep =  isSimSampleSizeStepAuto();

		/*Initialise confidence intervals and related parameters */
		/* The conf. int. left and right borders for the simulation part */
		double *pLeftBorder, *pRightBorder;
		int sample_size, max_sample_size, sample_size_step;
		initializeSHybridSimulationConfIntRelated( numberOfReachableBSCCs, numberOfNonTrivBSCCStates,
													numberOfReachableTrivialBSCCs, &pLeftBorder,
													&pRightBorder, &sample_size, &max_sample_size,
													&sample_size_step );

		/* Is a definite answer already available without simulating? */
		mc_result = checkConfIntSSHybrid( pReachableTrivialBSCCBitSet, pReachableNonTrivBSCCBitSet,
						pBSCCReachProbability, comparator, prob_bound, error_bound,
						pLeftBorder, pRightBorder, indiff_width, &leftBorder, &rightBorder);

		if( ( mc_result == TVL_NN ) && ( numberOfReachableBSCCs - numberOfReachableTrivialBSCCs > 0 ) ){
			/* Perform the main simulation cycle */
			do{
				/* The average length of cycles simulated in this step */
				int avg_cycle_length = 0;

				/* Perform simulations for all of the reachable \Psi BSCCs */
				simulateAllReachablePsiBSCCs( pStateSpace, &avg_cycle_length, pPsiBitSet, pExitRatesOfAllowedStates,
												pReachableNonTrivBSCCBitSet, numberOfReachableTrivialBSCCs,
												numberOfNonTrivBSCCStates, pSampleVecSteady, isDynamicSampleStep,
												&sample_size, &sample_size_step, max_sample_size, zeta_ss,
												pLeftBorder, pRightBorder );

				/* Compute the combined confidence borders of the numerical and simulation part */
				mc_result = checkConfIntSSHybrid( pReachableTrivialBSCCBitSet, pReachableNonTrivBSCCBitSet,
												pBSCCReachProbability, comparator, prob_bound, error_bound,
												pLeftBorder, pRightBorder, indiff_width, &leftBorder, &rightBorder);
			}while( ( mc_result == TVL_NN ) && ( sample_size < max_sample_size ) );
		}

		/* Report the number of states visited while this simulation run */
		/* Note that, since the samples are reused, we take into account all the states that */
		/* had to be visited in order to give the answer to this particular initial state */
		curr_index = -1; *pNumUsedObserv = 0;
		while( ( curr_index = get_idx_next_non_zero( pReachableNonTrivBSCCBitSet, curr_index ) ) != -1 ){
			*pNumUsedObserv += ( (PTSampleVec) pSampleVecSteady[curr_index])->num_visited_states;
		}

		/* Free memory */
		free(pLeftBorder);
		free(pRightBorder);
	} else {
		/* There is no reachable BSCC containing Psi states, therefore no simulation is needed. */
		/* The long run probability is zero! */
		leftBorder = rightBorder = 0.0;
		mc_result = ( doesComparisonHold( rightBorder, comparator, prob_bound) ) ? TVL_TT : TVL_FF ;
	}

	/* Update the pCiLeftBorders and pCiRightBorders */
	pCiLeftBorders[arr_index] = leftBorder;
	pCiRightBorders[arr_index] = rightBorder;

	/* Update the pYesBitSet or pNoBitSet if we have a definite answer */
	markYesNoSetEntree(mc_result, initial_state, pYesBitSet, pNoBitSet);

	/* Free memory */
	free( pBSCCReachProbability );
	free_bitset( pReachableTrivialBSCCBitSet );
	free_bitset( pReachableNonTrivBSCCBitSet );
}

/**
* This function is an implementation of the steadyStateHybrid function from
* the PhD Thesis of Ivan S. Zapreev. The function does the model
* checking of the steady-state operator on the provided CTMC and uses the
* runtime-simulation settings from the "simulation.h"
* @param pStateSpace the sparse matrix of the CTMC
* @param confidence this is the confidence we should use for simulations.
*
* NOTE: that since the formulas can be nested one can not use the
* overall confidence set in simulation.h
*
* @param pPsiBitSet the bitsets containing the Psi states for S<>p(Psi)
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param isSimOneInitState TRUE if we simulate one initial state only
* @param initial_state the initial state index if isSimOneInitState == TRUE
* @param pFNumUnbUntilCTMC the function for computing the unbounded-until
*				probabilities numerically (Hybrid mode).
* @param pFAllowedBSCCSearchCTMC the function that is used for searching the BSCCs that contain Psi states.
* @param error_bound the error bounf for the numerically computed reachability
*			probabilities (Hybrid mode).
*
* NOTE: Below we have pointers to the pointers, these are going to be the return values
*
* @param ppYesBitSet this bitset is going to be filled with the states which
*			satisfy the formula
* @param ppNoBitSet this bitset is going to be filled with the states which
*			do not satisfy the formula
* @param ppCiLeftBorders will contain the left conf. int. borders
* @param ppCiRightBorders will contain the right conf. int. borders
* @param pResultSize the pointer to the size of arrays *ppCiLeftBorders and *ppCiRightBorders
* @param pMaxNumUsedObserv the return-value variable that will store the maximum number
* (over all initial states) of used observations
*/
void modelCheckSteadyStateHybridCTMC( sparse * pStateSpace, const double * pCTMCRowSums ,
                const double confidence, const bitset * pPsiBitSet,
									bitset ** ppYesBitsetResult, bitset ** ppNoBitsetResult,
									double ** ppProbCILeftBorder, double ** ppProbCIRightBorder,
									int * pResultSize, const int comparator, const double prob_bound,
                const int initial_state, const BOOL isSimOneInitState_local,
									const TPFunctNumUnbUntilCTMC pFNumUnbUntilCTMC,
									const TPFunctAllowedBSCCSearchCTMC pFAllowedBSCCSearchCTMC,
									const double error_bound, unsigned int * pMaxNumUsedObserv) {
	/*Some technical variables*/
	int i;
	/* The number of BSCCs containing Psi states */
	int numberOfBSCCs;
	/* The number of non-trivial BSCCs containing Psi states */
	int numberOfNonTrivBSCCs = 0;
	/* The average number of states in non-trivial BSCCs containing Psi states */
	int numberOfNonTrivBSCCStates = 0;
	/* The number of trivial BSCCs containing Psi states */
	int numberOfTrivialBSCCs = 0;
	/* The structure, that stores the set of states belonging to non-trivial BSCC with array index i */
	bitset ** ppNonTrivBSCCBitSets = NULL;
	/* The bitset indicating the Psi states belonging to trivial BSCCs */
	bitset * pTrivialBSCCBitSet = NULL;
        /* The array that holds the indices of states belonging to non-trivial
           BSCCs */
	int * pNonTrivBSCCStates = NULL;
	/* The structure, that stores the probability to reach BSCC i from state j */
	double ** ppReachProbability;
	/* This vector will store the exit rates for the elements of pPsiBitSet*/
	double * pExitRatesOfAllowedStates = NULL;
	/* The state space dimension */
        const int n_states = mtx_rows(pStateSpace);

	/* The structure, which stores the samples */
	PTSampleVecSteady *pSampleVecSteady = NULL;

	/* The bitset that holds the possible initial states */
	bitset *pValidInitialStatesBitSet = get_new_bitset( n_states );
	/* All states are possible initial states */
	fill_bitset_one( pValidInitialStatesBitSet );

	/* Get the trivial and non=trivial BSCCs that contain \Psi states */
	ppNonTrivBSCCBitSets = getTrivialAndNonTrivialPsiBSCC( pStateSpace, pPsiBitSet, &pNonTrivBSCCStates,
															&pTrivialBSCCBitSet, &numberOfNonTrivBSCCs,
															&numberOfTrivialBSCCs, &numberOfBSCCs,
															pFAllowedBSCCSearchCTMC );

	/* Obtain the embedded DTMC and the vector of exit rates */
	/* WARNING: The self loops on states are not modified. Thus the proper way */
	/* to identify an absorbing state is to make sure that there are */
	/* no non-diagonal elements in the corresponding row. */
	pExitRatesOfAllowedStates = make_embedded_dtmc_out_of_rate_mtx_vs( pStateSpace, pCTMCRowSums, pNonTrivBSCCStates );
        if ( NULL == pExitRatesOfAllowedStates ) {
                exit(err_macro_18(err_CALLBY, "modelCheckSteadyStateHybridCTMC("
                        "%p[%dx%d],%p,%g,%p[%d],%p,%p,%p,%p,%p,%d,%g,%d,%d,_,_,"
                        "%g,%p)", (void*) pStateSpace, mtx_rows(pStateSpace),
                        mtx_cols(pStateSpace), (const void *) pCTMCRowSums,
                        confidence, (const void *) pPsiBitSet,
                        bitset_size(pPsiBitSet), (void *) ppYesBitsetResult,
                        (void *) ppNoBitsetResult, (void *) ppProbCILeftBorder,
                        (void *) ppProbCIRightBorder, (void *) pResultSize,
                        comparator, prob_bound, initial_state,
                        isSimOneInitState_local, error_bound,
                        (void *) pMaxNumUsedObserv, EXIT_FAILURE));
        }

	/* Prepare the steady state simulation samples vector for */
	/* each BSCC and choose regeneration states for each BSCC */
	/* Also count the number of the non-trivial BSCCs */
	pSampleVecSteady = prepareSSSimulationSample( pStateSpace, numberOfNonTrivBSCCs, ppNonTrivBSCCBitSets, &numberOfNonTrivBSCCStates );

    /* Calculate the average number of states in a BSCC */
	if( numberOfNonTrivBSCCs > 0 ){
		numberOfNonTrivBSCCStates = numberOfNonTrivBSCCStates / numberOfNonTrivBSCCs;
	}

	/* Compute the reachability probabilities for the BSCCs from each and every state */
	ppReachProbability = getBSCCReachabilityProbsNumerical(n_states, numberOfBSCCs, numberOfTrivialBSCCs,
															ppNonTrivBSCCBitSets, pTrivialBSCCBitSet,
															pFNumUnbUntilCTMC );

	/* Call the common model-checking procedure for all (one) state(s) */
    modelCheckStatesCommon( pStateSpace, confidence, ppYesBitsetResult, ppNoBitsetResult, ppProbCILeftBorder,
							ppProbCIRightBorder, pResultSize, comparator, prob_bound, initial_state,
                        isSimOneInitState_local, pValidInitialStatesBitSet,
                        pMaxNumUsedObserv,
							modelCheckOneStateSSHybridCTMC, pPsiBitSet, ppReachProbability, error_bound,
							pExitRatesOfAllowedStates, pSampleVecSteady, numberOfNonTrivBSCCs,
							numberOfTrivialBSCCs, numberOfNonTrivBSCCStates );

	/* Restore the state space */
        /* NOTE: The diagonal elements do not have to be restored since they
           have not been changed */
        if ( err_state_iserror(restore_rate_mtx_out_of_embedded_dtmc_vs(
                                pStateSpace, pNonTrivBSCCStates,
                                pExitRatesOfAllowedStates)) )
        {
                exit(err_macro_18(err_CALLBY, "modelCheckSteadyStateHybridCTMC("
                        "%p[%dx%d],%p,%g,%p[%d],%p,%p,%p,%p,%p,%d,%g,%d,%d,_,_,"
                        "%g,%p)", (void*) pStateSpace, mtx_rows(pStateSpace),
                        mtx_cols(pStateSpace), (const void *) pCTMCRowSums,
                        confidence, (const void *) pPsiBitSet,
                        bitset_size(pPsiBitSet), (void *) ppYesBitsetResult,
                        (void *) ppNoBitsetResult, (void *) ppProbCILeftBorder,
                        (void *) ppProbCIRightBorder, (void *) pResultSize,
                        comparator, prob_bound, initial_state,
                        isSimOneInitState_local, error_bound,
                        (void *) pMaxNumUsedObserv, EXIT_FAILURE));
        }

	/* Deallocate the used memory */
	free( pExitRatesOfAllowedStates ); pExitRatesOfAllowedStates = NULL;
	free( pNonTrivBSCCStates ); pNonTrivBSCCStates = NULL;
	free_bitset( pTrivialBSCCBitSet ); pTrivialBSCCBitSet = NULL;
	free_bitset( pValidInitialStatesBitSet ); pValidInitialStatesBitSet = NULL;

	for( i = 0; i < numberOfNonTrivBSCCs; i++ ){
		free_bitset( ppNonTrivBSCCBitSets[i] );
	}
	free( ppNonTrivBSCCBitSets ); ppNonTrivBSCCBitSets = NULL;

	for( i = 0; i < numberOfBSCCs; i++ ){
		free( ppReachProbability[i] );
	}
	free( ppReachProbability ); ppReachProbability = NULL;

	for( i = 0; i < numberOfNonTrivBSCCs; i++){
		freeSampleVectorSteady( pSampleVecSteady[i] );
	}
	free( pSampleVecSteady ); pSampleVecSteady = NULL;
}

/**
 * This method initializes the confidence intervals for our simulations and
 * also sets up some parameters such as the sample size (initial value), the
 * sample size step and the maximum sample size.
 * @param numberOfReachableBSCCs the total number of reachable BSCCs
 * @param numberOfNonTrivBSCCStates the average number of states in non-trivial BSCCs
 *                                  it is used to set up the sample step size in case of
 *                                  automatic choice for that value
 * @param numberOfReachableTrivialBSCCs The number of reachable but trivial BSCCs
 * @param pppSampleVecUntilOneArray the pointer to the array of left conf. int. sample
 *                                  pointers for the unbounded-until simulation.
 * @param pppSampleVecUntilTwoArray the pointer to the array of left conf. int. sample
 *                                  pointers for the unbounded-until simulation.
 * @param ppLeftBorderUU the pointer to the pLeftBorderUU parameter which will be initialized here
 * @param ppRightBorderUU the pointer to the pRightBorderUU parameter which will be initialized here
 * @param ppLeftBorderSS the pointer to the pLeftBorderSS parameter which will be initialized here
 * @param ppRightBorderSS the pointer to the pRightBorderSS parameter which will be initialized here
 * @param pSampleSizeUU the pointer to the sample_size_uu parameter which will be initialized here
 * @param pSimDepthUU the pointer to the sim_depth_uu parameter which will be initialized here
 * @param pSampleSizeSS the pointer to the sample_size_ss parameter which will be initialized here
 * @param pMaxSampleSizeUU the pointer to the max_sample_size_uu parameter which will be initialized here
 * @param pMaxSimDepthUU the pointer to the max_sim_depth_uu parameter which will be initialized here
 * @param pMaxSampleSizeSS the pointer to the max_sample_size_ss parameter which will be initialized here
 * @param pSampleSizeStepUU the pointer to the sample_size_step_uu parameter which will be initialized here
 *							Here the unbounded-until simulation only uses the user defined value.
 * @param pSimDepthStepUU the pointer to the sim_depth_step_uu parameter which will be initialized here
 * @param pSampleSizeStepSS the pointer to the sample_size_step_ss parameter which will be initialized here
 *							Here the steady-state simulation can use a user defined value or a dynamic one.
 */
static void initializeSPureSimulationConfIntRelated( const int initial_state, const int numberOfReachableBSCCs,
														const int numberOfNonTrivBSCCStates,
														const int numberOfReachableTrivialBSCCs,
														PTSampleVecUntil **pppSampleVecUntilOneArray,
														PTSampleVecUntil **pppSampleVecUntilTwoArray,
														double **ppLeftBorderUU, double **ppRightBorderUU,
														double **ppLeftBorderSS, double **ppRightBorderSS,
														int *pSampleSizeUU, int * pSimDepthUU, int * pSampleSizeSS,
														int *pMaxSampleSizeUU, int * pMaxSimDepthUU, int *pMaxSampleSizeSS,
														int *pSampleSizeStepUU, int *pSimDepthStepUU, int *pSampleSizeStepSS,
														const BOOL isAlwaysReachableBSCC ) {
	/* Some internal counter */
	int i;
	/* The sample vector arrays for reachability simulations */
	PTSampleVecUntil *ppSampleVecUntilOneArray, *ppSampleVecUntilTwoArray;

	/* The conf. int. left and right borders for the simulation part */
        double * pLeftBorderUU = (double *) calloc(
                        (size_t) numberOfReachableBSCCs, sizeof(double));
        double * pRightBorderUU = (double *) calloc(
                        (size_t) numberOfReachableBSCCs, sizeof(double));

	/* Initialize the sample sizes and conf. int. and etc the same way as for the hybrid case */
	initializeSHybridSimulationConfIntRelated( numberOfReachableBSCCs, numberOfNonTrivBSCCStates,
												numberOfReachableTrivialBSCCs, ppLeftBorderSS,
												ppRightBorderSS, pSampleSizeSS, pMaxSampleSizeSS,
												pSampleSizeStepSS );

	/* Here, for UU simulations, we mostly use the same sample-related values as for SS */
	/* This is done because we do not have distinguished setting in the tool's interface */
	/* to set all of these options separately. Another thing is that the sample size step */
	/* for UU simulation can not be dynamic as for SS, so we set it to the fixed value*/
	if( isAlwaysReachableBSCC ) {
		/* In case the is a BSCC that is always reachable from the actual initial state */
		/* We will not need to do any reachability (UU) simulations so we set the size to 0 */
		(* pSampleSizeUU) = 0;
	} else {
		/* In this case we will have to do reachability simulations */
		(* pSampleSizeUU) = (* pSampleSizeSS);
	}
	(* pMaxSampleSizeUU) = (* pMaxSampleSizeSS);
	(* pSampleSizeStepUU) = getSimSampleSizeStep();
	/* Initialise the sample vector arrays for the reachability simulations */
        ppSampleVecUntilOneArray = (PTSampleVecUntil *) calloc(
                (size_t) numberOfReachableBSCCs, sizeof(PTSampleVecUntil));
        ppSampleVecUntilTwoArray = (PTSampleVecUntil *) calloc(
                (size_t) numberOfReachableBSCCs, sizeof(PTSampleVecUntil));
	for( i = 0; i < numberOfReachableBSCCs; i++ ) {
		ppSampleVecUntilOneArray[i] = allocateSampleVectorUntil( (*pSampleSizeUU), initial_state );
		ppSampleVecUntilTwoArray[i] = allocateSampleVectorUntil( (*pSampleSizeUU), initial_state );
		if( isAlwaysReachableBSCC ) {
			/* If there is one BSCC to which we always go from the given initial state */
			pLeftBorderUU[i] = 1.0;
			pRightBorderUU[i] = 1.0;
		} else {
			/* If there is no one BSCC to which we always go from the given initial state */
			pLeftBorderUU[i] = 0.0;
			pRightBorderUU[i] = 1.0;
		}
	}

	/* Set the return argument values */
	(* pppSampleVecUntilOneArray) = ppSampleVecUntilOneArray;
	(* pppSampleVecUntilTwoArray) = ppSampleVecUntilTwoArray;

	(* ppLeftBorderUU) = pLeftBorderUU;
	(* ppRightBorderUU) = pRightBorderUU;

	(* pSimDepthUU) = getSimMinSimulationDepth();
	(* pMaxSimDepthUU) = getSimMaxSimulationDepth();
	(* pSimDepthStepUU) = getSimSimulationDepthStep();
}


/**
* This function is an implementation of the checkCIPure function from
* the PhD Thesis of Ivan S. Zapreev. The function does the recomputation of
* the confidence interval, afterwards checks it against the probability
* constraint p given by the formula "S<>p(Psi)".
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param error_bound the error bounf for the numerically computed reachability
*			probabilities (Hybrid mode).
* @param pLeftCiUUBorder the array containing the left conf. int. border for the
*                       reachability problem of every reachable \Psi BSCC.
* @param pRightCiUUBorder the array containing the right conf. int. border for the
*                       reachability problem of every reachable \Psi BSCC.
* @param pLeftCiSSBorder the array containing the left conf. int. border for the
*                        steady state problem of every reachable \Psi BSCC.
* @param pRightCiSSBorder the array containing the right conf. int. border for the
*                        steady state problem of every reachable \Psi BSCC.
* @param indiff_width the width of the indifference region
* @param numberOfReachableBSCCs the number of reachable \Psi BSCCs i.e. the size of
*                               the arrays pLeftCiUUBorder, pRightCiUUBorder,
*                               pLeftCiSSBorder, and pRightCiSSBorder.
*
* NOTE: The parameters below will contain the return values
*
* @param pLeftBorder the overall left conf. int. border for the actual initial state.
* @param pRightBorder the overall right conf. int. border for the actual initial state.
* @return returns one of: TVL_TT, TVL_FF, TVL_NN.
*/
static TV_LOGIC checkConfIntSSPure( const int comparator, const double prob_bound,
									double * pLeftCiUUBorder, double * pRightCiUUBorder,
									double * pLeftCiSSBorder, double * pRightCiSSBorder,
									const double indiff_width, const int numberOfReachableBSCCs,
									double * pLeftBorder, double * pRightBorder){
	/* The internal counter */
	int i;
	/* The confidence interval borders, which hold the end result finally */
	double LeftCiBorder = 0;
	double RightCiBorder = 0;

	for( i = 0; i < numberOfReachableBSCCs; i ++ ) {
		LeftCiBorder +=  pLeftCiUUBorder[i] * pLeftCiSSBorder[i];
		RightCiBorder += pRightCiUUBorder[i] * pRightCiSSBorder[i];
	}

	/* Check if we've gotten values < 0.0 and > 1.0 */
	if( LeftCiBorder < 0.0 ){
		LeftCiBorder = 0.0;
	}
	if( RightCiBorder > 1.0 ){
		RightCiBorder = 1.0;
	}

	/* Assign the return values for the combined confidence interval borders */
	( * pLeftBorder ) = LeftCiBorder;
	( * pRightBorder ) = RightCiBorder;

	/* Check the probability constraint */
	return checkBoundVSConfInt( comparator, prob_bound, LeftCiBorder, RightCiBorder, indiff_width );
}

/**
* Model checks the steady-state operator "S<>p(Psi)" for one initial state only.
* Here we take pure simulation approach, we do not use numerically computed probs.
* Also note that the BSCC samples are reused from one call of this method to another.
* Note that this is a universal procedure that assumes that we work with the embedded DTMC.
*
* WARNING: We simulate the states ASSUMING there are no self loops! The self loops of absorbing
* states have to be detected by the absence of other outgoing transitions.
*
* @param pStateSpace the sparse matrix of the embedded DTMC
* @param initial_state the initial state for simulations
* @param indiff_width the width of the indifference region
* @param confidence the desired confidence of the answer
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param arr_index indicates the index in arrays pCiLeftBorders and pCiRightBorders
*			that corresponds to the initial_state
*
* NOTE: The parameters below will contain the return values
*
* @param pCiLeftBorders for storing the left conf. int. border
* @param pCiRightBorders for storing the right conf. int. border
* @param pYesBitSet contains curr_state if it satisfies the formula
* @param pNoBitSet contains curr_state if it does not satisfy the formula
* @param pNumUsedObserv the pointer to the return-value variable that will store the
*		number of observations used for simulating this one state.
* @param pArguments contains the following (ordered!):
*		1. pPsiBitSet the bitset containing the Psi states
*		2. ppBSCCExistsReachableStates the array that stores bitsets (at index i) with states
			from which the BSCC_{i} is at least sometimes reachable.
*		3. ppBSCCAlwaysReachableStates the array that stores bitsets (at index i) with states
			from which the BSCC_{i} is always reachable.
*		4. pExitRatesOfAllStates the array of exit rates for the
*			"pPsiBitSet" states, note that its size corresponds
*			to the dimensions of pStateSpace
*		5. pSampleVecSteady the array of sample structures for the non-trivial BSCCs
*		6. numberOfNonTrivBSCCs the number of non-trivial BSCCs containing Psi states
*		7. numberOfTrivialBSCCs the number of trivial BSCCs consisting of good states
*		8. numberOfNonTrivBSCCStates the average number of states in the non-trivial \Psi BSCCs
*		9. pTransientBitSet the set of transient states of the model
*		10. ppNonTrivBSCCBitSets the array of bitsets each of which contains states of a BSCC with \Psi states
*		11. pTrivialBSCCBitSet the bitset of states each of which is a trivial Psi BSCC
 */
static void modelCheckOneStateSSPureCTMC( const sparse* pStateSpace, const int initial_state, const double indiff_width,
										const double confidence, const int comparator, const double prob_bound,
										const int arr_index, double *pCiLeftBorders, double *pCiRightBorders,
										bitset * pYesBitSet, bitset * pNoBitSet, unsigned int *pNumUsedObserv,
										va_list pArguments ) {
	/* Contains the pPsiBitSet */
	const bitset * const pPsiBitSet = va_arg( pArguments, bitset * );
	/* The arrays indexed by the BSCC number (first for trivial and then non-trivial BSCCs) */
	/* The sub array at index i contains states satisfying  \exists and \always true U BSCC_{s} */
	/* in ppBSCCExistsReachableStates[i] and ppBSCCAlwaysReachableStates[i] respectively */
	const bitset ** const ppBSCCExistsReachableStates = va_arg( pArguments, const bitset ** const );
	const bitset ** const ppBSCCAlwaysReachableStates = va_arg( pArguments, const bitset ** const );
	/* Contains the exit rates for the "pPsiBitSet" states */
	const double * const pExitRatesOfAllStates = va_arg( pArguments, double * );
	/* Contains a sample structure for every non-trivial BSCC */
	PTSampleVecSteady * pSampleVecSteady = va_arg( pArguments, PTSampleVecSteady * );
	/* The number of non-trivial BSCCs containing Psi states */
	const int numberOfNonTrivBSCCs = va_arg( pArguments, int );
	/* The number of trivial BSCCs containing Psi states */
	const int numberOfTrivialBSCCs = va_arg( pArguments, int );
	/* The average number of states in non-trivial BSCCs */
	const int numberOfNonTrivBSCCStates = va_arg( pArguments, int );
	/* Get the states that do not belong to any BSCC */
	const bitset * const pTransientBitSet = va_arg( pArguments, const bitset * const );
	/* The array of bitsets each of which contains states belonging to some*/
	/* non-trivial BSCC with \Psi states in it */
	const bitset ** const ppNonTrivBSCCBitSets = va_arg( pArguments, const bitset ** const  );
	/* The bitsets which contains states belonging to trivial BSCCs a \Psi state each */
	const bitset * const pTrivialBSCCBitSet = va_arg( pArguments, const bitset * const );

	/* Internally used iterator and index variables */
	int i, curr_index;
	/* Indicates that we have gotted an invalid conf. int. which we can not fix within the */
	/* given  max number of samples and so we need to terminate the model checking procedure */
	BOOL isInvalid = FALSE;
	/* The model checking result*/
	TV_LOGIC mc_result = TVL_NN;
	/* The combined conf. int. left and right borders */
	double leftBorder = 0.0, rightBorder = 1.0;
	/* The number of non-trivial BSCCs, that contain Psi states and */
	/* are some times reachable from the initial state */
	int numberOfReachableNonTrivBSCCs = 0;
	/* The number of trivial BSCCs, that contain Psi states and */
	/* are some times reachable from the initial state */
	int numberOfReachableTrivialBSCCs = 0;
	/* The number of overall BSCCs containing Psi states, that are reachable from the initial state */
	int numberOfReachableBSCCs = 0;

	/* The bitset indicating the trivial BSCCs containing Psi */
	/* states, which are some times/ always reachable */
	bitset * pReachableTrivialBSCCBitSet = NULL;
	/* The bitset indicating the non-trivial BSCCs containing Psi */
	/* states, which are some times/always reachable */
	bitset * pReachableNonTrivBSCCBitSet = NULL;

	/* Extract the reachability probabilities */
	BOOL isAlwaysReachableBSCC = extractReachabilityData( initial_state, numberOfTrivialBSCCs, numberOfNonTrivBSCCs,
														ppBSCCExistsReachableStates, ppBSCCAlwaysReachableStates,
														&pReachableTrivialBSCCBitSet, &pReachableNonTrivBSCCBitSet,
														&numberOfReachableTrivialBSCCs, &numberOfReachableNonTrivBSCCs);
	numberOfReachableBSCCs = (numberOfReachableTrivialBSCCs + numberOfReachableNonTrivBSCCs);
	/* printf("EXTRACTED THE REACHABLE PSI BSCCs FOR STATE %d\n", initial_state+1);
	printf("TOTAL NUMBER OF REACHABLE PSI BSCCs: %d\n", numberOfReachableBSCCs);
	printf("NUMBER OF REACHABLE TRIVIAL PSI BSCCs: %d\n", numberOfReachableTrivialBSCCs);
	printf("THE REACHABLE TRIVIAL PSI BSCCs: ");
	print_bitset_states(pReachableTrivialBSCCBitSet); printf("\n");
	printf("NUMBER OF REACHABLE NON-TRIVIAL PSI BSCCs: %d\n", numberOfReachableNonTrivBSCCs);
	printf("THE REACHABLE NON-TRIVIAL PSI BSCCs (index+numberOfTrivBSCCs): ");
	print_bitset_states(pReachableNonTrivBSCCBitSet); printf("\n"); */

	/* Do computations only if there is at least one reachable BSCC, that contains Psi states */
	if ( numberOfReachableBSCCs > 0 ) {
		/* Check if we need a dynamic change for the sample size step */
		const BOOL isDynamicSampleStepSS =  isSimSampleSizeStepAuto();
		/* The conf. int. left and right borders for the simulation parts and other simulation params */
		double *pLeftBorderUU = NULL, *pRightBorderUU = NULL, *pLeftBorderSS = NULL, *pRightBorderSS = NULL;
		int sample_size_uu, sim_depth_uu, sample_size_ss;
		int max_sample_size_uu, max_sim_depth_uu, max_sample_size_ss;
		int sample_size_step_uu, sim_depth_step_uu, sample_size_step_ss;
		PTSampleVecUntil *pSampleVecUntilOneArray = NULL, *pSampleVecUntilTwoArray = NULL;

		/* Compute the zeta values for the UU and SS simulations */
		double zeta_ss, zeta_uu;
		if( isAlwaysReachableBSCC ) {
			/* This case is the same as for the hybrid approach */
			zeta_ss = gsl_cdf_ugaussian_Pinv( 0.5 * ( 1 + ROOT(confidence, numberOfReachableNonTrivBSCCs) ) );
			/* And we do not perform any reachability simulation here at all*/
			zeta_uu = 0.0;
		} else {
			/*TODO: We could perhaps reduce the zeta values here if we would take into account */
			/* that there are some trivial BSCCs for which we do not need to do simulations */
			zeta_ss =  gsl_cdf_ugaussian_Pinv( 0.5 * ( 1 + ROOT( confidence , 2 * numberOfReachableBSCCs ) ) );
			zeta_uu = gsl_cdf_ugaussian_Pinv( ROOT( confidence, 4 * numberOfReachableBSCCs ) );
		}

		/* Initialise UU samples, confidence intervals, and other related parameters */
		initializeSPureSimulationConfIntRelated( initial_state, numberOfReachableBSCCs, numberOfNonTrivBSCCStates,
												numberOfReachableTrivialBSCCs, &pSampleVecUntilOneArray,
												&pSampleVecUntilTwoArray, &pLeftBorderUU, &pRightBorderUU,
												&pLeftBorderSS, &pRightBorderSS, &sample_size_uu, &sim_depth_uu,
												&sample_size_ss, &max_sample_size_uu, &max_sim_depth_uu, &max_sample_size_ss,
												&sample_size_step_uu, &sim_depth_step_uu, &sample_size_step_ss,
												isAlwaysReachableBSCC );

		/* Perform the main simulation cycle */
		do{
			/* The odd iteration marker */
			BOOL isOdd = FALSE;
			/* The average length of cycles simulated in this step */
			int avg_cycle_length = 0;

			/* Do UU simulation first, for all reachable BSCCs, if somewhere we get */
			/* invalid conf. int. then it means we have to terminate the algorithm. */
			/* This is because we must have tried to fix the conf. int. but failed */
			if( ! isAlwaysReachableBSCC ) {
				if( !( isInvalid = simulateAndComputeAllReachabilityConfInt( pStateSpace, pTransientBitSet, pTrivialBSCCBitSet,
									ppNonTrivBSCCBitSets, isOdd, pReachableTrivialBSCCBitSet, pReachableNonTrivBSCCBitSet,
									pSampleVecUntilOneArray, pSampleVecUntilTwoArray, pLeftBorderUU, pRightBorderUU, &sample_size_uu,
									max_sample_size_uu, sample_size_step_uu, &sim_depth_uu, max_sim_depth_uu, sim_depth_step_uu, zeta_uu ) ) ) {
					/* If we simulated the reachability problems and did not get invalid conf. int. then */
					/* check the current values of the conf. int. may be we can already give the answer */
					mc_result = checkConfIntSSPure( comparator, prob_bound, pLeftBorderUU, pRightBorderUU,
													pLeftBorderSS, pRightBorderSS, indiff_width, numberOfReachableBSCCs,
													&leftBorder, &rightBorder);
				}
			}

			/* If we successfully did the reachability simulation then move on */
			if( ! isInvalid ) {
				if( mc_result == TVL_NN ) {
					/* If not, then do steady state simulations an update their conf. int. */
					/* Perform simulations for all of the reachable \Psi BSCCs */
					simulateAllReachablePsiBSCCs( pStateSpace, &avg_cycle_length, pPsiBitSet, pExitRatesOfAllStates,
													pReachableNonTrivBSCCBitSet, numberOfReachableTrivialBSCCs,
													numberOfNonTrivBSCCStates, pSampleVecSteady, isDynamicSampleStepSS,
													&sample_size_ss, &sample_size_step_ss, max_sample_size_ss, zeta_ss,
													pLeftBorderSS, pRightBorderSS );

					/* Check the current values of the conf. int. may be we can give the answer now */
					mc_result = checkConfIntSSPure( comparator, prob_bound, pLeftBorderUU, pRightBorderUU,
													pLeftBorderSS, pRightBorderSS, indiff_width, numberOfReachableBSCCs,
													&leftBorder, &rightBorder);
				}
			}

			isOdd = ! isOdd;
		} while ( ( !isInvalid) && ( mc_result == TVL_NN ) && ( ( sample_size_uu < max_sample_size_uu ) ||
						( sim_depth_uu < max_sim_depth_uu ) || ( sample_size_ss < max_sample_size_ss ) ) );

		/* Report the number of states visited while this simulation run */
		/* Note that, since the samples are reused, we take into account all the states that */
		/* had to be visited in order to give the answer to this particular initial state */
		curr_index = -1; *pNumUsedObserv = 0;
		/* Cound the steady state simulation samples */
		while( ( curr_index = get_idx_next_non_zero( pReachableNonTrivBSCCBitSet, curr_index ) ) != -1 ){
			*pNumUsedObserv += ( (PTSampleVec) pSampleVecSteady[curr_index])->num_visited_states;
			/* printf("Adding observations of the SS sample %d, the new total is %d\n", curr_index, (*pNumUsedObserv) ); */
		}
		/* Count the reachability simulation samples */
		for( i = 0; i < numberOfReachableBSCCs; i++ ) {
			*pNumUsedObserv +=  ( (PTSampleVec) pSampleVecUntilOneArray[i] )->num_visited_states +
								( (PTSampleVec) pSampleVecUntilTwoArray[i] )->num_visited_states;
			/* printf("Adding observations of the UU samples %d, the new total is %d\n", i, (*pNumUsedObserv) ); */
		}

		/* Free the memory */
		for( i = 0; i < numberOfReachableBSCCs; i++ ){
			freeSampleVectorUntil( pSampleVecUntilOneArray[i] );
			freeSampleVectorUntil( pSampleVecUntilTwoArray[i] );
		}
		free( pSampleVecUntilOneArray ); pSampleVecUntilOneArray = NULL;
		free( pSampleVecUntilTwoArray ); pSampleVecUntilTwoArray = NULL;
		free( pLeftBorderUU ); pLeftBorderUU = NULL;
		free( pRightBorderUU ); pRightBorderUU = NULL;
		free( pLeftBorderSS ); pLeftBorderSS = NULL;
		free( pRightBorderSS ); pRightBorderSS = NULL;
	} else {
		/* There is no reachable BSCC containing Psi states, therefore no simulation is needed. */
		/* The long run probability is zero! */
		leftBorder = rightBorder = 0.0;
		mc_result = ( doesComparisonHold( rightBorder, comparator, prob_bound ) ) ? TVL_TT : TVL_FF ;
	}

	/* Update the pCiLeftBorders and pCiRightBorders and if needed the result */
	if( isInvalid ) {
		mc_result = TVL_NN;
		pCiLeftBorders[arr_index] = 0.0;
		pCiRightBorders[arr_index] = 1.0;
	} else {
		pCiLeftBorders[arr_index] = leftBorder;
		pCiRightBorders[arr_index] = rightBorder;
	}

	/* Update the pYesBitSet or pNoBitSet if we have a definite answer */
	markYesNoSetEntree( mc_result, initial_state, pYesBitSet, pNoBitSet );

	/* Free memory */
	free_bitset( pReachableTrivialBSCCBitSet ); pReachableTrivialBSCCBitSet = NULL;
	free_bitset( pReachableNonTrivBSCCBitSet ); pReachableNonTrivBSCCBitSet = NULL;
}

/********************************************************************************************************************/
/**************SOME IMPORTANT NOTES ABOUT THE IMPLEMENTATION OF SIMULATIONS FOR THE STEADY-OPERATOR******************/
/********************************************************************************************************************/
/*1. pTrivialBSCCBitSet -- the set of model states that are trivial \Psi BSCCs										*/
/*																													*/
/*	1.1 numberOfTrivialBSCCs -- number of elements in pTrivialBSCCBitSet											*/
/*																													*/
/*2. ppNonTrivBSCCBitSets -- the array of all non-trivial \Psi BSCCs												*/
/*																													*/
/*	2.1 numberOfNonTrivBSCCs -- the size of ppNonTrivBSCCBitSets													*/
/*																													*/
/*3. ppBSCCExistsReachableStates -- the array of bitsets each element of an array corresponds to a bitset.			*/
/*									Here each bitset lists set of states from which the corresponding BSCC			*/
/*									is always reachable. I.e. starting in the state we only end up in that BSCC.	*/
/*																													*/
/*	3.1 The first numberOfTrivialBSCCs bitset pointers correspond to the trivial									*/
/*		BSCCs. They are assumed to be ordered in the same way their states are										*/
/*		ordered in pTrivialBSCCBitSet.																				*/
/*																													*/
/*	3.2 The next numberOfNonTrivBSCCs bitset pointer correspond to the non-trivial									*/
/*		BSCCs. They are assumed to be ordered in the same way as in ppNonTrivBSCCBitSets.							*/
/*																													*/
/*4. ppBSCCAlwaysReachableStates -- the array of bitsets each element of an array corresponds to a bitset.			*/
/*									Here each bitset lists set of states from which the corresponding BSCC			*/
/*									can be reached. I.e. starting in the state we can end up in that BSCC.			*/
/*																													*/
/*	4.1 The first numberOfTrivialBSCCs bitset pointers correspond to the trivial									*/
/*		BSCCs. They are assumed to be ordered in the same way their states are										*/
/*		ordered in pTrivialBSCCBitSet.																				*/
/*																													*/
/*	4.2 The next numberOfNonTrivBSCCs bitset pointer correspond to the non-trivial									*/
/*		BSCCs. They are assumed to be ordered in the same way as in ppNonTrivBSCCBitSets.							*/
/*																													*/
/*5. pReachableTrivialBSCCBitSet -- the set of trivial BSCC indexes that are reachable from the given				*/
/*									initial states. These indexes correspond to BSCCs in pTrivialBSCCBitSet.		*/
/*									I.e. if we take all of the enabled states in pTrivialBSCCBitSet in their		*/
/*									natural order, and number them, these indexes will be possible elements			*/
/*									of pReachableTrivialBSCCBitSet.													*/
/*																													*/
/*	5.1 numberOfReachableTrivialBSCCs -- is the number of elements in pReachableTrivialBSCCBitSet					*/
/*																													*/
/*6. pReachableNonTrivBSCCBitSet -- the set of non-trivial BSCC indexes that are reachable from the given			*/
/*									initial states.  These indexes correspond to BSCCs in ppNonTrivBSCCBitSets.		*/
/*									I.e. it is the index of the bitset of ppNonTrivBSCCBitSets which contains		*/
/*									the states of this \Psi BSCC.													*/
/*																													*/
/*7. alwaysReachableBSCCIndex -- the index in the range
                        [0, (numberOfTrivialBSCCs + numberOfNonTrivBSCCs) )  */
/*                      it is set if from this state we always go to one
                        particular BSCC                                      */
/********************************************************************************************************************/

/**
* This function is an implementation of the steadyState function from
* the PhD Thesis of Ivan S. Zapreev. The functions do the model
* checking of unbounded-until operator on the provided CTMC and uses the
* runtime-simulation settings from the "simulation.h"
* @param pStateSpace the sparse matrix of the embedded DTMC with good and bad
*			states made absorbing
* @param confidence this is is confidence we should use for simulations.
*			Note that since the formulas can be nested one can not use the
*			overall confidence set in simulation.h
* @param pPsiBitSet the bitsets containing the Psi states for S<>p(Psi)
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @param isSimOneInitState TRUE if we simulate one initial state only
* @param initial_state the initial state index if isSimOneInitState == TRUE
* @param pExistsUnbUntilCTMC the function pointer that gives a template for the function
*                            that will do the reachability search. In general it is for the exists
*                            unbounded until operator (\exists \Phi U \Psi). So we look for the states
*                            from which there is a path to \Psi states via \Phi states.
* @param pAlwaysUnbUntilCTMC the function pointer that gives a template for the function
*                            that will do the reachability search. In general it is for the always
*                            unbounded until operator (\always \Phi U \Psi). So we look for the states
*                            from which all the paths go to \Psi states via \Phi states.
* @param pFAllowedBSCCSearchCTMC the function that is used for searching the BSCCs that contain Psi states.
*
* NOTE: Below we have pointers to the pointers, these are going to be the return values
*
* @param ppYesBitSet this bitset is going to be filled with the states which
*			satisfy the formula
* @param ppNoBitSet this bitset is going to be filled with the states which
*			do not satisfy the formula
* @param ppCiLeftBorders will contain the left conf. int. borders
* @param ppCiRightBorders will contain the right conf. int. borders
* @param pResultSize the pointer to the size of arrays *ppCiLeftBorders and *ppCiRightBorders
* @param pMaxNumUsedObserv the return-value variable that will store the maximum number
* (over all initial states) of used observations
*/
void modelCheckSteadyStatePureCTMC(sparse * pStateSpace,
                const double * pCTMCRowSums,
									const double confidence, const bitset * pPsiBitSet,
									bitset ** ppYesBitsetResult, bitset ** ppNoBitsetResult,
									double ** ppProbCILeftBorder, double ** ppProbCIRightBorder,
									int * pResultSize, const int comparator, const double prob_bound,
                const int initial_state, const BOOL isSimOneInitState_local,
									const TPFunctExistUnbUntil pExistsUnbUntilCTMC,
									const TPFunctAlwaysUnbUntil pAlwaysUnbUntilCTMC,
									const TPFunctAllowedBSCCSearchCTMC pFAllowedBSCCSearchCTMC,
									unsigned int * pMaxNumUsedObserv ) {
	/*Some technical variables*/
	int i;
	/* The number of BSCCs containing Psi states */
	int numberOfBSCCs;
	/* The number of non-trivial BSCCs containing Psi states */
	int numberOfNonTrivBSCCs = 0;
	/* The average number of states in non-trivial BSCCs containing Psi states */
	int numberOfNonTrivBSCCStates = 0;
	/* The number of trivial BSCCs containing Psi states */
	int numberOfTrivialBSCCs = 0;
	/* The structure, that stores the set of states belonging to non-trivial BSCC with array index i */
	bitset ** ppNonTrivBSCCBitSets = NULL;
	/* The bitset indicating the Psi states belonging to trivial BSCCs */
	bitset * pTrivialBSCCBitSet = NULL;
	/* This vector will store the exit rates for all the states since */
	/* we do pure simulations but not only simulate \Psi BSCCs */
	double * pExitRatesOfAllStates = NULL;
	/* The state space dimension */
        const int n_states = mtx_rows(pStateSpace);

	/* The structure, which stores the samples */
	PTSampleVecSteady *pSampleVecSteady = NULL;

	/* These are the arrays that store BSCC reachability information, they are */
	/* indexed by the naturaly ordered trivial BSCCs and then by the non-trivial */
	/* BSCCs ordered as in ppNonTrivBSCCBitSets */
	bitset ** ppBSCCAlwaysReachableStates = NULL;
	bitset ** ppBSCCExistsReachableStates = NULL;

	/* Compute the transient states, i.e. states that do not belong to any BSCC */
	/* These are needed for the UU simulation, i.e. reaching BSCCs */
	bitset * pTransientBitSet = getTransientStates( pStateSpace, pFAllowedBSCCSearchCTMC );
	/* printf("TRANSIENT STATES: "); print_bitset_states( pTransientBitSet ); printf("\n"); */

	/* The bitset that holds the possible initial states */
	bitset *pValidInitialStatesBitSet = get_new_bitset( n_states );
	/* All states are possible initial states */
	fill_bitset_one( pValidInitialStatesBitSet );

	/* Find out the BSCCs that contain \Psi states including the trivial ones */
	ppNonTrivBSCCBitSets = getTrivialAndNonTrivialPsiBSCC( pStateSpace, pPsiBitSet, NULL,
															&pTrivialBSCCBitSet, &numberOfNonTrivBSCCs,
															&numberOfTrivialBSCCs, &numberOfBSCCs,
															pFAllowedBSCCSearchCTMC );
	/* printf("NUMBER OF PSI BSCCs: %d\n", numberOfBSCCs);
	printf("NUMBER OF TRIVIAL PSI BSCCs: %d\n", numberOfTrivialBSCCs);
	printf("TRIVIAL PSI BSCCs: "); print_bitset_states( pTrivialBSCCBitSet ); printf("\n");
	printf("NUMBER OF NON-TRIVIAL PSI BSCCs: %d\n", numberOfNonTrivBSCCs);
	for( i=0; i < numberOfNonTrivBSCCs; i++ ) {
		printf("NON-TRIVIAL PSI BSCCs [%d]: ", i); print_bitset_states( ppNonTrivBSCCBitSets[i] ); printf("\n");
	} */

	/* Compute the sets of states from which we always, sometimes and never go some BSCC with \Phi states */
        ppBSCCExistsReachableStates = (bitset **) calloc((size_t) numberOfBSCCs,
                        sizeof(bitset *));
        ppBSCCAlwaysReachableStates = (bitset **) calloc((size_t) numberOfBSCCs,
                        sizeof(bitset *));
	getReachabilityTrivialAndNonTrivialPsiBSCC( pStateSpace, numberOfBSCCs, pTrivialBSCCBitSet,
												ppNonTrivBSCCBitSets, numberOfNonTrivBSCCs,
												ppBSCCExistsReachableStates, ppBSCCAlwaysReachableStates,
												pExistsUnbUntilCTMC, pAlwaysUnbUntilCTMC );
	/* for( i=0; i < numberOfBSCCs; i++ ) {
		printf("EXIST A PATH TO BSCC [%d] from: ", i); print_bitset_states( ppBSCCExistsReachableStates[i] ); printf("\n");
		printf("ALL PATHs TO BSCC [%d] from: ", i); print_bitset_states( ppBSCCAlwaysReachableStates[i] ); printf("\n");
	 } */

	/* Obtain the embedded DTMC and the vector of exit rates */
	/* WARNING: The self loops on states are not modified. Thus the proper way */
	/* to identify an absorbing state is to make sure that there are */
	/* no non-diagonal elements in the corresponding row. */
	/* printf("MAKING EMBEDDED DTMC\n"); */
	pExitRatesOfAllStates = make_embedded_dtmc_out_of_rate_mtx_all( pStateSpace, pCTMCRowSums );
        if ( NULL == pExitRatesOfAllStates ) {
                exit(err_macro_17(err_CALLBY, "modelCheckSteadyStatePureCTMC("
                        "%p[%dx%d],%p,%g,%p[%d],%p,%p,%p,%p,%p,%d,%g,%d,%d,_,_,"
                        "_,%p)", (void*) pStateSpace, mtx_rows(pStateSpace),
                        mtx_cols(pStateSpace), (const void *) pCTMCRowSums,
                        confidence, (const void *) pPsiBitSet,
                        bitset_size(pPsiBitSet), (void *) ppYesBitsetResult,
                        (void *) ppNoBitsetResult, (void *) ppProbCILeftBorder,
                        (void *) ppProbCIRightBorder, (void *) pResultSize,
                        comparator, prob_bound, initial_state,
                        isSimOneInitState_local,
                        (void *) pMaxNumUsedObserv, EXIT_FAILURE));
        }

	/* Prepare the steady state simulation samples vector for */
	/* each BSCC and choose regeneration states for each BSCC */
	/* Also count the number of the non-trivial BSCCs */
	/* printf("PREPARING SS SIMULATION SAMPLES\n"); */
	pSampleVecSteady = prepareSSSimulationSample( pStateSpace, numberOfNonTrivBSCCs, ppNonTrivBSCCBitSets, &numberOfNonTrivBSCCStates );

    /* Calculate the average number of states in a BSCC */
	if( numberOfNonTrivBSCCs > 0 ){
		numberOfNonTrivBSCCStates = numberOfNonTrivBSCCStates / numberOfNonTrivBSCCs;
	}

	/* Call the common model checking procedure algorithm */
	/* printf("CALLING THE MAIN SIMULATION PROCEDURE\n"); */
    modelCheckStatesCommon( pStateSpace, confidence, ppYesBitsetResult, ppNoBitsetResult, ppProbCILeftBorder,
							ppProbCIRightBorder, pResultSize, comparator, prob_bound, initial_state,
                        isSimOneInitState_local, pValidInitialStatesBitSet,
                        pMaxNumUsedObserv,
							modelCheckOneStateSSPureCTMC, pPsiBitSet, ppBSCCExistsReachableStates,
							ppBSCCAlwaysReachableStates, pExitRatesOfAllStates, pSampleVecSteady,
							numberOfNonTrivBSCCs, numberOfTrivialBSCCs, numberOfNonTrivBSCCStates,
							pTransientBitSet, ppNonTrivBSCCBitSets, pTrivialBSCCBitSet );

	/* Restore the state space */
	/* NOTE: The diagonal elements do not have to be restored since they were not changed */
	/* printf("REVERTING THE EMBEDDED DTMC BACK INTO A CTMC\n"); */
        if ( err_state_iserror(restore_rate_mtx_out_of_embedded_dtmc_all(
                                pStateSpace, pExitRatesOfAllStates)) )
        {
                exit(err_macro_17(err_CALLBY, "modelCheckSteadyStatePureCTMC("
                        "%p[%dx%d],%p,%g,%p[%d],%p,%p,%p,%p,%p,%d,%g,%d,%d,_,_,"
                        "_,%p)", (void*) pStateSpace, mtx_rows(pStateSpace),
                        mtx_cols(pStateSpace), (const void *) pCTMCRowSums,
                        confidence, (const void *) pPsiBitSet,
                        bitset_size(pPsiBitSet), (void *) ppYesBitsetResult,
                        (void *) ppNoBitsetResult, (void *) ppProbCILeftBorder,
                        (void *) ppProbCIRightBorder, (void *) pResultSize,
                        comparator, prob_bound, initial_state,
                        isSimOneInitState_local,
                        (void *) pMaxNumUsedObserv, EXIT_FAILURE));
        }

	/* Clean the memory: BSCC's bitsets, samples and etc */
	/* printf("CLEANING MEMORY\n"); */
	free( pExitRatesOfAllStates ); pExitRatesOfAllStates = NULL;
	free_bitset( pTrivialBSCCBitSet ); pTrivialBSCCBitSet = NULL;
	free_bitset( pValidInitialStatesBitSet ); pValidInitialStatesBitSet = NULL;

	for( i = 0; i < numberOfNonTrivBSCCs; i++ ){
		free_bitset( ppNonTrivBSCCBitSets[i] );
	}
	free( ppNonTrivBSCCBitSets ); ppNonTrivBSCCBitSets = NULL;

	for( i = 0; i < numberOfBSCCs; i++ ){
		free_bitset( ppBSCCExistsReachableStates[i] );
		free_bitset( ppBSCCAlwaysReachableStates[i] );
	}
	free( ppBSCCExistsReachableStates ); ppBSCCExistsReachableStates = NULL;
	free( ppBSCCAlwaysReachableStates ); ppBSCCAlwaysReachableStates = NULL;

	for( i = 0; i < numberOfNonTrivBSCCs; i++){
		freeSampleVectorSteady( pSampleVecSteady[i] );
	}
	free( pSampleVecSteady ); pSampleVecSteady = NULL;

	free_bitset( pTransientBitSet ); pTransientBitSet = NULL;
}
