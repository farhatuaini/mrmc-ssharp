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
*	Source description: This is a source file for the simulation engine
*	here we intend to define the common functions and data structures.
*/

#include "simulation_common.h"

#include "simulation.h"
#include "simulation_utils.h"

#include <gsl/gsl_cdf.h>
#include <math.h>

/****************************************************************************/
/********THE COMMON FUNCTION FOR INVOCATION THE SIMULATION PROCEDURE ********/
/****************************************************************************/

/**
* This function is responsible for invoking the model procedure either
* for one or for all initial states the invoked finction for model checking
* one initial state is passed as an argument pFMCOneState
* @param pStateSpace the state-space matrix
* @param confidence this is is confidence we should use for simulations.
*	Note that since the formulas can be nested one can not use the
*	overall confidence set in simulation.h
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
* @param pValidInitialStates the bitset containing all allowed initial
* @param pMaxNumUsedObserv the return-value variable that will store the maximum number
*			(over all initial states) of used observations
* @param pFMCOneState the pointer to a function which model checks one initial state
* WARNING: The parameter pFMCOneState should be the last one before "..."
* @param ... the extra arguments needed for a specific model-checking procedure
*/
void modelCheckStatesCommon( const sparse * pStateSpace, const double confidence, bitset ** ppYesBitSet,
				bitset ** ppNoBitSet, double ** ppCiLeftBorders, double ** ppCiRightBorders,
				int * pResultSize, const int comparator, const double prob_bound,
                                const int initial_state,
                                const BOOL isSimOneInitState_local,
				const bitset * pValidInitialStates, unsigned int *pMaxNumUsedObserv,
				TPFunctMCOneState pFMCOneState, ... ){
	/* Will contain the extra arguments of the function pFMCOneState */
	va_list fMCOneStateArguments;
	/* Obtain the width of the indifference region */
	const double indiff_width = getSupIndifferenceWidth();
	/* Will store the number of used observations */
	unsigned int numUsedObserv = 0;

	/* Allocate the resulting bitsets */
        const int state_space_size = mtx_rows(pStateSpace);
	*ppYesBitSet = get_new_bitset( state_space_size );
	*ppNoBitSet = get_new_bitset( state_space_size );

	/* Reset the maximum number of used observations. */
	*pMaxNumUsedObserv = 0;

	if( ( *ppYesBitSet != NULL ) && ( *ppYesBitSet != NULL ) ){
		/* Depending on whether there is one initial state only or not. */
                if ( isSimOneInitState_local ) {
			/* Check for the valid initial state */
			IF_SAFETY( ( initial_state >= 0) || ( initial_state < state_space_size ) )
				/* If it is a valid state */
				if( get_bit_val( pValidInitialStates, initial_state ) ){
					/* Allocate the return data */
					*pResultSize = 1;
                                        *ppCiLeftBorders = (double *) calloc(
                                                (size_t) 1, sizeof(double));
                                        *ppCiRightBorders = (double *) calloc(
                                                (size_t) 1, sizeof(double));

					/* Obtain the extra arguments of the function pFMCOneState */
					va_start( fMCOneStateArguments, pFMCOneState );

					/* Reset the observations counter for one initial state */
					numUsedObserv = 0;

					/* Model check the state */
					pFMCOneState( pStateSpace, initial_state, indiff_width, confidence,
							comparator, prob_bound, 0, *ppCiLeftBorders, *ppCiRightBorders,
							*ppYesBitSet, *ppNoBitSet, &numUsedObserv, fMCOneStateArguments );

					/* Since there is just one initial state: maxNumUsedObserv = numUsedObserv */
					*pMaxNumUsedObserv = numUsedObserv;

					/* Complete the argument list usage */
					va_end( fMCOneStateArguments );
				}else{
					*pResultSize = 0;
				}
			ELSE_SAFETY
				printf("ERROR: An ivalid initial state %d, should be >= %d and <= %d.\n",
					initial_state+1, 1, state_space_size );
                                exit(EXIT_FAILURE);
			ENDIF_SAFETY
		} else {
			int curr_state = -1;

			/* Allocate the return data */
			*pResultSize = state_space_size;
                        *ppCiLeftBorders = (double *) calloc(
                                (size_t) state_space_size, sizeof(double));
                        *ppCiRightBorders = (double *) calloc(
                                (size_t) state_space_size, sizeof(double));

			/* Work out all the transient states */
			while( ( curr_state = get_idx_next_non_zero( pValidInitialStates, curr_state  )) != -1 ){
				/* Obtain the extra arguments of the function pFMCOneState */
				va_start( fMCOneStateArguments, pFMCOneState );

				/* Reset the observations counter for one initial state */
				numUsedObserv = 0;

				/* Model check the state */
				pFMCOneState( pStateSpace, curr_state, indiff_width, confidence, comparator,
						prob_bound, curr_state, *ppCiLeftBorders, *ppCiRightBorders,
						*ppYesBitSet, *ppNoBitSet, &numUsedObserv, fMCOneStateArguments );

				/* Take the muximum between the number of observations used now and before*/
				if( numUsedObserv > *pMaxNumUsedObserv ){
					*pMaxNumUsedObserv = numUsedObserv;
				}

				/* printf("The number of observations used to simulate state %d is %u.\n", curr_state+1, numUsedObserv ); */

				/* Complete the argument list usage */
				va_end( fMCOneStateArguments );
			}
		}
	} else {
		printf("ERROR: Failed to allocate data structures for the simulation engine.\n");
                exit(EXIT_FAILURE);
	}
}

/****************************************************************************/
/********THE COMMON FUNCTION FOR SIMULATING UNBOUNDED UNTIL CTMC/DTMC********/
/****************************************************************************/

/**
* Model checks the unbounded-until operator for one initial state only.
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
* @param pTransAndGoodStatesArguments contain the two bitsets in the following order:
*		1. pTransientStates the bitsets containing all the transient states
*		2. pGoodStates the bitsets containing all the good absorbing states
*/
void modelCheckOneStateUUCommon( const sparse* pStateSpace, const int initial_state, const double indiff_width,
					const double confidence, const int comparator, const double prob_bound,
					const int arr_index, double *pCiLeftBorders, double *pCiRightBorders,
					bitset * pYesBitSet, bitset * pNoBitSet, unsigned int *pNumUsedObserv,
					va_list pTransAndGoodStatesArguments ){
	/* Contain the Transient states and also the Good absorbing states respectively */
	bitset * pTransientStates = va_arg( pTransAndGoodStatesArguments, bitset * );
	bitset * pGoodStates = va_arg( pTransAndGoodStatesArguments, bitset * );

	TV_LOGIC mc_result = TVL_NN;
	double leftBorderOne = 0.0, rightBorderOne = 1.0, leftBorderTwo = 0.0, rightBorderTwo = 1.0;

	/* Compute the general confidence zeta value: */
	const double gen_conf_zeta = gsl_cdf_ugaussian_Pinv( sqrt( confidence ) );
	/* NOTE: The confidence we have is used to derive z for the left and the */
	/* right border BUT the TRUE/FALSE answer is based only on one border of */
	/* the conf. int., i.e. we could use gsl_cdf_ugaussian_Pinv( confidence ) */
	/* for model checking, but then again we have a different confidence for the */
	/* produced conf. int. this is bad because then the user does not know that. */

	/* Initialize the sample size parameters */
	int sample_size = getSimMinSampleSize();
	const int sample_size_step = getSimSampleSizeStep();
	const int max_sample_size = getSimMaxSampleSize();

	/* Initialize the simulation depth parameters */
	int simulation_depth = getSimMinSimulationDepth();
	const int simulation_depth_step = getSimSimulationDepthStep();
	const int max_simulation_depth = getSimMaxSimulationDepth();

	/* The odd/even iteration and invalid conf.int. indicators */
	BOOL isOdd = TRUE, isInvalidCI = FALSE;

	/* Allocate and simulate the initial sample vectors */
	PTSampleVecUntil pSampleVecUntilOne = allocateSampleVectorUntil( sample_size, initial_state );
	PTSampleVecUntil pSampleVecUntilTwo = allocateSampleVectorUntil( sample_size, initial_state );

	/* Simulate observations in depth */
	simulateSampleVectorUnbUntilDTMC( pStateSpace, pSampleVecUntilOne, sample_size, simulation_depth,
										pGoodStates, pTransientStates );
	simulateSampleVectorUnbUntilDTMC( pStateSpace, pSampleVecUntilTwo, sample_size, simulation_depth,
										pGoodStates, pTransientStates );

	/* Perform the main simulation cycle */
	do{
		/* Choose what to increment  */
		if( ( isOdd || ( sample_size >= max_sample_size ) ) &&
			( simulation_depth < max_simulation_depth ) && !isInvalidCI ){
			/*One could also check for the following: */
			/*( pSampleVecUntilOne->sum_trans + pSampleVecUntilTwo->sum_trans ) > 0 */
			/*but it seems to cause infinite looping some times and*/
			/*there is not real need for checking this condition, because */
			/*in this case we mostly skip the iteration any ways.*/
			increment( &simulation_depth, simulation_depth_step, max_simulation_depth );
		} else {
			increment( &sample_size, sample_size_step, max_sample_size );
		}

		/* Simulate the samples */
		simulateSampleVectorUnbUntilDTMC( pStateSpace, pSampleVecUntilOne, sample_size, simulation_depth,
											pGoodStates, pTransientStates );
		simulateSampleVectorUnbUntilDTMC( pStateSpace, pSampleVecUntilTwo, sample_size, simulation_depth,
											pGoodStates, pTransientStates );

		/* Compute the Conf Int borders */
		computeBordersUU( gen_conf_zeta, &leftBorderOne, &rightBorderOne, pSampleVecUntilOne, AGRESTI_COULL_CONF_INT );
		computeBordersUU( gen_conf_zeta, &leftBorderTwo, &rightBorderTwo, pSampleVecUntilTwo, AGRESTI_COULL_CONF_INT );

		/* Test the produced conf. int. for having a non-zero intersection */
		/* Note that leftBorderOne <= rightBorderOne and leftBorderTwo <= rightBorderTwo */
		isInvalidCI = ( rightBorderOne < leftBorderTwo ) || ( rightBorderTwo < leftBorderOne);

		/* Check the probability constraint */
		/* WARNING: We can not do any sort of sample swapping here! */
		if( !isInvalidCI ){
			mc_result = checkBoundVSConfInt( comparator, prob_bound, leftBorderOne, rightBorderTwo, indiff_width );
		}
		isOdd = ! isOdd;
	}while( ( mc_result == TVL_NN ) &&
			( ( ( simulation_depth < max_simulation_depth ) && ! isInvalidCI ) ||
				( sample_size < max_sample_size ) ) );

	/* Update the pCiLeftBorders and pCiRightBorders */
	pCiLeftBorders[arr_index] = leftBorderOne;
	pCiRightBorders[arr_index] = rightBorderTwo;

	/* printf("Result: %i, borders are: [%lf, %lf], prob bound %d %lf.\n", mc_result, */
	/*	pCiLeftBorders[arr_index], pCiRightBorders[arr_index], comparator, prob_bound); */

	/* Update the pYesBitSet or pNoBitSet if we have a definite answer */
	markYesNoSetEntree(mc_result, initial_state, pYesBitSet, pNoBitSet);

	/* printf("mc_result = %d, sample_size = %d and simulation_depth = %d\n", */
	/*	mc_result, sample_size, simulation_depth); */
	/* printSampleVectorUntil(pSampleVecUntilOne); */
	/* printSampleVectorUntil(pSampleVecUntilTwo); */

	/* Report the number of states visited while this simulation run */
	*pNumUsedObserv = ( (PTSampleVec) pSampleVecUntilOne )->num_visited_states +
				( (PTSampleVec) pSampleVecUntilTwo )->num_visited_states;

	/* Free the sample vectors */
	freeSampleVectorUntil( pSampleVecUntilOne );
	freeSampleVectorUntil( pSampleVecUntilTwo );
}


/**
 * This function allows to simulate one step and compute the next conf. int. for the reachbility
 * problem. We use this method when model checking steady-state operator with pure simulations.
 * Note: When a confidence interval is produced the method automatically tries to increase the
 * sample size this is done until we either get a good conf. int. or we run out of samples.
 * @param pStateSpace the state space
 * @param pGoodStates the states we want to reach, these are states of some \Psi BSCC.
 * @param pTransientStates the set of all transient states in our model.
 * @param zeta the zeta value computed with proper confidence.
 * @param pSampleVecUntilOne the left sample we extend and use
 * @param pSampleVecUntilTwo the right sample we extend and use
 * @param sample_size the provided (possibly new) sample size
 * @param sample_size_step the size by which we increase sample_size in case we get an invalid conf. int.
 * @param max_sample_size the size until which we can increase sample_size when we get an invalid conf. int.
 * @param simulation_depth the provided (possibly new) simulation depth size
 * @return true if we have gotten an invalid conf. int. that we could not fix even using the max
 *			number of samples. in this case the model checking procedure has to be terminated.
 */
static BOOL simulateAndComputeOneReachabilityConfInt(const sparse * const pStateSpace, const bitset * const pGoodStates,
												const bitset * const pTransientStates, const double zeta,
												PTSampleVecUntil pSampleVecUntilOne, PTSampleVecUntil pSampleVecUntilTwo,
												double * pLeftConfIntBorder, double * pRightConfIntBorder,
												int sample_size, const int sample_size_step, const int max_sample_size,
												int simulation_depth ) {
	/* The temporary holders for the conf. int. values for each of the vectors */
	double leftBorderOne = 0.0, rightBorderOne = 1.0, leftBorderTwo = 0.0, rightBorderTwo = 1.0;
	/* The indicator of the invalidity of the conf. int. */
	BOOL isInvalidCI = FALSE;

	/* Do simulation and compute the conf int borders, in case the conf. int. is invalid, then */
	/* we have to repeat iterations, increasing the sample size, until it becomes valid again. */
	do {
		/* Simulate the samples */
		simulateSampleVectorUnbUntilDTMC( pStateSpace, pSampleVecUntilOne, sample_size, simulation_depth,
											pGoodStates, pTransientStates );
		simulateSampleVectorUnbUntilDTMC( pStateSpace, pSampleVecUntilTwo, sample_size, simulation_depth,
											pGoodStates, pTransientStates );

		/* Compute the conf. int. borders */
		computeBordersUU( zeta, &leftBorderOne, &rightBorderOne, pSampleVecUntilOne, AGRESTI_COULL_CONF_INT );
		computeBordersUU( zeta, &leftBorderTwo, &rightBorderTwo, pSampleVecUntilTwo, AGRESTI_COULL_CONF_INT );

		/* Test the produced conf. int. for having a non-zero intersection */
		/* Note that leftBorderOne <= rightBorderOne and leftBorderTwo <= rightBorderTwo */
		if( ( isInvalidCI = ( rightBorderOne < leftBorderTwo ) || ( rightBorderTwo < leftBorderOne) ) ) {
			/* NOTE: Here we increment the local variable only, this means that in the calling method */
			/* the value of sample_size will remain the same, we did it to avoid large samples in case */
			/* we do reachability for one BSCC then we get an invalid interval and then we have to */
			/* increase the sample size until the moment the interval is valid again. The problem is*/
			/* with this approach all of the following BSCCs will have to have larger sample size */
			/* even if it is not needed. The drawback of our current approach is that sometimes */
			/* we can have all of the sample sizes set to maximum but we will not know that in the */
			/* calling method. In this case we will do iterations that will be empty, but still */
			/* the algorithm should eventually terminate, so we do not worry for this overhead */
			increment( &sample_size, sample_size_step, max_sample_size );
		}
	}while( isInvalidCI );

	/* Store the obtained values for the conf. int. if it is valid, otherwise keep the old one */
	if( ! isInvalidCI ) {
		(* pLeftConfIntBorder) = leftBorderOne;
		(* pRightConfIntBorder) = rightBorderTwo;
	}

	return isInvalidCI;
}

/**
 * This function simulated the unbounded reachability samples for the model checking
 * steady-state operator with a pure simulation method.
 * NOTE: The increase of sample size and simulation depth are guided by the common values
 * of (*pSampleSizeUU) and (*pSimDepthUU) at the same time, the actual size (depth) of
 * the samples can be larger (smaller) but this should be fine because the depth is defined
 * by the time we reach the BSCC with \Psi states and the sample size can be increased if we
 * get an invalid conf. int. for the reachability problem. In the latter case if the size
 * (depth) is larger (smaller) then we will just skip increase of the sample size
 * (increase the simulation depth).
 * @param pStateSpace the state space
 * @param pTransientBitSet the transient states of the model
 * @param ppNonTrivBSCCBitSets the non-trivial BSCC's (their states) with \Psi states in them
 * @param pTrivialBSCCBitSet the bitset that contains trivial \PSI BSCCs (single-state BSCCs)
 * @param isOdd indicates whether it is an odd iteration or not
 * @param pReachableTrivialBSCCBitSet contains the indices of reachable trivial
 *                                    \Psi BSCCs
 * @param pReachableNonTrivBSCCBitSet contains the indices of reachable
 *                                    non-trivial BSCCs with \Psi states.
 * @param pSampleVecUntilOneArray the array of the number of reachable BSCCs with \Psi states
 *                                which stores samples for computing the left conf. int. borders
 * @param pSampleVecUntilTwoArray the array of the number of reachable BSCCs with \Psi states
 *                                which stores samples for computing the right conf. int. borders
 * @param pLeftBorderUU the array with the left conf. int. borders that have to be computed here.
 *                      This array has the size == number of reachable \PSI BSCCs
 * @param pRightBorderUU the array with the right conf. int. borders that have to be computed here.
 *                       This array has the size == number of reachable \PSI BSCCs
 * @param pSampleSizeUU the pointer to the value that stores the current sample size common
 *						for each sample in pSampleVecUntilOneArray and pSampleVecUntilTwoArray
 * @param max_sample_size_uu the maximum allowed value for (*pSampleSizeUU)
 * @param sample_size_step_uu the step by which we increase (*pSampleSizeUU)
 * @param pSimDepthUU the pointer to the variable that stores current simulation depth common
 *                    for each sample in pSampleVecUntilOneArray and pSampleVecUntilTwoArray
 * @param max_sim_depth_uu the maximum allowed value for (*pSimDepthUU)
 * @param sim_depth_step_uu the step by which we increase (*pSimDepthUU)
 * @return false if we could not make one of the reachability conf. int. to be non-invalid
 */
BOOL simulateAndComputeAllReachabilityConfInt( const sparse * const pStateSpace, const bitset * const pTransientBitSet,
                const bitset * const pTrivialBSCCBitSet,
                const bitset * const * const ppNonTrivBSCCBitSets,
												const BOOL isOdd, const bitset * const pReachableTrivialBSCCBitSet,
												const bitset * const pReachableNonTrivBSCCBitSet, PTSampleVecUntil * pSampleVecUntilOneArray,
												PTSampleVecUntil * pSampleVecUntilTwoArray, double * pLeftBorderUU, double * pRightBorderUU,
												int * pSampleSizeUU, const int max_sample_size_uu, const int sample_size_step_uu,
												int * pSimDepthUU, const int max_sim_depth_uu, const int sim_depth_step_uu, const double zeta ) {
	/* Counters and iterator varibale */
        int state_index_local, bscc_index, conf_int_sample_index;
	/* The temporary bitset */
        bitset * pTmpBitSet = get_new_bitset(mtx_rows(pStateSpace));
	/* The indicator of the invalidity of the conf. int. */
	BOOL isInvalidCI = FALSE;

	/* Increment the sample size/simulation depth, Choose what to increment */
	if( ( isOdd || ( (*pSampleSizeUU) >= max_sample_size_uu ) ) && ( (*pSimDepthUU) < max_sim_depth_uu ) ){
		increment( pSimDepthUU, sim_depth_step_uu, max_sim_depth_uu );
	} else {
		increment( pSampleSizeUU, sample_size_step_uu, max_sample_size_uu );
	}

	/* First process the trivial BSCCs, i.e. get the simulations for them and */
	/* then compute the conf. int. Terminate if we got an invalid conf. int. */
        state_index_local = state_index_NONE;
        conf_int_sample_index = 0;
        bscc_index = 0;
	/* Go through the state indices that correspond to trivial \Psi BSCCs */
	/* If isInvalidCI becomes true then it means that we have used all available */
	/* sample sizes for the conf. int. for reaching some BSCCs but this interval */
	/* remains to be invalid, this means that any further simulations are worthless */
        while ( (state_index_local = get_idx_next_non_zero(pTrivialBSCCBitSet,
                                        state_index_local)) != state_index_NONE
                                && ! isInvalidCI )
        {
		if( get_bit_val( pReachableTrivialBSCCBitSet, bscc_index ) ) {
			/* If this BSCC is reachable then create the set of states that correspond to it */
                        set_bit_val(pTmpBitSet, state_index_local, BIT_ON);

			/* So simulation here and compute the conf. int. for the reachability problem */
			isInvalidCI = simulateAndComputeOneReachabilityConfInt( pStateSpace, pTmpBitSet, pTransientBitSet, zeta,
																	pSampleVecUntilOneArray[ conf_int_sample_index ],
																	pSampleVecUntilTwoArray[ conf_int_sample_index ],
																	(& pLeftBorderUU[ conf_int_sample_index ] ),
																	(& pRightBorderUU[ conf_int_sample_index ] ), (*pSampleSizeUU),
																	sample_size_step_uu, max_sample_size_uu, (*pSimDepthUU) );
			/* Clean this BSCC states set */
                        set_bit_val(pTmpBitSet, state_index_local, BIT_OFF);
		}
		conf_int_sample_index++;
		bscc_index++;
	}

	/* Second process the non-trivial BSCCs, i.e. get the simulations for them and */
	/* then compute the conf. int. Terminate if we got an invalid conf. int. */
	bscc_index = -1;
	while( ( ( bscc_index = get_idx_next_non_zero( pReachableNonTrivBSCCBitSet, bscc_index ) ) != -1 ) && !isInvalidCI ) {
		/* So simulation here and compute the conf. int. for the reachability problem */
		isInvalidCI = simulateAndComputeOneReachabilityConfInt( pStateSpace, ppNonTrivBSCCBitSets[bscc_index],
																pTransientBitSet, zeta,
																pSampleVecUntilOneArray[ conf_int_sample_index ],
																pSampleVecUntilTwoArray[ conf_int_sample_index ],
																(& pLeftBorderUU[ conf_int_sample_index ] ),
																(& pRightBorderUU[ conf_int_sample_index ] ), (*pSampleSizeUU),
																sample_size_step_uu, max_sample_size_uu, (*pSimDepthUU) );
		conf_int_sample_index++;
	}

	/* Free the unneeded memory */
	free_bitset( pTmpBitSet ); pTmpBitSet = NULL;

	return isInvalidCI;
}

/**
 * This method allows to prepare the sample vectors for each of the BSCCs
 * that contain \Psi states and thus might need to be simulated.
 * @param pStateSpace the actual state space, i.e. the sparse matrix representing the model.
 * @param numberOfNonTrivBSCCs the number of non-trivial BSCCs that contain \Psi states
 * @param ppNonTrivBSCCBitSets the set of states belonging to non-trivial BSCC with array index i
 * @param pNumberOfNonTrivBSCCStates the pointer to the variable storing the number of the non-
 *									trivial BSCCs, this is an implicit return parameter
 * @return the vector of sample for each non-trivial BSCC initialised with the regeneration state
 */
PTSampleVecSteady * prepareSSSimulationSample(const sparse * pStateSpace, const int numberOfNonTrivBSCCs,
												bitset ** const ppNonTrivBSCCBitSets, int * pNumberOfNonTrivBSCCStates ) {
	/*Simple counter variables for the "for" loops*/
	int i,j;
	/*The state space dimension, i.e. the number of rows(columns) in the space matrix*/
        const int n_states = mtx_rows(pStateSpace);
	/* The multiplier for the choice of the sample-based regeneration point */
	const int sampleBSCCDimensionMultiplier = getSampleBSCCDimensionMultiplier();
	/* The number of states in the BSCC currently handled */
	int numberOfBSCCStates;
	/*This variable stores the initial (regeneration) state for every non trivial BSCC*/
	int initial_state_sim;
	/* The array that holds the number state occurrences while simulating to fix the regeneration state */
	int * pSimStateOccurrence;
	/*The counter for the number of non-trivial BSCCs*/
	int numberOfNonTrivBSCCStates = 0;

    /* The structure, which stores the samples */
    /* NOTE: We only create a sample vector for non-trivial BSCCs */
        PTSampleVecSteady * pSampleVecSteady = (PTSampleVecSteady *) calloc(
                (size_t) numberOfNonTrivBSCCs, sizeof (PTSampleVecSteady));

    /* Allocate an empty sample vector for every reachable non-trivial BSCC containing Psi states */
    /* This array of sample vectors contains the samples corresponding to the non-trivial BSCCs */
    /* The order is important, i.e. it is the same order of BSCCs as in ppNonTrivBSCCBitSets */
    for(i = 0; i < numberOfNonTrivBSCCs; i++) {
        numberOfBSCCStates = count_non_zero( ppNonTrivBSCCBitSets[i] );
        /* Count the overall number of states in non-trivial BSCCs */
        numberOfNonTrivBSCCStates += numberOfBSCCStates;
        /* If pure regenerative method is used, the initial state for the regenerative method */
        /* is chosen "randomly", here: the BSCC state with the smallest state number (index) */
        if(isRegMethodModePure()){
            initial_state_sim = get_idx_next_non_zero( ppNonTrivBSCCBitSets[i], -1 );
            /* Allocate & initialise the structure that holds the sample vectors */
            pSampleVecSteady[i] = allocateSampleVectorSteady(0, initial_state_sim);
        } else {
            if(isRegMethodModeHeuristic()) {
                /* If heuristic regenerative method is used, the initial state is chosen according to */
                /* the heuristic: simulate a sample for every BSCC with length */
                /* sampleBSCCDimensionMultiplier x number of states in the BSCC, taking the most visited */
                /* state within this sample as regeneration state */
                initial_state_sim = get_idx_next_non_zero( ppNonTrivBSCCBitSets[i], -1 );
                pSimStateOccurrence = (int *) calloc((size_t) n_states,
                                sizeof(int));
                pSimStateOccurrence[initial_state_sim]++;
                /* Simulate a sample of length sampleBSCCDimensionMultiplier x number of states in the BSCC */
                /* Store the number of occurrences for every state */
                for(j = 0;j < numberOfBSCCStates * sampleBSCCDimensionMultiplier;j++){
                    initial_state_sim = computeNextState(pStateSpace,
                                                initial_state_sim);
                    pSimStateOccurrence[initial_state_sim]++;
                }
                /* Choose the state as regeneration state, that was most visited during the sample */
                for(j = 0;j < n_states;j++){
                    if(pSimStateOccurrence[j] > pSimStateOccurrence[initial_state_sim]){
                        initial_state_sim = j;
                    }
                }

                /* Allocate & initialise the structure that holds the sample vectors */
                pSampleVecSteady[i] = allocateSampleVectorSteady(0, initial_state_sim);
                /* Free memory */
                free(pSimStateOccurrence);
            }
        }
    }
    /* Return the number of non-trivial BSCCs */
    (*pNumberOfNonTrivBSCCStates) = numberOfNonTrivBSCCStates;
    /* Return the initialised vector of SS samples */
    return pSampleVecSteady;
}
