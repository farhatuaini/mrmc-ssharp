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
*	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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
*	here we intend to define the common functions and data structures.
*/

#ifndef SIMULATION_COMMON_H
#define SIMULATION_COMMON_H

#include "sample_vec.h"
#include "sparse.h"
#include "rand_num_generator.h"

#include <stdarg.h>

	/****************************************************************************/
	/********THE COMMON FUNCTION FOR INVOCATION THE SIMULATION PROCEDURE ********/
	/****************************************************************************/

	/**
	* This is a template method for all the method that model check one single initial state
	* @param pStateSpace the state-space matrix
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
	* @param pNumUsedObserv the return-value variable for storing the number of used observations
	* @param args the extra arguments needed for a specific model-checking procedure
	*/

	typedef void ( * TPFunctMCOneState ) ( const sparse* pStateSpace, const int initial_state, const double indiff_width,
						const double confidence, const int comparator, const double prob_bound,
						const int arr_index, double *pCiLeftBorders, double *pCiRightBorders,
						bitset * pYesBitSet, bitset * pNoBitSet, unsigned int *pNumUsedObserv, va_list args );

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
        extern
	void modelCheckStatesCommon( const sparse * pStateSpace, const double confidence, bitset ** ppYesBitSet,
					bitset ** ppNoBitSet, double ** ppCiLeftBorders, double ** ppCiRightBorders,
					int * pResultSize, const int comparator, const double prob_bound,
					const int initial_state, const BOOL isSimOneInitState,
					const bitset * pValidInitialStates, unsigned int *pMaxNumUsedObserv,
					TPFunctMCOneState pFMCOneState, ... );

	/****************************************************************************/
	/*****THE COMMON FUNCTION FOR SIMULATING CSL/PCTL operators on CTMC/DTMC*****/
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
        extern
	void modelCheckOneStateUUCommon( const sparse* pStateSpace, const int initial_state, const double indiff_width,
						const double confidence, const int comparator, const double prob_bound,
						const int arr_index, double *pCiLeftBorders, double *pCiRightBorders,
						bitset * pYesBitSet, bitset * pNoBitSet, unsigned int *pNumUsedObserv,
						va_list pTransAndGoodStatesArguments );

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
        extern
	PTSampleVecSteady * prepareSSSimulationSample(const sparse * pStateSpace, const int numberOfNonTrivBSCCs,
													bitset ** const ppNonTrivBSCCBitSets, int * pNumberOfNonTrivBSCCStates );

	/**
	 * This function computes to which state we will go.
         * @param pM                the sparse matrix describing the state space
         * @param current_obs_state the current state index
	 * WARNING: The self loops are not taken into account, we expect there to be NONE
	 * @return the next state
	 */
        extern state_index computeNextState(
                        /*@observer@*/ /*@sef@*/ const sparse * pM,
                        /*@sef@*/ state_index current_obs_state);

#       define computeNextState(pM,current_obs_state) \
                /* If we are not in an absorbing state */ \
                (0 != mtx_next_num((pM), (current_obs_state)) \
                        ? /* Compute to what state we will go */ \
                          generateRandNumberDiscrete( \
                                (pM)->valstruc[(current_obs_state)].col, \
                                (pM)->valstruc[(current_obs_state)].val, \
                                mtx_next_num((pM), (current_obs_state))) \
                        : (current_obs_state))

	/**
         * This function simulates the unbounded reachability samples for model
         * checking the
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
	 * @param pTrivialBSCCBitSet the bitset that contains trivial \PSI BSCCs (single-state BSCCs)
	 * @param ppNonTrivBSCCBitSets the non-trivial BSCC's (their states) with \Psi states in them
	 * @param isOdd indicates whether it is an odd iteration or not
         * @param pReachableTrivialBSCCBitSet contains the indices of reachable
         *                                    trivial \Psi BSCCs
         * @param pReachableNonTrivBSCCBitSet contains the indices of reachable
         *                                    non-trivial BSCCs with \Psi states
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
        extern
	BOOL simulateAndComputeAllReachabilityConfInt( const sparse * const pStateSpace, const bitset * const pTransientBitSet,
                        const bitset * const pTrivialBSCCBitSet,
                        const bitset * const * const ppNonTrivBSCCBitSets,
													const BOOL isOdd, const bitset * const pReachableTrivialBSCCBitSet,
													const bitset * const pReachableNonTrivBSCCBitSet, PTSampleVecUntil * pSampleVecUntilOneArray,
													PTSampleVecUntil * pSampleVecUntilTwoArray, double * pLeftBorderUU, double * pRightBorderUU,
													int * pSampleSizeUU, const int max_sample_size_uu, const int sample_size_step_uu,
													int * pSimDepthUU, const int max_sim_depth_uu, const int sim_depth_step_uu, const double zeta );

#endif
