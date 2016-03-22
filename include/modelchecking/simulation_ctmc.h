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
*	here we intend to define the CTMC-related functions and data
*	structures.
*/

#ifndef SIMULATION_CTMC_H
#define SIMULATION_CTMC_H

#include "simulation_utils.h"

	/****************************************************************************/
	/**********************THE FUNCTIONS FOR SIMULATING CTMCs********************/
	/****************************************************************************/

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
        extern void modelCheckUnboundedUntilCTMC(sparse * pStateSpace,
                                                        const double *
                                                        pCTMCRowSums,
							const double confidence, const bitset *pPhiBitSet,
							const bitset * pPsiBitSet, bitset ** ppYesBitsetResult,
							bitset ** ppNoBitsetResult, double ** ppProbCILeftBorder,
							double ** ppProbCIRightBorder, int * pResultSize,
							const int comparator, const double prob_bound,
                                                        const int initial_state,
                                                        const BOOL
                                                        isSimOneInitState_local,
							unsigned int * pMaxNumUsedObserv );

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
        extern void modelCheckTimeIntervalUntilCTMC(sparse * pStateSpace,
                                                        const double *
                                                        pCTMCRowSums,
							const double confidence, const bitset *pPhiBitSet,
							const bitset * pPsiBitSet, const double left_time_bound,
							const double right_time_bound, bitset ** ppYesBitsetResult,
							bitset ** ppNoBitsetResult, double ** ppProbCILeftBorder,
							double ** ppProbCIRightBorder, int * pResultSize,
							const int comparator, const double prob_bound,
                                                        const int initial_state,
                                                        const BOOL
                                                        isSimOneInitState_local,
							unsigned int * pMaxNumUsedObserv );

	/**
	* This function is an implementation of the steadyStateHybrid function from
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
	extern void modelCheckSteadyStateHybridCTMC( sparse * pStateSpace, const double * pCTMCRowSums ,
                                                        const double confidence,
                                                        const bitset*pPsiBitSet,
							bitset ** ppYesBitsetResult, bitset ** ppNoBitsetResult,
							double ** ppProbCILeftBorder, double ** ppProbCIRightBorder,
							int * pResultSize, const int comparator, const double prob_bound,
                                                        const int initial_state,
                                                        const BOOL
                                                        isSimOneInitState_local,
							const TPFunctNumUnbUntilCTMC pFNumUnbUntilCTMC,
							const TPFunctAllowedBSCCSearchCTMC pFBSCCSearchCTMC,
							const double error_bound, unsigned int * pMaxNumUsedObserv );


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
        extern void modelCheckSteadyStatePureCTMC(sparse * pStateSpace,
                                                        const double *
                                                        pCTMCRowSums,
							const double confidence, const bitset * pPsiBitSet,
							bitset ** ppYesBitsetResult, bitset ** ppNoBitsetResult,
							double ** ppProbCILeftBorder, double ** ppProbCIRightBorder,
							int * pResultSize, const int comparator, const double prob_bound,
                                                        const int initial_state,
                                                        const BOOL
                                                        isSimOneInitState_local,
							const TPFunctExistUnbUntil pExistsUnbUntilCTMC,
							const TPFunctAlwaysUnbUntil pAlwaysUnbUntilCTMC,
							const TPFunctAllowedBSCCSearchCTMC pFAllowedBSCCSearchCTMC,
							unsigned int * pMaxNumUsedObserv );

#endif
