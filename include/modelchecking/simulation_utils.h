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
*	here we intend to define the utility functions and data structures.
*/

#ifndef SIMULATION_UTILS_H
#define SIMULATION_UTILS_H

#include "sample_vec.h"
#include "sparse.h"

	/**
	* This is the function pointer that gives a template for the function
	* that will do numerical model checking of the unbounded reachability
	* property. In principle we are interested only in the Eventually until.
	* As it is be used in the Hybrid model checking of the steady-state
	* operator of CSL (on CTMCs). Nevertheless, for the optimization reasons,
	* we keep two parameters of the unbounded until.
	* @param bitset the set of allowed states.
	* @param bitset the set of states we want to reach.
	*/
	typedef double * (* TPFunctNumUnbUntilCTMC ) ( const bitset *, const bitset * );

	/**
	* This is the function pointer that gives a template for the function
	* that will do the BSCC search for model checking steady-state operators
	* of CSL in simulation mode.
	* @param sparse the sparse matrix
	* @param bitset the set of good states
	* @param bitset the bitset indicating all states belonging to non-trivial BSCCs
	* @param int the number of non-trivial BSCCs containing good states
	*/
	typedef bitset * (* TPFunctAllowedBSCCSearchCTMC ) ( const sparse *, const bitset *, bitset ***, int * );

	/**
	* This is the function pointer that gives a template for the function
	* that will do the reachability search. In general it is for the exists
	* unbounded until operator (\exists \Phi U \Psi). So we look for the states
	* from which there is a path to \Psi states via \Phi states.
	* @param sparse the sparse matrix
	* @param bitset the set of \Phi states
	* @param bitset the set of \Psi states
	* @return set of states from which we can reach \Psi states via \Phi states
	*/
	typedef bitset * (* TPFunctExistUnbUntil) (const sparse *, const bitset *, const bitset *);

	/**
	* This is the function pointer that gives a template for the function
	* that will do the reachability search. In general it is for the always
	* unbounded until operator (\always \Phi U \Psi). So we look for the states
	* from which all the paths go to \Psi states via \Phi states.
	* @param sparse the sparse matrix
	* @param bitset the set of \Phi states
	* @param bitset the set of \Psi states
	* @param bitset the set of states for which (\exists \Phi U \Psi) holds
	* @return set of states from which we always reach \Psi states via \Phi states
	*/
	typedef bitset * (* TPFunctAlwaysUnbUntil) (const sparse *, const bitset *, const bitset *, const bitset * );

	/****************************************************************************/
	/**************CHECK THE CONF. INT. AGAINST THE PROB. CONSTRAINT*************/
	/****************************************************************************/

	/**
	* This method sorts out the pAlwaysBitSet and pNeverBitSet states and updates the sets
	* pYesBitSet and pNoBitSet with respect to the probability constraint
	* "comparator prob_bound".
	* NOTE: both pAlwaysBitSet and pNeverBitSet are allowed to be NULL
	* @param isSimOneInitState TRUE if we simulate one initial state only
	* @param pYesBitSet this bitset is going to be filled with the states which
	*			satisfy the formula
	* @param pNoBitSet this bitset is going to be filled with the states which
	*			do not satisfy the formula
	* @param pAlwaysBitSet the states that satisfy the formula with probability >= 1.0
	* @param pNeverBitSet the states that satisfy the formula with probability <= 0.0
	* @param initial_state the initial state index if isSimOneInitState == TRUE
	* @param comparator the comparator, one of:
	*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
	*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
	* @param prob_bound the probability bound
	*/
        extern void considerAlwaysAndNeverStates(const BOOL
                                                isSimOneInitState_local,
						bitset * pYesBitSet, bitset * pNoBitSet,
						const bitset * pAlwaysBitSet, const bitset * pNeverBitSet,
						const int initial_state, const int comparator,
						const double prob_bound );

	/**
	* This function updated the pYesBitSet and pNoBitSet bitsets with respect
	* to the model checking result mc_result for the initial_state state.
	* @param mc_result the model checking result in the three valued logic
	* @param initial_state the initial state for which model checking was done
	* @param pYesBitSet the initial_state'th bit of this set will be set on
	*			if mc_result == TVL_TT
	* @param pNoBitSet the initial_state'th bit of this set will be set on
	*			if mc_result == TVL_FF
	*/
	extern inline void markYesNoSetEntree(const TV_LOGIC mc_result, const int initial_state,
					bitset * pYesBitSet, bitset * pNoBitSet);

	/**
	* This method checks if "probability comparator prob_bound" holds.
	* @param probability the left side of the constraint, i.e. the probability
	* @param comparator the comparator, one of:
	*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
	*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
	* @param prob_bound the probability bound
	* @return TRUE if "probability comparator prob_bound" holds, otherwise FALSE
	*/
	extern inline BOOL doesComparisonHold( const double probability, const int comparator, const double prob_bound );

	/**
	* This is an implementation of the checkBoundVSConfInt function from the PhD thesis
	* of Ivan S. Zapreev, see Chapter 6.
	* Plus if "ciRightBorder - ciLeftBorder >= indiff_width" then the method returns TVL_NN
	* @param comparator the comparator, one of:
	*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
	*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
	* @param prob_bound the probability bound
	* @param ciLeftBorder the left border of the conf. int.
	* @param ciRightBorder the right border of the conf. int.
	* @param indiff_width the indifference-region width
	* @return returns one of: TVL_TT, TVL_FF, TVL_NN
	*/
        extern
	TV_LOGIC checkBoundVSConfInt( const int comparator, const double prob_bound, const double ciLeftBorder,
					const double ciRightBorder, const double indiff_width );

	/**
	* This function checks the precondition of the conf. int. like it is done in checkBoundVSConfInt(...)
	* Namely we check if the width of the confidence interval ( ciRightBorder - ciLeftBorder ) for the state
	* initial_state is >= indiff_width. Then the sample size is insufficient and we have an error state which
	* we indicate by setting isIndiffErrStates to TRUE and adding the error state to the pErrStatesBitset bitset.
	*
	* NOTE: Just in case we also check if the conf. int. is invalid i.e. ( ciRightBorder - ciLeftBorder ) < 0
	*
	* @param ciLeftBorder the left border of the conf int.
	* @param ciRightBorder the right border of the conf. int.
	* @param initial_state the initial state for which the conf. int. was computed
	* @param indiff_width the supremum of the allowed conf. int. width
	* @param pIsIndiffErrStates the pointer to the boolean variable that should be TRUE if an error state exists
	* @param pErrStatesBitset the bitset containing the indexes of error states
	*/
        extern
	inline void checkErrorConfIntState(const double ciLeftBorder,
			const double ciRightBorder, const int initial_state,
			const double indiff_width, BOOL * pIsIndiffErrStates,
			bitset * pErrStatesBitset);

	/****************************************************************************/
	/**************THE INCREASE FUNCTION WITH THE STEP AND MAX BOUND*************/
	/****************************************************************************/

	/**
	* This method does the following:
	* 1. (* pTarget) += inc_step;
	* 2. if (* pTarget) > max_val then
	*	(* pTarget) = max_val;
	* @param pTarget the pointer to the variable that will be increased.
	* @param inc_step the increment step
	* @param max_val the maximum value
	*/
	extern inline void increment( int * pTarget, const int inc_step, const int max_val );

	/****************************************************************************/
	/**************THE CONF. INT. COMPUTING FUNCTIONS FOR DIFF SAMPLES***********/
	/****************************************************************************/

	/* The standard conf int as resulted from the central limit theorem */
#       define STANDARD_CONF_INT 0
	/* The improved version of the standard conf. int. for Bernoulli Trials */
#       define AGRESTI_COULL_CONF_INT 1

	/**
	* This method computes different conf intervals based on the provided values
	* Note: this method allows to use two different parameters for the left and
	* the right conf. int. borders. This is needed for computing the intervals
	* in case of the unbounded-reachability problem.
	* @param conf_zeta the zeta value derived from the confidence
	* @param pCiLeftBorder the pointer to the left conf. int. border
	* @param pCiRightBorder the pointer to the right conf. int. border
	* @param left_param the parameter for computing the left conf. int. border
	* @param right_param the parameter for computing the right conf. int. border
	* @param theConfIntType indicate the kind of the conf. int. to be computed:
	*				STANDARD_CONF_INT, AGRESTI_COULL_CONF_INT, etc.
	*/
	extern inline void computeBordersConfInt(const double conf_zeta, double * pCiLeftBorder,
                                        double *pCiRightBorder, int sample_size,
                                        int left_param, int right_param,
					const int theConfIntType);

	/**
        * This method computes the conf. int. borders for the given sample, i.e.
	* we compute and set: (*pCiLeftBorder), (*pCiRightBorder).
	* NOTE: The borders that do not fit into the interval R[0,1] are truncated.
	* @param conf_zeta the zeta value derived from the confidence
	* @param pCiLeftBorder the pointer to the left conf. int. border
	* @param pCiRightBorder the pointer to the right conf. int. border
	* @param pSampleVecUntil the sample we use to derive the conf int
	* @param theConfIntType indicate the kind of the conf. int. to be computed:
	*				STANDARD_CONF_INT, AGRESTI_COULL_CONF_INT, etc.
	*/
	extern inline void computeBordersUU( const double conf_zeta, double * pCiLeftBorder,
					double * pCiRightBorder, const PTSampleVecUntil pSampleVecUntil,
					const int theConfIntType );

	/**
	* This method computes the conf. int. borders of steady-state simulation for
	* the given sample, i.e. we compute and set: (*pCiLeftBorder), (*pCiRightBorder).
	* For more information see PhD Thesis of Ivan S. Zapreev.
	* NOTE: The borders that do not fit into the interval R[0,1] are truncated.
	* @param reg_cycles the number of samples, i.e. regeneration cycles calculated so far
	* @param conf_zeta the zeta value derived from the confidence
	* @param pCiLeftBorder the pointer to the left conf. int. border
	* @param pCiRightBorder the pointer to the right conf. int. border
	* @param pSampleVecSteady the sample we use to derive the conf int
	*/
	extern inline void computeBordersConfIntSS( const int reg_cycles, const double conf_zeta, double * pCiLeftBorder,
		double * pCiRightBorder, const PTSampleVecSteady pSampleVecSteady );

	/****************************************************************************/
	/**********THE FUNCTION THAT SIMULATES THE SAMPLE FOR UNBOUNDED UNTIL********/
	/****************************************************************************/

	/**
	* This function is an implementation of the extendSamples function from the PhD Thesis
	* of Ivan S. Zapreev. It either add new observations to the sample and simulates them
	* or it increases the simulation depth of every observation in the sample.
	* NOTE: This method can increase the sample size and the simulation depth at the same time.
	* WARNING: We simulate the states ASSUMING there are no self loops!
	* @param pStateSpace the sparse matrix of the embedded DTMC with good ad bad
	*			states made absorbing
	* @param pSampleVecUntil the sample vector obtained on the previous iteration
	* @param new_sample_size the new sample size
	* @param new_simulation_depth the new simulation depth
	* @param pGoodStates the bitsets containing all the good absorbing states
	* @param pTransientStates the bitsets containing all the transient states
	*/
	extern void simulateSampleVectorUnbUntilDTMC( const sparse * pStateSpace, PTSampleVecUntil pSampleVecUntil,
							const int new_sample_size, const int new_simulation_depth,
							const bitset * pGoodStates, const bitset * pTransientStates );

	/**
        * Computes the reachability probabilities for the trivial and
        * non-trivial
	* BSCCs. In this step, we treat all states as initial because of the global
	* numerical-model checking procedure.
	* I.e. we compute the probabilities P(true U BSCC_{i})
	* @param n_states the number of states in the model
	* @return the array of the size equal to the number of BSCCs. Each element of the
	* @param numberOfBSCCs the total number of considered BSCCs
	* @param numberOfTrivialBSCCs the number of trivial BSCCs from the total number of
	* @param ppNonTrivBSCCBitSets the array storing pointers to bitsets each of which
	* contains states of the corresponding non-trivial BSCCs
	* @param pTrivialBSCCBitSet indicates states each of which is a trivial BSCC
	* BSCCs given by numberOfBSCCs
	* @param pFNumUnbUntilCTMC the function for computing the unbounded-until
	*							probabilities numerically (Hybrid mode).
	* @return an array of arrays each index in the enclosing array is the index of the BSCC
	* the data in the related sub array is the reachability probabilities for this BSCC from
	* the initial states given in pInitialStatesBitSet
	*/
	extern double ** getBSCCReachabilityProbsNumerical( const int n_states, const int numberOfBSCCs,
														const int numberOfTrivialBSCCs,
														bitset ** const ppNonTrivBSCCBitSets,
														const bitset * pTrivialBSCCBitSet,
														const TPFunctNumUnbUntilCTMC pFNumUnbUntilCTMC );

	/**
	 * This function allows to find trivial and non-trivial BSCCs that contain \Psi states
	 * @param pStateSpace the state space
	 * @param pPsiBitSet the \Psi states
	 * @param ppNonTrivBSCCStates the pointer to an array that by the end of the function
	 *								will store the array of non-trivial BSCC states
         * @param ppTrivialBSCCBitSet the pointer to a bitset pointer, it will
         *                      be set to a pointer to a bitset containing
         *                      states of the trivial BSCCs
	 * @param pNumberOfNonTrivBSCCs this is the pointer to a variable that will store the number
	 *								of non-trivial BSCCs
	 * @param pNumberOfTrivialBSCCs this is the pointer to a variable that will store the number
	 *								of trivial BSCCs
	 * @param pNumberOfBSCCs this is the pointer to a variable that will store the total number
	 *								of BSCCs
	 * @param pFAllowedBSCCSearchCTMC the pointer to a function that allows to search for BSCCs
         * @return the array of bitset pointers, each of which contains the
         *                      states of the
	 *			corresponding non-trivial BSCC.
	 */
	extern bitset ** getTrivialAndNonTrivialPsiBSCC( const sparse * pStateSpace, const bitset * pPsiBitSet,
													int ** ppNonTrivBSCCStates, bitset ** ppTrivialBSCCBitSet,
													int * pNumberOfNonTrivBSCCs, int *pNumberOfTrivialBSCCs,
													int * pNumberOfBSCCs, const TPFunctAllowedBSCCSearchCTMC pFAllowedBSCCSearchCTMC );

	/**
        * This method allows to initialize ppBSCCExistsReachableStates and
        * ppBSCCAlwaysReachableStates.
	* In essence for every BSCC it returns the set of states from which it is sometimes and always reachable
	* @param pStateSpace the state space
	* @param numberOfBSCCs the total number of BSCCs with \Psi states
	* @param pTrivialBSCCBitSet the list of states which are trivial BSCCs (single \Psi states)
	* @param ppNonTrivBSCCBitSets the array of non-trivial BSCCS with \Psi states
	* @param numberOfNonTrivBSCCs the number of non-trivial BSCCs with \Psi states
	* @param ppBSCCExistsReachableStates the return parameter, where for each BSCC
	* we will store the bitset of states from which can reach the BSCC.
	* @param ppBSCCAlwaysReachableStates the return parameter, where for each BSCC
	* we will store the bitset of states from which always reach the BSCC.
	* @param pExistsUnbUntilCTMC the function pointer that gives a template for the function
	*                            that will do the reachability search. In general it is for the exists
	*                            unbounded until operator (\exists \Phi U \Psi). So we look for the states
	*                            from which there is a path to \Psi states via \Phi states.
	* @param pAlwaysUnbUntilCTMC the function pointer that gives a template for the function
	*                            that will do the reachability search. In general it is for the always
	*                            unbounded until operator (\always \Phi U \Psi). So we look for the states
	*                            from which all the paths go to \Psi states via \Phi states.
	*/
	extern void getReachabilityTrivialAndNonTrivialPsiBSCC( const sparse * pStateSpace, const int numberOfBSCCs,
															bitset * const pTrivialBSCCBitSet,
															bitset ** const ppNonTrivBSCCBitSets,
															const int numberOfNonTrivBSCCs,
															bitset ** const ppBSCCExistsReachableStates,
															bitset ** const ppBSCCAlwaysReachableStates,
															const TPFunctExistUnbUntil pExistsUnbUntilCTMC,
															const TPFunctAlwaysUnbUntil pAlwaysUnbUntilCTMC );

	/**
	 * This method takes the initial state and checks which BSCCs are reachable with which probabilities.
	 * This information is stored in the implicit return arguments of the method.
	 * @param initial_state the initial state index to check reachability probabilities for
	 * @param numberOfTrivialBSCCs the total number of trivial BSCCs
	 * @param numberOfNonTrivBSCCs the total number of non-trivial BSCCs
	 * @param ppReachProbability the reachability probabilities from all initial states to BSCCs with Psi states
	 *                             e.g. ppReachProbability[0][0] contains the probability to reach BSCC 0 from state 1
	 *                             ppReachProbability[0][1] contains the probability to reach BSCC 0 from state 2 etc.
	 * @param ppBSCCReachProbability the pointer to the array that will be allocated and filled with all of the
	 *                               reachability probabilities for trivial and non-trivial BSCCs
	 * @param ppReachableTrivialBSCCBitSet the pointer to the pointer for a bitset containing the reachable trivial BSCCs
	 *                                     The bit set is allocated and filled in here.
	 * @param ppReachableNonTrivBSCCBitSet the pointer to the pointer for a bitset containing the reachable non-trivial BSCCs
	 *                                     The bit set is allocated and filled in here.
	 * @param pNumberOfReachableTrivialBSCCs the pointer to the value that by the end of the method call will store the
	 *                                       number of reachable (trivial) BSCCs
	 * @param pNumberOfReachableNonTrivBSCCs the pointer to the value that by the end of the method call will store the
	 *                                       number of reachable (non-trivial) BSCCs
	 * @return the total number of reachable BSCCs
	 */
	extern int extractReachabilityProbs( const int initial_state, const int numberOfTrivialBSCCs, const int numberOfNonTrivBSCCs,
											double ** ppReachProbability, double ** ppBSCCReachProbability,
											bitset ** ppReachableTrivialBSCCBitSet, bitset ** ppReachableNonTrivBSCCBitSet,
											int * pNumberOfReachableTrivialBSCCs, int *pNumberOfReachableNonTrivBSCCs );

	/**
	 * This method extracts the reachable BSCCs for the given initial state.
	 * @param initial_state the initial state index
	 * @param numberOfTrivialBSCCs the number of trivial BSCCs
	 * @param numberOfNonTrivBSCCs the number of non-trivial BSCCs
	 * @param ppBSCCExistsReachableStates the array of bitset pointers (with index i), each ppBSCCExistsReachableStates[i]
	 *                                    is the set of states from which BSCC_{i} can be reached
	 * @param ppBSCCAlwaysReachableStates the array of bitset pointers (with index i), each ppBSCCExistsReachableStates[i]
	 *                                    is the set of states from which BSCC_{i} is always reached
	 * @param ppExistsReachableTrivialBSCCBitSet the pointer to a bitset pointer what will contain all the trivial BSCCs
	 *                                           indices that can be reached from the given initial states
	 * @param ppExistsReachableNonTrivBSCCBitSet the pointer to a bitset pointer what will contain all the non-trivial
	 *                                           BSCC indices that can be reached from the given initial states
	 * @param pNumberOfExistsReachableTrivialBSCCs the number of BSCCs in ppExistsReachableTrivialBSCCBitSet
	 * @param pNumberOfExistsReachableNonTrivBSCCs the number of BSCCs in ppExistsReachableNonTrivBSCCBitSet
	 * @param the total number of reachable BSCCs
	 * @return TRUE if there is just one BSCC that is always reached from the actual initial state.
	 */
        extern
	BOOL extractReachabilityData( const int initial_state, const int numberOfTrivialBSCCs, const int numberOfNonTrivBSCCs,
                        const bitset * const *const ppBSCCExistsReachableStates,
                        const bitset * const *const ppBSCCAlwaysReachableStates,
								bitset ** ppReachableTrivialBSCCBitSet, bitset ** ppReachableNonTrivBSCCBitSet,
								int * pNumberOfExistsReachableTrivialBSCCs, int * pNumberOfExistsReachableNonTrivBSCCs );

	/**
	 * This function allows to search for all transient states in the model
         * TODO: May be its current performance is not optimal.
	 * @param pStateSpace the state space
	 * @param pFAllowedBSCCSearchCTMC the function that allows to search for BSCCs
	 * @return the pointer to the bitset (size of state space) that contains all transient states
	 */
	extern bitset * getTransientStates( const sparse * pStateSpace, const TPFunctAllowedBSCCSearchCTMC pFAllowedBSCCSearchCTMC );

#endif
