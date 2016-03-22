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
*	here we intend to define the utility functions and data structures.
*/

#include "simulation_utils.h"

#include "parser_to_tree.h"
#include "simulation_common.h"

#include <math.h>

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
inline void considerAlwaysAndNeverStates( const BOOL isSimOneInitState,
					bitset * pYesBitSet, bitset * pNoBitSet,
					const bitset * pAlwaysBitSet, const bitset * pNeverBitSet,
					const int initial_state, const int comparator,
					const double prob_bound ){
	/* Check if the probability constraint is satisfied for probabilities 1.0 or 0.0 */
	BOOL isAlwaysSetGood = doesComparisonHold( 1.0, comparator, prob_bound );
	BOOL isNeverSetGood =  doesComparisonHold( 0.0, comparator, prob_bound );

	if( isSimOneInitState ){
		/* If we simulate just one initial state than all the other data is not iteresting */
		if( ( pAlwaysBitSet != NULL ) && get_bit_val( pAlwaysBitSet, initial_state ) ){
			if( isAlwaysSetGood ){
				set_bit_val( pYesBitSet, initial_state, BIT_ON );
			} else {
				set_bit_val( pNoBitSet, initial_state, BIT_ON );
			}
		} else {
			if( ( pNeverBitSet != NULL ) && get_bit_val( pNeverBitSet, initial_state ) ){
				if( isNeverSetGood ){
					set_bit_val( pYesBitSet, initial_state, BIT_ON );
				} else {
					set_bit_val( pNoBitSet, initial_state, BIT_ON );
				}
			}
		}
	} else {
		/* All pAlwaysBitSet states satisfy the formula with probability >= 1.0 */
		if( pAlwaysBitSet != NULL ){
			if( isAlwaysSetGood ){
				or_result( pAlwaysBitSet, pYesBitSet );
			} else {
				or_result( pAlwaysBitSet, pNoBitSet );
			}
		}
		/* All pNeverBitSet states satisfy the formula with probability <= 0.0 */
		if( pNeverBitSet != NULL ){
			if( isNeverSetGood ){
				or_result( pNeverBitSet, pYesBitSet );
			} else {
				or_result( pNeverBitSet, pNoBitSet );
			}
		}
	}
}

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
inline void markYesNoSetEntree(const TV_LOGIC mc_result, const int initial_state,
				bitset * pYesBitSet, bitset * pNoBitSet){
	switch( mc_result ){
		case TVL_TT:
			set_bit_val( pYesBitSet, initial_state, BIT_ON);
			break;
		case TVL_FF:
			set_bit_val( pNoBitSet, initial_state, BIT_ON);
			break;
                case TVL_NN:
                        break;
                default:
                        fprintf(stderr, "markYesNoSetEntree: "
                                        "illegal parameter\n");
                        exit(EXIT_FAILURE);
	}
}

/**
* This method checks if "probability comparator prob_bound" holds.
* @param probability the left side of the constraint, i.e. the probability
* @param comparator the comparator, one of:
*	COMPARATOR_SF_GREATER, COMPARATOR_SF_GREATER_OR_EQUAL,
*	COMPARATOR_SF_LESS, COMPARATOR_SF_LESS_OR_EQUAL
* @param prob_bound the probability bound
* @return TRUE if "probability comparator prob_bound" holds, otherwise FALSE
*/
inline BOOL doesComparisonHold( const double probability, const int comparator, const double prob_bound ){
	switch( comparator ){
		case COMPARATOR_SF_GREATER:
			return ( probability > prob_bound ? TRUE : FALSE );
		case COMPARATOR_SF_GREATER_OR_EQUAL:
			return ( probability >= prob_bound ? TRUE : FALSE );
		case COMPARATOR_SF_LESS:
			return ( probability < prob_bound ? TRUE : FALSE );
		case COMPARATOR_SF_LESS_OR_EQUAL:
			return ( probability <= prob_bound ? TRUE : FALSE );
		default:
			printf("ERROR: An unknown comparator type %d.\n", comparator);
                        exit(EXIT_FAILURE);
	}
}

static inline TV_LOGIC checkBoundVSConfIntHelper(BOOL theTRUECond, BOOL theFALSECond){
	TV_LOGIC tvl_result = TVL_NN;
	if( theTRUECond ){
		tvl_result = TVL_TT;
	} else {
		if( theFALSECond ){
			tvl_result = TVL_FF;
		}
	}
	return tvl_result;
}

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
TV_LOGIC checkBoundVSConfInt( const int comparator, const double prob_bound, const double ciLeftBorder,
				const double ciRightBorder, const double indiff_width ){
	TV_LOGIC tvl_result = TVL_NN;

	/* First check if it is a valid conf. int. */
	if( ( ciRightBorder - ciLeftBorder ) < indiff_width ){
		/* Distinguish the cases based on the comparator type */
		switch( comparator ){
			case COMPARATOR_SF_LESS_OR_EQUAL:
				tvl_result = checkBoundVSConfIntHelper( prob_bound >= ciRightBorder,
									prob_bound < ciLeftBorder );
				break;
			case COMPARATOR_SF_LESS:
				tvl_result = checkBoundVSConfIntHelper( prob_bound > ciRightBorder,
									prob_bound <= ciLeftBorder );
				break;
			case COMPARATOR_SF_GREATER_OR_EQUAL:
				tvl_result = checkBoundVSConfIntHelper( prob_bound <= ciLeftBorder,
									prob_bound > ciRightBorder );
				break;
			case COMPARATOR_SF_GREATER:
				tvl_result = checkBoundVSConfIntHelper( prob_bound < ciLeftBorder,
									prob_bound >= ciRightBorder );
				break;
			default:
				printf("ERROR: An unexpected comparator type: %d.\n", comparator);
                                exit(EXIT_FAILURE);
		}
	}

	return tvl_result;
}

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
inline void checkErrorConfIntState(const double ciLeftBorder,
		const double ciRightBorder, const int initial_state,
		const double indiff_width, BOOL * pIsIndiffErrStates,
		bitset * pErrStatesBitset)
{
	if( ( ciRightBorder - ciLeftBorder ) < 0 || ( ( ciRightBorder - ciLeftBorder ) >= indiff_width ) ){
		*pIsIndiffErrStates = TRUE;
		set_bit_val(pErrStatesBitset, initial_state, BIT_ON);
	}
}

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
inline void increment( int * pTarget, const int inc_step, const int max_val ){
	IF_SAFETY( pTarget != NULL )
		(* pTarget) += inc_step;
		if( (* pTarget) > max_val ){
			(* pTarget) = max_val;
		}
	ELSE_SAFETY
		printf("ERROR: Trying to increment a value defined by a NULL pointer.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/****************************************************************************/
/**************THE INCREASE FUNCTION WITH THE STEP AND MAX BOUND*************/
/****************************************************************************/

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
inline void computeBordersConfInt(const double conf_zeta, double * pCiLeftBorder,
                                        double *pCiRightBorder, int sample_size,
                                        int left_param, int right_param,
					const int theConfIntType){
	double left_mean, right_mean, left_variance_sqrt, right_variance_sqrt;
	double work_sample_size = sample_size; /* For some conf. int.the sample size is changed */
	const double conf_zeta_squared = conf_zeta * conf_zeta;

	/* Compute the conf. int. borders */
	switch( theConfIntType ){
		case STANDARD_CONF_INT:
			left_mean = left_param / work_sample_size;
			right_mean = right_param / work_sample_size;
			left_variance_sqrt = sqrt( left_mean * ( sample_size - left_param ) / ( sample_size - 1.0 ) );
			right_variance_sqrt = sqrt( right_mean * ( sample_size - right_param ) / ( sample_size - 1.0 ) );
			break;
		case AGRESTI_COULL_CONF_INT:
			work_sample_size += conf_zeta_squared;
			left_mean = ( left_param + 0.5 * conf_zeta_squared ) / work_sample_size;
			right_mean = ( right_param + 0.5 * conf_zeta_squared ) / work_sample_size;
			left_variance_sqrt = sqrt( left_mean * ( 1.0 - left_mean ) );
			right_variance_sqrt = sqrt( right_mean * ( 1.0 - right_mean ) );
			break;
		default :
			printf("ERROR: An unexpected conf. int. type: %d.\n", theConfIntType);
                        exit(EXIT_FAILURE);
	}

	/* Compute the left border */
	( * pCiLeftBorder ) = left_mean - conf_zeta * left_variance_sqrt / sqrt(work_sample_size);
	/* Compute the right border */
	( * pCiRightBorder ) = right_mean + conf_zeta * right_variance_sqrt / sqrt(work_sample_size);

	/* Check if we've gotten values < 0.0 and > 1.0 */
	if( ( * pCiLeftBorder ) < 0.0 ){
		( * pCiLeftBorder ) = 0.0;
	}
	if( ( * pCiRightBorder ) > 1.0 ){
		( * pCiRightBorder ) = 1.0;
	}
}

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
inline void computeBordersUU( const double conf_zeta, double * pCiLeftBorder, double * pCiRightBorder,
				const PTSampleVecUntil pSampleVecUntil, const int theConfIntType){
	IF_SAFETY( pSampleVecUntil != NULL )
		IF_SAFETY( ( pCiLeftBorder != NULL ) && ( pCiRightBorder != NULL ) )
			/* Compute the conf. int. for Bernoulli trials */
			computeBordersConfInt( conf_zeta, pCiLeftBorder, pCiRightBorder,
                                pSampleVecUntil->sampleVecBase.curr_sample_size,
                                pSampleVecUntil->sum_good,
                                pSampleVecUntil->sum_good
                                                + pSampleVecUntil->sum_trans,
									theConfIntType );
		ELSE_SAFETY
			printf("ERROR: The left/right conf. int. border is passed by the NULL pointer.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: The sample vector is passed by the NULL pointer.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

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
inline void computeBordersConfIntSS( const int reg_cycles, const double conf_zeta, double * pCiLeftBorder,
					double * pCiRightBorder, const PTSampleVecSteady pSampleVecSteady ){
	int i;
	double point_estimate, half_width_ci;
	double good_mean = 0; double time_mean = 0;
	double point_estimate_variance;
	double good_variance = 0; double time_variance = 0; double time_good_covariance = 0;
	double good_helper, time_helper;

	IF_SAFETY( pSampleVecSteady != NULL )
		IF_SAFETY( ( pCiLeftBorder != NULL ) && ( pCiRightBorder != NULL ) )

			/* Compute the mean value for "good" states and time */
			for( i = 0; i < reg_cycles; i++){
				good_mean += pSampleVecSteady->pGoodObservationsVec[i];
				time_mean += pSampleVecSteady->pTimeObservationsVec[i];
			}
			good_mean = good_mean / reg_cycles;
			time_mean = time_mean / reg_cycles;

			/* Compute the variances */
			for( i = 0; i < reg_cycles; i++){
				good_helper = pSampleVecSteady->pGoodObservationsVec[i] - good_mean;
				time_helper = pSampleVecSteady->pTimeObservationsVec[i] - time_mean;

                                /* Multiplying a double with itself is faster
                                   than calling pow(...,2.0). David N. Jansen.*/
                                good_variance += good_helper * good_helper;
                                time_variance += time_helper * time_helper;
				time_good_covariance += good_helper * time_helper;
			}
			good_variance = good_variance / (reg_cycles - 1);
			time_variance = time_variance / (reg_cycles - 1);
			time_good_covariance = time_good_covariance / (reg_cycles - 1);

			/* Compute the point estimate */
			point_estimate = good_mean / time_mean;
			/* Compute the overall variance */
                        point_estimate_variance = good_variance
                                + point_estimate * (-2 * time_good_covariance
                                                + point_estimate*time_variance);

			/* Compute the half width of the conf. int. */
                        half_width_ci = conf_zeta
                                        * sqrt((double) point_estimate_variance)
                                        / (time_mean*sqrt((double) reg_cycles));
			/* Compute the left border */
			( * pCiLeftBorder ) = point_estimate - half_width_ci;
			/* Compute the right border */
			( * pCiRightBorder ) = point_estimate + half_width_ci;

			/* Check if we've gotten values < 0.0 and > 1.0 */
			if( ( * pCiLeftBorder ) < 0.0 ){
				( * pCiLeftBorder ) = 0.0;
			}
			if( ( * pCiRightBorder ) > 1.0 ){
				( * pCiRightBorder ) = 1.0;
			}
		ELSE_SAFETY
			printf("ERROR: The left/right conf. int. border is passed by the NULL pointer.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: The sample vector is passed by the NULL pointer.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/****************************************************************************/
/**********THE FUNCTION THAT SIMULATES THE SAMPLE FOR UNBOUNDED UNTIL********/
/****************************************************************************/

/* These are the same two values and they indicate that the simulation of */
/* new observations should start from the first one and also that all the */
/* observation in the sample were processed */
#define END_OF_OBSERVATIONS -1
#define FOR_ALL_OBSERVATIONS END_OF_OBSERVATIONS

/**
* This function is supposed to perform "depth_steps_needed" successive simulations
* for every observation after and including the "initial_obs_index"
* (if initial_obs_index != -1). The fields of "pSampleVecUntil", such as
*	sum_good, sum_trans, curr_simulation_depth, pTransStateInd, pObservationsVec.
* are updated correspondingly.
* WARNING: We simulate the states ASSUMING there are no self loops!
* @param pStateSpace the sparse matrix of the embedded DTMC with good ad bad
*			states made absorbing
* @param pSampleVecUntil the sample vector obtained on the previous iteration
* @param initial_obs_index the index of the observation we start simulating from
*				in the array pSampleVecUntil->pObservationsVec
* WARNING: For simulating all observation should be set to -1
* @param depth_steps_needed the simulation depth, i.e. the number of times we
*				simulate each observation in the array pObservationsVec
* @param pGoodStates the bitsets containing all the good absorbing states
* @param pTransientStates the bitsets containing all the transient states
*/
static void simulateUnbUntilSampleDTMC( const sparse * pStateSpace, PTSampleVecUntil pSampleVecUntil,
					const int initial_obs_index, const int depth_steps_needed,
					const bitset * pGoodStates, const bitset *pTransientStates ){
        /* Note that the second parameter to get_idx_next_non_zero() is the
           index of the previously found bitset element. */
	/* Therefore if we are not searching from the beginning of the bitset we have to use  */
	/* "initial_obs_index - 1" as the initial index */
	int curr_obs_idx = ( initial_obs_index == -1 ? initial_obs_index : initial_obs_index - 1 );
	bitset * pTransStateInd = pSampleVecUntil->pTransStateInd;
	/* Will store the number of newly visited states */
	unsigned int newlyVisitedStates = 0;

	/* In case we are not adding new observations but just simulating all further */
	if( initial_obs_index == FOR_ALL_OBSERVATIONS ){
		/* Update the current simulation depth */
		pSampleVecUntil->curr_simulation_depth += depth_steps_needed;
	}

	/* For all the transient-state observation in the sample */
	while( ( curr_obs_idx = get_idx_next_non_zero( pTransStateInd, curr_obs_idx  )) != END_OF_OBSERVATIONS ){
		/* Take the current value of the curr_obs_idx's observation */
		int current_obs_state = pSampleVecUntil->pObservationsVec[curr_obs_idx];
		int extra_simulation_depth = 0;

		/* Iterate until we are in the absorbing state or we are on the right depth */
		while( ( extra_simulation_depth < depth_steps_needed ) && get_bit_val( pTransientStates, current_obs_state) ){
			/* WARNING: We simulate the state ASSUMING there are no self loops! */
			/* WARNING: We simulate a pTransientStates state, i.e. there are outgoing */
			/* transitions in this state that are not self-loop transitions */
                        current_obs_state = computeNextState(pStateSpace,
                                        current_obs_state);

			/* The next state was simulated */
			newlyVisitedStates++;

			/* Update the loop-condition variable */
			extra_simulation_depth += 1;
		}

		/* Update the observation pObservationsVec[curr_obs_idx] */
		pSampleVecUntil->pObservationsVec[curr_obs_idx] = current_obs_state;

		/* Update the sum_good, sum_trans and pTransStateInd fields */
		if( ! get_bit_val( pTransientStates, current_obs_state) ){
			pSampleVecUntil->sum_trans -= 1;
                        set_bit_val(pTransStateInd, curr_obs_idx, BIT_OFF);
			if( get_bit_val( pGoodStates, current_obs_state) ){
				pSampleVecUntil->sum_good += 1;
			}
		}
	}

	/* printSampleVectorUntil(pSampleVecUntil); */

	/* Update the number of states visited while preparing this sample */
	( (PTSampleVec) pSampleVecUntil)->num_visited_states += newlyVisitedStates;
}

/**
* This function is an implementation of the "extendSamples" function from the PhD Thesis
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
void simulateSampleVectorUnbUntilDTMC( const sparse * pStateSpace, PTSampleVecUntil pSampleVecUntil,
										const int new_sample_size, const int new_simulation_depth,
										const bitset * pGoodStates, const bitset * pTransientStates ){
	IF_SAFETY( pSampleVecUntil != NULL )
		/* Cast to the parent structure which contains initial_state and curr_sample_size */
		PTSampleVec pSampleVecUntilBase = (PTSampleVec) pSampleVecUntil;

		IF_SAFETY( ( pGoodStates != NULL ) && ( pTransientStates != NULL ) )
			IF_SAFETY( pStateSpace != NULL )
				int depth_steps_needed;

				/* If we need to extend the sample size */
				if( pSampleVecUntilBase->curr_sample_size < new_sample_size ){
					const int old_sample_size = pSampleVecUntilBase->curr_sample_size;

					/* Extend the sample size */
					extendSampleVectorUntil( pSampleVecUntil, new_sample_size );

					/* Do simulations for the newly added observations */
					simulateUnbUntilSampleDTMC( pStateSpace, pSampleVecUntil, old_sample_size,
									pSampleVecUntil->curr_simulation_depth,
									pGoodStates, pTransientStates );
				}

				/* If we need to increase the simulation depth */
				depth_steps_needed = new_simulation_depth - pSampleVecUntil->curr_simulation_depth;
				if( depth_steps_needed > 0 ){
					/* Do simulations for all the observations */
					simulateUnbUntilSampleDTMC( pStateSpace, pSampleVecUntil, FOR_ALL_OBSERVATIONS,
									depth_steps_needed, pGoodStates, pTransientStates );
				}
			ELSE_SAFETY
				printf("ERROR: The sparse matrix is NULL.\n");
                                exit(EXIT_FAILURE);
			ENDIF_SAFETY
		ELSE_SAFETY
			printf("ERROR: The good/transient states are undefined.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: The sample that has to be simulated further is NULL.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}


/**
* Computes the reachability probabilities for the trivial and non-trivial
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
double ** getBSCCReachabilityProbsNumerical( const int n_states, const int numberOfBSCCs,
											const int numberOfTrivialBSCCs,
											bitset ** const ppNonTrivBSCCBitSets,
											const bitset * pTrivialBSCCBitSet,
											const TPFunctNumUnbUntilCTMC pFNumUnbUntilCTMC ) {
	/* Technical variables */
    int index = -1, i = 0;
    bitset * pTmpBitSet = NULL;

    /* Initialise the structure, that stores the probability to reach BSCC i from state j */
    /* The first numberOfTrivialBSCCs elements will contain pointers to arrays of */
    /* reachability probabilities for the trivial BSCCs, in the same order they are */
    /* located in the pTrivialBSCCBitSet array then the probabilities for BSCCs */
    /* stored in the ppNonTrivBSCCBitSets array */
        double ** ppReachProbability = (double**) calloc((size_t) numberOfBSCCs,
                        sizeof (double *));

    /* Create True bitset*/
	bitset *pTRUEBitSet = get_new_bitset( n_states );
	fill_bitset_one( pTRUEBitSet );

    /* For the trivial BSCCs containing Psi states */
    pTmpBitSet = get_new_bitset( n_states );
    while((index = get_idx_next_non_zero( pTrivialBSCCBitSet, index ) ) != -1){
        set_bit_val( pTmpBitSet, index, BIT_ON );
        /* WARNING: We have to be careful here because the state space is an implicit parameter */
        /* We use the state space accessible from the get_state_space() function of runtime.h */
        ppReachProbability[i] = pFNumUnbUntilCTMC( pTRUEBitSet, pTmpBitSet );
        set_bit_val( pTmpBitSet, index, BIT_OFF );
        i++;
    }
    free_bitset(pTmpBitSet);
    pTmpBitSet = NULL;

    /* For the non-trivial BSCCs containing Psi states */
    for(i = numberOfTrivialBSCCs;i < numberOfBSCCs;i++){
        ppReachProbability[i] = pFNumUnbUntilCTMC( pTRUEBitSet, ppNonTrivBSCCBitSets[i - numberOfTrivialBSCCs] );
    }

    /* Clean the unneeded memory */
	free_bitset( pTRUEBitSet ); pTRUEBitSet = NULL;

    return ppReachProbability;
}

/**
 * This function allows to find trivial and non-trivial BSCCs that contain \Psi states
 * @param pStateSpace the state space
 * @param pPsiBitSet the \Psi states
 * @param ppNonTrivBSCCStates the pointer to an array that by the end of the function
 *								will store the array of non-trivial BSCC states.
 *								If the pointer is NULL then the result is not computed.
 * @param ppTrivialBSCCBitSet the pointer to a bitset pointer, it will be initialised with
 *								a pointer to a bitset containing states of the trivial BSCCa
 * @param pNumberOfNonTrivBSCCs this is the pointer to a variable that will store the number
 *								of non-trivial BSCCs
 * @param pNumberOfTrivialBSCCs this is the pointer to a variable that will store the number
 *								of trivial BSCCs
 * @param pNumberOfBSCCs this is the pointer to a variable that will store the total number
 *								of BSCCs
 * @param pFAllowedBSCCSearchCTMC the pointer to a function that allows to search for BSCCs
 * @return the array of bitset pointers, each of which is contains the states of the
 *			corresponding non-trivial BSCC.
 */
bitset ** getTrivialAndNonTrivialPsiBSCC( const sparse * pStateSpace, const bitset * pPsiBitSet, int ** ppNonTrivBSCCStates,
											bitset ** ppTrivialBSCCBitSet, int * pNumberOfNonTrivBSCCs, int *pNumberOfTrivialBSCCs,
											int * pNumberOfBSCCs, const TPFunctAllowedBSCCSearchCTMC pFAllowedBSCCSearchCTMC ) {
	/* Some technical variables */
	int i = 0;
	/* The structure, that stores the set of states belonging to non-trivial BSCC with array index i */
	bitset ** ppNonTrivBSCCBitSets = NULL;
	/* The bitset indicating the Psi states belonging to trivial BSCCs */
	bitset * pTrivialBSCCBitSet = NULL;
	/* This vvariable will store the number of the non-trivial BSCCs */
	int numberOfNonTrivBSCCs = 0;
	/* This variable will hold the number of trivial BSCCs */
	int numberOfTrivialBSCCs = 0;
	/* This variable will store the total number of BSCCs */
	int numberOfBSCCs = 0;
	/*This variable will store an array of indexes of states belonging to non-trivial BSCCs*/
	int * pNonTrivBSCCStates = NULL;

    /* Find BSCCs that contain Psi states, divided into the set of trivial */
    /* and non-trivial BSCCs */
	pTrivialBSCCBitSet = pFAllowedBSCCSearchCTMC( pStateSpace, pPsiBitSet, &ppNonTrivBSCCBitSets, &numberOfNonTrivBSCCs);
    numberOfTrivialBSCCs = count_non_zero( pTrivialBSCCBitSet );
    /* Compute the total number of BSCCs that contain Psi states */
    numberOfBSCCs = numberOfNonTrivBSCCs + numberOfTrivialBSCCs;

    /* If we need to return the state of all the non-trivial \Psi BSCCs */
    if( ppNonTrivBSCCStates != NULL ) {
		/* Create bitset, that in the end will hold all states belonging to non-trivial BSCCs */
                bitset * pNonTrivBSCCStatesBitSet = get_new_bitset(
                                                        mtx_rows(pStateSpace));
		/* For the non-trivial BSCCs containing Psi states */
		for(i = 0; i < numberOfNonTrivBSCCs; i++){
			/* Create a bitset that holds all states belonging to non-trivial BSCCs */
			or_result( ppNonTrivBSCCBitSets[i], pNonTrivBSCCStatesBitSet );
		}

		/* Fill the array with the indices of states belonging to non-trivial BSCCs */
		pNonTrivBSCCStates = count_set( pNonTrivBSCCStatesBitSet );
		(*ppNonTrivBSCCStates) = pNonTrivBSCCStates;

		/* Clean memory */
		free_bitset( pNonTrivBSCCStatesBitSet ); pNonTrivBSCCStatesBitSet = NULL;
    }

    /* Initialise the remaining return argument values and return the result */
    (*ppTrivialBSCCBitSet) = pTrivialBSCCBitSet;
    (*pNumberOfNonTrivBSCCs) = numberOfNonTrivBSCCs;
    (*pNumberOfTrivialBSCCs) = numberOfTrivialBSCCs;
    (*pNumberOfBSCCs) = numberOfBSCCs;

	return ppNonTrivBSCCBitSets;
}

/**
* This method allows to initialized ppBSCCExistsReachableStates and ppBSCCAlwaysReachableStates.
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
void getReachabilityTrivialAndNonTrivialPsiBSCC( const sparse * pStateSpace, const int numberOfBSCCs,
													bitset * const pTrivialBSCCBitSet,
													bitset ** const ppNonTrivBSCCBitSets,
													const int numberOfNonTrivBSCCs,
													bitset ** const ppBSCCExistsReachableStates,
													bitset ** const ppBSCCAlwaysReachableStates,
													const TPFunctExistUnbUntil pExistsUnbUntilCTMC,
													const TPFunctAlwaysUnbUntil pAlwaysUnbUntilCTMC ) {
	/* Some internal indices and counters */
        int state_index_local, index, current_index;
	/* The number of trivial BSCCs */
	const int numberOfTrivBSCCs = ( numberOfBSCCs - numberOfNonTrivBSCCs );
	/* The state space dimension */
        const int n_states = mtx_rows(pStateSpace);
	/* The temporary bitsets needed for the trivial BSCCs */
	bitset * pTmpBitSet = get_new_bitset( n_states );
	bitset * pTRUEBitSet = get_new_bitset( n_states );
	fill_bitset_one( pTRUEBitSet );

	/* First get always and exists reachability analysis for the trivial BSCCs */
        index = 0; state_index_local = state_index_NONE;
        while( (state_index_local = get_idx_next_non_zero(pTrivialBSCCBitSet,
                                state_index_local)) != state_index_NONE )
        {
                /*printf("GETTING A/E REACHABILITY SETS FOR A TRIVIAL PSI BSCC "
                        "%d STATE: %d\n", index+1, state_index_local+1); */
        set_bit_val(pTmpBitSet, state_index_local, BIT_ON);
	ppBSCCExistsReachableStates[index] = pExistsUnbUntilCTMC( pStateSpace, pTRUEBitSet, pTmpBitSet );
	ppBSCCAlwaysReachableStates[index] = pAlwaysUnbUntilCTMC( pStateSpace, pTRUEBitSet, pTmpBitSet,
																	ppBSCCExistsReachableStates[index] );
        set_bit_val(pTmpBitSet, state_index_local, BIT_OFF);
	index++;
    }
    free_bitset(pTmpBitSet); pTmpBitSet = NULL;

	/* Second get always and exists reachability analysis for the non-trivial BSCCs  */
    current_index = numberOfTrivBSCCs;
    for( index = 0; index < numberOfNonTrivBSCCs; index++, current_index++ ) {
		/* printf("GETTING A/E REACHABILITY SETS FOR A NON TRIVIAL PSI "
                        "BSCC %d (index + numberOfTrivBSCCs) STATES:", index+1);
		print_bitset_states( ppNonTrivBSCCBitSets[index] ); printf("\n"); */
	ppBSCCExistsReachableStates[current_index] = pExistsUnbUntilCTMC( pStateSpace, pTRUEBitSet, ppNonTrivBSCCBitSets[index] );
	ppBSCCAlwaysReachableStates[current_index] = pAlwaysUnbUntilCTMC( pStateSpace, pTRUEBitSet, ppNonTrivBSCCBitSets[index],
																			ppBSCCExistsReachableStates[current_index] );
    }

    /* Free the remaining data */
    free_bitset(pTRUEBitSet); pTRUEBitSet = NULL;
}

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
int extractReachabilityProbs( const int initial_state, const int numberOfTrivialBSCCs, const int numberOfNonTrivBSCCs,
										double ** ppReachProbability, double ** ppBSCCReachProbability,
										bitset ** ppReachableTrivialBSCCBitSet, bitset ** ppReachableNonTrivBSCCBitSet,
										int * pNumberOfReachableTrivialBSCCs, int *pNumberOfReachableNonTrivBSCCs ) {
	/* Some internal counters */
	int i;
	/* The pBSCCReachProbability will store all the probabilities of reaching */
	/* all the BSCCs from the given initial state. Some of them are going to be zero */
	/* The indexes of non zero elements of pBSCCReachProbability are going to be stored */
	/* in the next two bitsets. Note that they are indexed independently, i.e. the indexes */
	/* of both bitsets will start from zero, and to get the pBSCCReachProbability index from */
	/* pReachableNonTrivBSCCBitSet we will need to add numberOfTrivialBSCCs to it */
	/* The bitset indicating the trivial BSCCs containing Psi states, which are also reachable */
	bitset * pReachableTrivialBSCCBitSet = get_new_bitset( numberOfTrivialBSCCs );
	/* The bitset indicating the non-trivial BSCCs containing Psi states, which are also reachable */
	bitset * pReachableNonTrivBSCCBitSet = get_new_bitset( numberOfNonTrivBSCCs );
	/* The number of non-trivial BSCCs, that contain Psi */
	/* states and are reachable from the initial state */
	int numberOfReachableNonTrivBSCCs = 0;
	/* The number of trivial BSCCs, that contain Psi states and are reachable from the */
	/* initial state */
	int numberOfReachableTrivialBSCCs = 0;
	/* The number of BSCCs containing Psi states */
	const int numberOfBSCCs = numberOfNonTrivBSCCs + numberOfTrivialBSCCs;
	/* The array holding the reachability probabilities from the initial state */
        double *pBSCCReachProbability = (double*) calloc((size_t) numberOfBSCCs,
                        sizeof(double));

	/* Having reachability probabilities in place, we know what BSCCs, containing Psi states, */
	/* are reachable with non-zero probability. */
	/* WARNING: There might be problems with precision, i.e. 0.0 probability might be represented */
        /* as a very very small, but non-zero value, like
           0.000000000000000000000000000234234234. This */
	/* will have to be checked, at present we assume that it JUST works. */

	/* The BSCC reachability probabilities for every BSCC i*/
	/* For the trivial BSCCs found */
	for( i = 0; i < numberOfTrivialBSCCs; i++ ){
		/* If BSCC i is reachable from the initial state, increment reachable BSCC counter */
		if( ( pBSCCReachProbability[i] = ppReachProbability[i][initial_state] ) > 0.0 ){
			numberOfReachableTrivialBSCCs++;
			set_bit_val( pReachableTrivialBSCCBitSet, i, BIT_ON );
		}
	}

	/* For the non-trivial BSCCs found */
	for( i = numberOfTrivialBSCCs; i < numberOfBSCCs; i++ ){
		/* If non-trivial BSCC i is reachable from the initial state, increment reachable non-trivial */
		/* BSCC counter */
		if( ( pBSCCReachProbability[i] = ppReachProbability[i][initial_state] ) > 0.0 ){
			numberOfReachableNonTrivBSCCs++;
			set_bit_val(pReachableNonTrivBSCCBitSet, i - numberOfTrivialBSCCs, BIT_ON);
		}
	}

	/* Return the explicit and implicit results */
	(* ppBSCCReachProbability) = pBSCCReachProbability;
	(* ppReachableTrivialBSCCBitSet) = pReachableTrivialBSCCBitSet;
	(* ppReachableNonTrivBSCCBitSet) = pReachableNonTrivBSCCBitSet;
	(* pNumberOfReachableTrivialBSCCs) = numberOfReachableTrivialBSCCs;
	(* pNumberOfReachableNonTrivBSCCs) = numberOfReachableNonTrivBSCCs;

	return (numberOfReachableTrivialBSCCs + numberOfReachableNonTrivBSCCs);
}

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
BOOL extractReachabilityData( const int initial_state, const int numberOfTrivialBSCCs, const int numberOfNonTrivBSCCs,
                const bitset * const * const ppBSCCExistsReachableStates,
                const bitset * const * const ppBSCCAlwaysReachableStates,
							bitset ** ppReachableTrivialBSCCBitSet, bitset ** ppReachableNonTrivBSCCBitSet,
							int * pNumberOfExistsReachableTrivialBSCCs, int * pNumberOfExistsReachableNonTrivBSCCs ) {
	/* Internal storage variables */
	bitset * pReachableTrivialBSCCBitSet;
	bitset * pReachableNonTrivBSCCBitSet;
	BOOL isAlwaysReachableBSCC = FALSE;
	int numberOfExistsReachableTrivialBSCCs = 0;
	int numberOfExistsReachableNonTrivBSCCs = 0;
	int index;

	/* Process the trivial \Psi BSCCs first, these are the 0 to numberOfTrivialBSCCs */
	/* elements in ppBSCCAlwaysReachableStates and ppBSCCExistsReachableStates */
	pReachableTrivialBSCCBitSet = get_new_bitset( numberOfTrivialBSCCs );
	/* printf("GETTING THE TRIVIAL BSCCs (count = %d) REACHABLE FROM THE STATE %d\n", numberOfTrivialBSCCs, initial_state); */
	for( index = 0; index < numberOfTrivialBSCCs; index++ ){
		/* If the trivial BSCC index is some times reachable from state initial_state */
		if( get_bit_val( ppBSCCExistsReachableStates[index], initial_state ) ) {
			/* Set this BSCCs on */
			set_bit_val( pReachableTrivialBSCCBitSet, index, BIT_ON );
			numberOfExistsReachableTrivialBSCCs++;
			/* If the trivial BSCC index is always reachable from state initial_state */
			if( get_bit_val( ppBSCCAlwaysReachableStates[index], initial_state ) ) {
				/* There can be no other BSCCs reachable from the state so we quit */
				isAlwaysReachableBSCC = TRUE;
				break;
			}
		}
	}

        /* Process the non-trivial BSCCs second, these are the
           numberOfTrivialBSCCs to numberOfNonTrivBSCCs */
	/* elements in ppBSCCAlwaysReachableStates and ppBSCCExistsReachableStates */
	pReachableNonTrivBSCCBitSet = get_new_bitset( numberOfNonTrivBSCCs );
	if( !isAlwaysReachableBSCC  ) {
		/* printf("GETTING THE NON-TRIVIAL BSCCs (count = %d) REACHABLE FROM THE STATE %d\n", numberOfTrivialBSCCs, initial_state); */
		/* If there is no always reachable trivial BSCC then we can process the non-trivial ones, if any */
		for( index = numberOfTrivialBSCCs; index < ( numberOfTrivialBSCCs + numberOfNonTrivBSCCs ); index++ ) {
			/* If the trivial BSCC index is some times reachable from state initial_state */
			if( get_bit_val( ppBSCCExistsReachableStates[index], initial_state ) ) {
				/* Set this BSCCs on */
				set_bit_val( pReachableNonTrivBSCCBitSet, index - numberOfTrivialBSCCs, BIT_ON );
				numberOfExistsReachableNonTrivBSCCs++;
				/* If the non-trivial BSCC index is always reachable from state initial_state */
				if( get_bit_val( ppBSCCAlwaysReachableStates[index], initial_state ) ){
					/* There can be no other BSCCs reachable from the state so we quit */
					isAlwaysReachableBSCC = TRUE;
					break;
				}
			}
		}
	}

	/* Set the return parameter values */
	(* ppReachableTrivialBSCCBitSet) = pReachableTrivialBSCCBitSet;
	(* ppReachableNonTrivBSCCBitSet) = pReachableNonTrivBSCCBitSet;
	(* pNumberOfExistsReachableTrivialBSCCs) = numberOfExistsReachableTrivialBSCCs;
	(* pNumberOfExistsReachableNonTrivBSCCs) = numberOfExistsReachableNonTrivBSCCs;

	return isAlwaysReachableBSCC;
}

/**
 * This function allows to search for all transient states in the model
 * TODO: May be it's current performance is not optimal.
 * @param pStateSpace the state space
 * @param pFAllowedBSCCSearchCTMC the function that allows to search for BSCCs
 * @return the pointer to the bitset (size of state space) that contains all transient states
 */
bitset * getTransientStates( const sparse * pStateSpace, const TPFunctAllowedBSCCSearchCTMC pFAllowedBSCCSearchCTMC ) {
	/* Internal couters and indices */
	int i;
	/* The state space size */
        const int n_states = mtx_rows(pStateSpace);
	/* The structure, that stores the set of states belonging to non-trivial BSCC with array index i */
	bitset ** ppNonTrivBSCCBitSets = NULL;
	/* This variable will store the number of the non-trivial BSCCs */
	int numberOfNonTrivBSCCs = 0;

	/* The bitset that will contain all the transient states */
	bitset * pTransientStatesBitSet = NULL;

	/* The bitset with all of the states */
	bitset *pAllStatesBitSet = get_new_bitset( n_states );
	fill_bitset_one( pAllStatesBitSet );

	/* After this statement pTransientStatesBitSet will store trivial BSCCs */
	pTransientStatesBitSet = pFAllowedBSCCSearchCTMC( pStateSpace, pAllStatesBitSet, &ppNonTrivBSCCBitSets, &numberOfNonTrivBSCCs);
	/* Deallocate unneeded memory */
	free_bitset( pAllStatesBitSet );

	/* Now we add non-trivial BSCC states by using in-place or operation and deallocate memory */
	for( i = 0; i < numberOfNonTrivBSCCs; i ++ ) {
		or_result( ppNonTrivBSCCBitSets[i] , pTransientStatesBitSet );
		free_bitset( ppNonTrivBSCCBitSets[i] );
	}
	free( ppNonTrivBSCCBitSets ); ppNonTrivBSCCBitSets = NULL;

	/* Do in-place of the set that now contains all BSCC states */
	not_result( pTransientStatesBitSet );

	/* Return the set with all transient states */
	return pTransientStatesBitSet;
}
