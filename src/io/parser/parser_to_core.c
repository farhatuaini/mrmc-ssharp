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
*	Source description: This is an intermediate interface between
*	the parser and the core model checkier.
*/

#include "parser_to_core.h"

#include "core_to_core.h"
#include "simulation.h"
#include "help.h"

#include "runtime.h"

#include <time.h>
#include <math.h>

/*******************************************************************/
/*******************This is the timer implementation****************/
/*******************************************************************/

/* This variable stores the time when model checking has started */
static clock_t start_time = 0;
static BOOL timer_started = FALSE;

/**
* Start the timer.
* NOTE: Here the start of the timer should be done only in case we have an atomic
* proposition, 'tt' or 'ff' formula, because model checking starts from the
* leaves.
* WARNING: BE CAREFUL USING THIS METHOD. THERE SHOULD BE NO CONCURRENT USE OF
* THE TIMER FROM DIFFERENT ALGORITHMS AND SUCH.
* For example it is allowed to use this method for monitoring elapsed time of
* model checking and
* formula independent lumping, which do not overlap in time!
*/
void startTimer(void) {
	/* If timer has not been started yet */
        if ( ! timer_started ) {
                timer_started = TRUE;
                start_time = clock();
	}
}

/**
* Stop the timer and return the elapsed time in milliseconds.
* NOTE: Here the stop of the timer should be done only in case the top level
* state formula procedure is called. I.e. stateformula NEWLINE.
* WARNING: BE CAREFUL USING THIS METHOD. THERE SHOULD BE NO CONCURRENT USE OF
* THE TIMER FROM DIFFERENT ALGORITHMS AND SUCH.
* For example it is allowed to use this method for monitoring elapsed time of
* model checking and
* formula independent lumping, which do not overlap in time!
* @return elapsed time in milliseconds
*/
unsigned long stopTimer(void) {
        clock_t stop_time = clock();

	/* If timer has been started */
        if ( ! timer_started ) {
                return 0L;
        }

        timer_started = FALSE;
        return ((unsigned long) (stop_time - start_time) * 1000
                        + CLOCKS_PER_SEC / 2) / CLOCKS_PER_SEC;
}

/*******************************************************************/
/************The general purpose methods of the interface***********/
/*******************************************************************/

/**
* This procedure does printing of the error messages stored in the help.h file.
* @param help_msg_type the help message type, one of:
*	HELP_GENERAL_MSG_TYPE, HELP_COMMON_MSG_TYPE, HELP_REWARDS_MSG_TYPE,
*	HELP_SIMULATION_MSG_TYPE, HELP_LOGIC_MSG_TYPE.
*/
inline void printHelpMessage(const int help_msg_type){
	switch( help_msg_type ){
		case HELP_GENERAL_MSG_TYPE:
			printf("%s%s", HELP_GENERAL_MSG1, HELP_GENERAL_MSG2);
			break;
		case HELP_COMMON_MSG_TYPE:
			printf("%s%s", HELP_COMMON_MSG1, HELP_COMMON_MSG2);
			break;
		case HELP_REWARDS_MSG_TYPE:
			printf("%s", HELP_REWARDS_MSG);
			break;
		case HELP_SIMULATION_MSG_TYPE:
			printf("%s%s%s%s", HELP_SIMULATION_MSG1, HELP_SIMULATION_MSG2, HELP_SIMULATION_MSG3, HELP_SIMULATION_MSG4);
			break;
		case HELP_LOGIC_MSG_TYPE:
			switch( isRunMode(ANY_MODEL_MODE) ){
				case CTMC_MODE:
					printf("%s",HELP_CSL_MSG);
					break;
				case DTMC_MODE:
					printf("%s",HELP_PCTL_MSG);
					break;
				case DMRM_MODE:
					printf("%s",HELP_PRCTL_MSG);
					break;
				case CMRM_MODE:
					printf("%s",HELP_CSRL_MSG);
					break;
				case CTMDPI_MODE:
					printf("%s",HELP_RESTRICTED_CSL_MSG);
					break;
				default :
					printf("ERROR: An unexpected logic type.\n");
                                        exit(EXIT_FAILURE);
			}
			break;
		default :
			printf("ERROR: An unexpected help-message type %d.\n", help_msg_type);
                        exit(EXIT_FAILURE);
	}
}

/**
* Prints the formula tree
*/
void printResultFormulaTree(void) {
	PTFTypeRes pFTypeRes = (PTFTypeRes) get_formula_tree_result();
	if( pFTypeRes != NULL){
		dumpFormulaTree( pFTypeRes );
	} else {
		printf("WARNING: There is NO formula tree to print.\n");
	}
}

/**
* This method is used for printing state probability if any.
* In case get_result_probs_size() == 0, get_result_probs() == NULL
* or the index is out of bounds, the error messages are printed.
* @param index: the state index >= 1
*/
inline void printResultingStateProbability( const int index ){
	PTFTypeRes pFTypeRes = (PTFTypeRes) get_formula_tree_result();
	PTCompStateF pCompStateF; PTFTypeRes pFTypeResSubForm;
	/* Convert the user-level state index into the internal index */
	const int internal_state_index = index - 1;

	if( ( pFTypeRes != NULL ) && ( pFTypeRes->formula_type == COMPARATOR_SF ) ){
		pCompStateF = (PTCompStateF) pFTypeRes;
		pFTypeResSubForm = (PTFTypeRes) pCompStateF->unary_op.pSubForm;
		/* Check if the simulations are done here */
		if( pFTypeResSubForm->doSimHere ){
			/* WARNING: Here we assume that the conf. int. borders are computed exactly */
			print_state_prob( internal_state_index, pFTypeResSubForm->prob_result_size,
						pFTypeResSubForm->pProbCILeftBorder, LEFT_CI_RESULT_STR, 0.0, NULL );
			print_state_prob( internal_state_index, pFTypeResSubForm->prob_result_size,
						pFTypeResSubForm->pProbCIRightBorder, RIGHT_CI_RESULT_STR, 0.0, NULL );
		} else {
			print_state_prob( internal_state_index, pFTypeResSubForm->prob_result_size,
					pFTypeResSubForm->pProbRewardResult, LEFT_RIGHT_RESULT_STR,
					pFTypeResSubForm->error_bound, pFTypeResSubForm->pErrorBound );
		}
	}else{
		printf("WARNING: There are NO results to print.\n");
	}
}

/**
* This method is used for printing state satisfyability result if any.
* I.e. if the state with the given index satisfyes the previously checked
* stateformula. In case get_result_bitset() == NULL or the index is out
* of bounds, the error messages are printed.
* @param index: the state index >= 1
*/
inline void printResultingStateSatisfyability( const int index ){
	PTFTypeRes pFTypeRes = (PTFTypeRes) get_formula_tree_result();

	if( pFTypeRes != NULL ){
		if( pFTypeRes->doSimHere || pFTypeRes->doSimBelow ){
			print_state( pFTypeRes->pYesBitsetResult, index - 1, YES_STATES_STR );
			print_state( pFTypeRes->pNoBitsetResult, index - 1, NO_STATES_STR );
		} else {
			print_state( pFTypeRes->pYesBitsetResult, index - 1, YES_NO_STATES_STR );
		}
	}else{
		printf("WARNING: There are NO results to print.\n");
	}
}

/**
* This procedure sets the results of the stateformula into the runtime.c
* using the set_result_bitset(void *).
* This method is also used for printing the resulting probabilities if
* any and the set of states that satisfy the formula.
* NOTE: Printing is done only in case when the call of isPrintingOn()
* returns TRUE, get_result_probs() != NULL or get_result_bitset() != NULL
*/
void printFormulaResults(void) {
	PTFTypeRes pFTypeRes = (PTFTypeRes) get_formula_tree_result();
	if( isPrintingOn() && ( pFTypeRes != NULL ) ){
		printFormulaBitsetResultParams(pFTypeRes);
	}
}

/**
* This method is responsibe for the actual bottom-up formula-tree traversal and
* model checking. It also monitores the possibility of employing the simulation
* engine, i.e. on every branch of the formula tree only one formula can be model
* checked using simulations. The simulations in case of nesting formulas is not
* supported because it seems to be inefficient. Note that we only do simulations
* for L, S, U, operators of CSL (planned for PCTL).
* @param pFormulaTreeRootNode the root node of the formula tree
*/
static inline void doFormulaTreeTraversalAndModelChecking( void * pFormulaTreeRootNode ){
	doFormulaTreeTraversal( pFormulaTreeRootNode, FALSE, TRUE, FALSE, NULL, NULL, modelCheckAtomicFormula,
				modelCheckUnaryOperator, modelCheckLongSteadyFormula,
				modelCheckBinaryOperator, modelCheckComparatorFormula,
				modelCheckPureRewardFormula, modelCheckNextFormula,
				modelCheckUntilFormula);
}

/***************************************************************************/
/*************THE METHODS NEEDED TO CHECK THE FORMULA TREE FOR**************/
/*******************WHERE WE NEED OR CAN USE SIMULATION*********************/
/***************************************************************************/

static inline BOOL markAtomicFormulaSim( BOOL UNUSED(before), PTAtomicF UNUSED(pAtomicF) ){
	/* The are no subformulas */
	return FALSE;
}

static inline BOOL markUnaryOperatorSim( BOOL UNUSED(before), PTUnaryOp UNUSED(pUnaryOp) ){
	/* The subformulas can be simulated, so keep looking */
	return TRUE;
}

static inline BOOL markLongSteadyFormulaSim( BOOL UNUSED(before), PTLongSteadyF pLongSteadyF ){
	PTFTypeRes pFTypeRes = (PTFTypeRes) pLongSteadyF;
	pFTypeRes->doSimHere = TRUE;
	/* The subformulas should not be simulated, so stop looking */
	return FALSE;
}

static inline BOOL markBinaryOperatorSim( BOOL UNUSED(before), BOOL UNUSED(between), PTBinaryOp UNUSED(pBinaryOp) ){
	/* The subformulas can be simulated, so keep looking */
	return TRUE;
}

static inline BOOL markComparatorFormulaSim( BOOL UNUSED(before), PTCompStateF UNUSED(pCompStateF) ){
	/* The subformulas can be simulated, so keep looking */
	return TRUE;
}

static inline BOOL markPureRewardFormulaSim( BOOL UNUSED(before), PTPureRewardF UNUSED(pPureRewardF) ){
	BOOL isRecNeeded = FALSE;
	/* In case of the one initial-state we obviously can not apply */
	/* simulation to the subformula of the next operator! */
	/* For all initial-state simulation we could, but only in case */
	/* for every initial state we guarantee the definite answer (TT,FF). */
	/* Since we can not provide that, we have to exclude this possibility. */
	/*if( isSimOneInitState() ){
		isRecNeeded = FALSE;
	}*/
	return isRecNeeded;
}

static inline BOOL markNextFormulaSim( BOOL UNUSED(before), PTNextF UNUSED(pNextF) ){
	BOOL isRecNeeded = FALSE;
	/* In case of the one initial-state we obviously can not apply */
	/* simulation to the subformula of the next operator! */
	/* For all initial-state simulation we could, but only in case */
	/* for every initial state we guarantee the definite answer (TT,FF). */
	/* Since we can not provide that, we have to exclude this possibility. */
	/*if( isSimOneInitState() ){
		isRecNeeded = FALSE;
	}*/
	return isRecNeeded;
}

static inline BOOL markUntilFormulaSim( BOOL UNUSED(before), BOOL UNUSED(between), PTUntilF pUntilF ){
	PTFTypeRes pFTypeRes = (PTFTypeRes) pUntilF;
	/* Distinguish between the cases where we know hod to do */
	/* simulation and where we do not. */
	switch(pUntilF->binary_op.binary_type){
		case UNTIL_PF_UNB:
		case UNTIL_PF_TIME:
			pFTypeRes->doSimHere = TRUE;
			break;
		case UNTIL_PF_TIME_REWARD:
			/* We do not know how to simulate this one yet */
			break;
		default:
			printf("ERROR: An unknown type '%d' of UNTIL_PF operator.\n", pUntilF->binary_op.binary_type);
                        exit(EXIT_FAILURE);
	}
	/* The subformulas should not be simulated, so stop looking */
	return FALSE;
}

/***************************************************************************/
/************THE METHODS NEEDED TO INDICATE THE FORMULA BRANCES*************/
/*******************WHERE WE NEED OR CAN USE SIMULATION*********************/
/***************************************************************************/

static inline void markTargetFromUnaryOpSim( void * pCurrentNode, PTUnaryOp pUnaryOp ){
	PTFTypeRes pFTypeRes = (PTFTypeRes) pCurrentNode;
	PTFTypeRes pFTypeResSubForm = (PTFTypeRes) pUnaryOp->pSubForm;

	/* If it has to be done below, then we mark this branch */
	pFTypeRes->doSimBelow = pFTypeResSubForm->doSimHere || pFTypeResSubForm->doSimBelow;

	/* If we are going to have simulation with one initial state we better remember it. */
	if( ( pFTypeRes->doSimHere || pFTypeRes->doSimBelow ) && isSimOneInitState() ){
		pFTypeRes->initial_state = getSimInitialState();
		pFTypeRes->isSimOneInitState = TRUE;
	}
}

static inline void markTargetFromBinaryOpSim( void * pCurrentNode, PTBinaryOp pBinaryOp ){
	PTFTypeRes pFTypeRes = (PTFTypeRes) pCurrentNode;
	PTFTypeRes pFTypeResSubFormL = (PTFTypeRes) pBinaryOp->pSubFormL;
	PTFTypeRes pFTypeResSubFormR = (PTFTypeRes) pBinaryOp->pSubFormR;

	/* If it had to be done below, then we mark this branch */
	pFTypeRes->doSimBelow = pFTypeResSubFormL->doSimHere || pFTypeResSubFormL->doSimBelow ||
					pFTypeResSubFormR->doSimHere || pFTypeResSubFormR->doSimBelow;

	/* If we are going to have simulation with one initial state we better remember it. */
	if( ( pFTypeRes->doSimHere || pFTypeRes->doSimBelow ) && isSimOneInitState() ){
		pFTypeRes->initial_state = getSimInitialState();
		pFTypeRes->isSimOneInitState = TRUE;
	}
}

static inline BOOL markBranchAtomicFormulaSim( BOOL UNUSED(before), PTAtomicF UNUSED(pAtomicF) ){
	/* The end of branch - do nothing */
	return FALSE;
}

static inline BOOL markBranchUnaryOperatorSim( BOOL UNUSED(before), PTUnaryOp pUnaryOp ){
	markTargetFromUnaryOpSim( pUnaryOp, pUnaryOp );
	return FALSE;
}

static inline BOOL markBranchLongSteadyFormulaSim( BOOL before, PTLongSteadyF pLongSteadyF ){
	/* WARNING: Works as long as TLongSteadyF has TUnaryOp as the first field */
        markBranchUnaryOperatorSim(before, &pLongSteadyF->unary_op);
	return FALSE;
}

static inline BOOL markBranchBinaryOperatorSim( BOOL UNUSED(before), BOOL UNUSED(between), PTBinaryOp pBinaryOp ){
	markTargetFromBinaryOpSim( pBinaryOp, pBinaryOp );
	return FALSE;
}

static inline BOOL markBranchComparatorFormulaSim( BOOL UNUSED(before), PTCompStateF pCompStateF ){
	markTargetFromUnaryOpSim( pCompStateF, &(pCompStateF->unary_op) );
	return FALSE;
}

static inline BOOL markBranchPureRewardFormulaSim( BOOL UNUSED(before), PTPureRewardF pPureRewardF ){
	markTargetFromUnaryOpSim( pPureRewardF, &(pPureRewardF->unary_op) );
	return FALSE;
}

static inline BOOL markBranchNextFormulaSim( BOOL UNUSED(before), PTNextF pNextF ){
	markTargetFromUnaryOpSim( pNextF, &(pNextF->unary_op) );
	return FALSE;
}

static inline BOOL markBranchUntilFormulaSim( BOOL UNUSED(before), BOOL UNUSED(between), PTUntilF pUntilF ){
	markTargetFromBinaryOpSim( pUntilF, &(pUntilF->binary_op) );
	return FALSE;
}

/***************************************************************************/
/************THE METHODS NEEDED TO DERIVE THE CONFIDENCE LEVELS*************/
/****************FOR THE SUBFORMULAS WHERE USE SIMULATION*******************/
/***************************************************************************/

static inline BOOL propConfFromTragetToUnaryOpSim( void * pCurrentNode, PTUnaryOp pUnaryOp ){
	BOOL searchSubTree = FALSE;
	PTFTypeRes pFTypeRes = (PTFTypeRes) pCurrentNode;
	PTFTypeRes pFTypeResSubForm = (PTFTypeRes) pUnaryOp->pSubForm;
	if( pFTypeResSubForm->doSimHere || pFTypeResSubForm->doSimBelow ){
		/* For negation, braced formula, comparator, pure reward and next */
		/* operators the confidence of the sub-formula stays the same. */
		/* For the negation operator it is not so obvious consider: */
		/* When we model check A using simulation with confidence 1-beta */
		/* then the algorithm provides us the conf. int. we base our TT,FF,NN answers on. */
		/* The definite TT (A),FF (!A) answers are given with the confidence at least 1-beta */
		pFTypeResSubForm->confidence = pFTypeRes->confidence;
		searchSubTree = TRUE;
	}
	return searchSubTree;
}

static inline BOOL propConfAtomicFormulaSim( BOOL UNUSED(before), PTAtomicF UNUSED(pAtomicF) ){
	/* No need to assign any confidence because this is the tree leaf */
	return FALSE;
}

static inline BOOL propConfUnaryOperatorSim( BOOL UNUSED(before), PTUnaryOp pUnaryOp ){
	return propConfFromTragetToUnaryOpSim( pUnaryOp, pUnaryOp );
}

static inline BOOL propConfLongSteadyFormulaSim( BOOL UNUSED(before), PTLongSteadyF UNUSED(pLongSteadyF) ){
	/* If we are here it means that below we definitely do not do any simulations. */
	return FALSE;
}

static inline BOOL propConfBinaryOperatorSim( BOOL UNUSED(before), BOOL UNUSED(between), PTBinaryOp pBinaryOp ){
	double confSubForm;
	BOOL searchSubTree = FALSE;
	PTFTypeRes pFTypeRes = (PTFTypeRes) pBinaryOp;
	PTFTypeRes pFTypeResSubFormL = (PTFTypeRes) pBinaryOp->pSubFormL;
	PTFTypeRes pFTypeResSubFormR = (PTFTypeRes) pBinaryOp->pSubFormR;
	if( ( pFTypeResSubFormL->doSimHere || pFTypeResSubFormL->doSimBelow ) &&
		( pFTypeResSubFormR->doSimHere || pFTypeResSubFormR->doSimBelow ) ){
		/* If the left and the right subformulas need simulation then thigs */
		/* does not depend on the binary operator type, but we just want to */
		/* avoid having an unknown pBinaryOp-binary_type */
		switch( pBinaryOp->binary_type ){
			case BINARY_OP_SF_OR:
			case BINARY_OP_SF_IMPLIES:
			case BINARY_OP_SF_AND:
				/* For all the cases AvB, A^B, A=>B (== !AvB) */
				/* If we have independently concluded A and B with conf 1-alpha each */
				/* Then the implied conclusions AvB, A^B, A=>B will have confidence */
				/* (1-alpha)^2 because they are an implication of independent A and B. */
				confSubForm = sqrt(pFTypeRes->confidence);
				pFTypeResSubFormL->confidence = confSubForm;
				pFTypeResSubFormR->confidence = confSubForm;
				break;
			default:
				printf("ERROR: An unknown type of binary operator %d, while conf. level propogation.\n",pBinaryOp->binary_type);
                                exit(EXIT_FAILURE);
		}
		searchSubTree = TRUE;
	}else{
		if( pFTypeResSubFormL->doSimHere || pFTypeResSubFormL->doSimBelow ){
			/* If only the left subformula needs simulation */
			/* Then the confidence is derived as it is */
			pFTypeResSubFormL->confidence = pFTypeRes->confidence;
			searchSubTree = TRUE;
		}else{
			if( pFTypeResSubFormR->doSimHere || pFTypeResSubFormR->doSimBelow ){
				/* If only the right subformula needs simulation */
				/* Then the confidence is derived as it is */
				pFTypeResSubFormR->confidence = pFTypeRes->confidence;
				searchSubTree = TRUE;
			}
		}
	}
	return searchSubTree;
}

static inline BOOL propConfComparatorFormulaSim( BOOL UNUSED(before), PTCompStateF pCompStateF ){
	return propConfFromTragetToUnaryOpSim( pCompStateF, &(pCompStateF->unary_op) );
}

static inline BOOL propConfPureRewardFormulaSim( BOOL UNUSED(before), PTPureRewardF pPureRewardF ){
	return propConfFromTragetToUnaryOpSim( pPureRewardF, &(pPureRewardF->unary_op) );
}

static inline BOOL propConfNextFormulaSim( BOOL UNUSED(before), PTNextF pNextF ){
	return propConfFromTragetToUnaryOpSim( pNextF, &(pNextF->unary_op) );
}

static inline BOOL propConfUntilFormulaSim( BOOL UNUSED(before), BOOL UNUSED(between), PTUntilF UNUSED(pUntilF) ){
	/* If we are here it means that below we definitely do not do any simulations. */
	return FALSE;
}

/**
* This method is responsibe for the bottom-up formula-tree traversal and
* deriving the right confidences for the subformulas, as well as restricting
* the use of simulation engine for nested formulas and in case of single state
* simulations.
* @param pFormulaTreeRootNode the root node of the formula tree
*/
static inline void doFormulaTreeTraversalForSimulations( void * pFormulaTreeRootNode ){
	PTFTypeRes pFTypeRes = (PTFTypeRes) pFormulaTreeRootNode;
	/* Set the confidence level of the root node, i.e. the general */
	/* confidencelevel, as defined in the runtime.c settings */
	pFTypeRes->confidence = getSimGeneralConfidence();
	pFTypeRes->indiff_width = getSupIndifferenceWidth();

	/* First find out and mark the operator nodes that we will simulate */
	doFormulaTreeTraversal( pFormulaTreeRootNode, TRUE, FALSE, FALSE, NULL, NULL, markAtomicFormulaSim,
				markUnaryOperatorSim, markLongSteadyFormulaSim,
				markBinaryOperatorSim, markComparatorFormulaSim,
				markPureRewardFormulaSim, markNextFormulaSim,
				markUntilFormulaSim);
	/* Then mark the tree branches from the before marked operators to the root. */
	/* This is needed to derive the confidence values for model checking subformulas. */
	doFormulaTreeTraversal( pFormulaTreeRootNode, FALSE, TRUE, FALSE, NULL, NULL, markBranchAtomicFormulaSim,
				markBranchUnaryOperatorSim, markBranchLongSteadyFormulaSim,
				markBranchBinaryOperatorSim, markBranchComparatorFormulaSim,
				markBranchPureRewardFormulaSim, markBranchNextFormulaSim,
				markBranchUntilFormulaSim);
	/* Derive the proper confidence values. */
	doFormulaTreeTraversal( pFormulaTreeRootNode, TRUE, FALSE, FALSE, NULL, NULL, propConfAtomicFormulaSim,
				propConfUnaryOperatorSim, propConfLongSteadyFormulaSim,
				propConfBinaryOperatorSim, propConfComparatorFormulaSim,
				propConfPureRewardFormulaSim, propConfNextFormulaSim,
				propConfUntilFormulaSim);
}

/**
* This method is responsible for cleaning the old model-checking formula tree.
*/
void clearOldModelCheckingResults(void) {
	freeFormulaTree( get_formula_tree_result() );
	set_formula_tree_result( NULL );
}

/**
* This method does bottom-up traversal of the formula tree
* along with model checking the formulas. It invokes printing
* the results if needed and is also responsible for outputting
* the model-checking time.
* @param pFormulaTreeRootNode the formula tree root node
*/
inline void modelCheckFormulaTree( void * pFormulaTreeRootNode ){
        unsigned
	long elapsed_time = 0;

	/* Clean the old results */
	clearOldModelCheckingResults();

	/* Star the timer */
	startTimer();

	/* In case of simulation engine on we have to traverse the tree */
	/* in order to assign proper confidence levels to the operator nodes. */
	if( isSimulationOn() ){
		doFormulaTreeTraversalForSimulations( pFormulaTreeRootNode );
	}

	/* Invokes the model checking procedure */
	doFormulaTreeTraversalAndModelChecking( pFormulaTreeRootNode );

	/* Stop the timer */
	elapsed_time = stopTimer();

	/* Set the formula tree with the model checking results in it */
	set_formula_tree_result( pFormulaTreeRootNode );

	/* Print the results if required */
	printFormulaResults();

        printf("The Total Elapsed Model-Checking Time is %lu milli sec(s).\n",
                        elapsed_time);
}
