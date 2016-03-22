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
*	Source description: This is an intermediate interface between
*	the parser and the core model checkier.
*/

#ifndef PARSER_TO_CORE
#define PARSER_TO_CORE

/**
* This method invokes recursive model-checking procedure of the parsed formula
* @param pResultBitset the pointer to the formula tree, that has to be model checked
*/
extern
inline void modelCheckFormula( void* pResultBitset );

/**
* This procedure does printing of the error messages stored in the help.h file.
* @param help_msg_type the help message type, one of:
*	HELP_GENERAL_MSG_TYPE, HELP_COMMON_MSG_TYPE, HELP_REWARDS_MSG_TYPE,
*	HELP_SIMULATION_MSG_TYPE, HELP_LOGIC_MSG_TYPE.
*/
extern
inline void printHelpMessage(const int help_msg_type);

/**
* Prints the formula tree
*/
extern void printResultFormulaTree(void);

/**
* This method is used for printing state probability if any.
* In case get_result_probs_size() == 0, get_result_probs() == NULL
* or the index is out of bounds, the error messages are printed.
* @param index: the state index >= 1
*/

extern void printResultingStateProbability(const int index_local);

/**
* This method is used for printing state satisfyability result if any.
* I.e. if the state with the given index satisfyes the previously checked
* stateformula. In case get_result_bitset() == NULL or the index is out
* of bounds, the error messages are printed.
* @param index: the state index >= 1
*/
extern void printResultingStateSatisfyability(const int index_local);

/**
* This procedure sets the results of the stateformula into the runtime.c
* using the set_result_bitset(void *).
* This method is also used for printing the resulting probabilities if
* any and the set of states that satisfy the formula.
* NOTE: Printing is done only in case when the call of isPrintingOn()
* returns TRUE, get_result_probs() != NULL or get_result_bitset() != NULL
*/
extern void printFormulaResults(void);

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
extern void startTimer(void);

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
extern unsigned long stopTimer(void);

/**
* This method is responsible for cleaning the old model-checking formula tree.
*/
extern void clearOldModelCheckingResults(void);

/**
* This method does bottom-up traversal of the formula tree
* along with model checking the formulas. It invokes printing
* the results if neede and is also responsiblce for outputting
* the model-checking time.
* @param pFormulaTreeRootNode the formula tree root node
*/
extern
inline void modelCheckFormulaTree( void * pFormulaTreeRootNode );

#endif
