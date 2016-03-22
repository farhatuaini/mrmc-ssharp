/**
*	WARNING: Do Not Remove This Section
*
*       $LastChangedRevision: 415 $
*       $LastChangedDate: 2010-12-18 17:21:05 +0100 (Sa, 18. Dez 2010) $
*       $LastChangedBy: davidjansen $
*
*	MRMC is a model checker for discrete-time and continuous-time Markov reward models.
*	It supports reward extensions of PCTL and CSL (PRCTL and CSRL), and allows for the
*	automated verification of properties concerning long-run and instantaneous rewards
*	as well as cumulative rewards.
*
*	Copyright (C) The University of Twente, 2004-2006.
*	Authors: Maneesh Khattri, Ivan Zapreev
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
*       Main contact:
*               Lehrstuhl für Informatik 2, RWTH Aachen University
*               Ahornstrasse 55, 52074 Aachen, Germany
*               E-mail: info@mrmc-tool.org
*
*       Old contact:
*		Formal Methods and Tools Group, University of Twente,
*		P.O. Box 217, 7500 AE Enschede, The Netherlands,
*		Phone: +31 53 4893767, Fax: +31 53 4893247,
*		E-mail: mrmc@cs.utwente.nl
*
*	Source description: Store global variable: state_space, labelling, result_bitset.
*	Uses: DEF: sparse.h, bitset.h, label.h, runtime.h
*		LIB: sparse.c, bitset.c, label.c, runtime.c
*		Definition of runtime - runtime.h
*/

#ifndef RUNTIME_H
#define RUNTIME_H

#include "partition.h"
#include "sparse.h"
#include "mdp_sparse.h"


/* The runtime flag settings */
#define BLANK_MODE      ((unsigned) 0x00000000)
/*This one is used to test whether a model has been set.
*I.e. ANY_MODEL_MODE == CTMC_MODE | DTMC_MODE | DMRM_MODE | CMRM_MODE | CTMDPI_MODE */
#define ANY_MODEL_MODE  ((unsigned) 0x0000008F) /* 00000000 10001111 */
#define CTMC_MODE       ((unsigned) 0x00000001) /* 00000000 00000001 */
#define DTMC_MODE       ((unsigned) 0x00000002) /* 00000000 00000010 */
#define DMRM_MODE       ((unsigned) 0x00000004) /* 00000000 00000100 */
#define CMRM_MODE       ((unsigned) 0x00000008) /* 00000000 00001000 */
#define F_IND_LUMP_MODE ((unsigned) 0x00000010) /* 00000000 00010000 */
#define F_DEP_LUMP_MODE ((unsigned) 0x00000020) /* 00000000 00100000 */
#define SIMULATION_MODE ((unsigned) 0x00000040) /* 00000000 01000000 */
#define CTMDPI_MODE     ((unsigned) 0x00000080) /* 00000000 10000000 */
#define TEST_VMV_MODE   ((unsigned) 0x80000000) /* ????????????????? */

/* ToDo: Make this in one bitset, all of them are just flags */
#define GJ 11  /*Gauss-Jacobi*/
#define GS 12  /*Gauss-Seidel*/
#define US 13  /* uniformization Sericola */
#define UQS 14 /* uniformization Qureshi & Sanders */
#define DTV 15 /* discretization Tijms & Veldman */
#define REC 16 /* recursive version of BSCC search */
#define NON_REC 17 /* non-recursive version of BSCC search */

/* The comparator status */
#define C_LESS 1
#define C_LESS_EQUAL 2
#define C_GREATER 3
#define C_GREATER_EQUAL 4

/* The CTMDPI transient algorithm */
#define CTMDPI_HD_UNI_METHOD 0 /* method for uniform, HD scheduler */
#define CTMDPI_HD_NON_UNI_METHOD 1  /* method for nonuniform, HD scheduler */
#define CTMDPI_HD_AUTO_METHOD 2 /* automatic choice, HD scheduler */

/************************************************************************************/
/******************************THE RUN-MODE ACCESS FUNCTIONS*************************/
/************************************************************************************/

/**
* This method sets the mode in which we run the mcc tool
* It sets a flag to a bitset which can have the following flags
* on or off:  CTMC_MODE, DTMC_MODE, DMRM_MODE, CMRM_MODE etc.
*/
extern void addRunMode(unsigned int flag);

/**
* This method clears part of the mode in which we run the mcc tool
* It clears the flag from a bitset which can have the following flags
* on or off:  CTMC_MODE, DTMC_MODE, DMRM_MODE, CMRM_MODE etc.
*/
extern void clearRunMode(unsigned int flag);

/**
* This method tests the mode in which we run the mcc tool
* It returns true(The logic ID) if the given flag is set
* in a mode bitset which can have the following flags on
* or off: CTMC_MODE, DTMC_MODE, DMRM_MODE, CMRM_MODE etc.
*
* Example: CTMC_MODE == isRunMode(CTMC_MODE) in case we are
* in CTMC_MODE otherwise BLANK_MODE == isRunMode(CTMC_MODE)
*/
extern unsigned int isRunMode(unsigned int flag);

/**
* This method (using the mode local file variable and the ANY_MODEL_MODE
* macro value) tests whether the logic has been already set. Returns
* True if yes, otherwise, False.
*/
extern unsigned int isRunModeSet(void);

/************************************************************************************/
/****************************THE FORMULA-TREE ACCESS FUNCTIONS***********************/
/************************************************************************************/

/**
* Sets the formula tree after it has been modelchecked
* @param void* the formula tree with the results.
*/
extern void set_formula_tree_result(void *);

/**
* Gets the formula tree after it has been modelchecked
* @return the formula tree with the results.
*/
extern void * get_formula_tree_result(void);

/************************************************************************************/
/****************************THE STATE-SPACE ACCESS FUNCTIONS************************/
/************************************************************************************/

/**
* This function sets the state space for global access.
* @param space the state space
* WARNING: If space == NULL then the row sums are freed using the free_row_sums() method.
*/
extern void set_state_space(sparse *);

/**
* Get the currently set state space.
*
* WARNING: Can be used only in DTMC, CTMC, DMRM, and CMRM modes!
*
* @return the pointer to the sparse matrix containing the current state.
*/
extern sparse * get_state_space(void);

/**
* Globally access the MDP state space.
* @return the pointer to the sparse matrix of the MDP while im CTMDPI_MODE.
*/
extern NDSparseMatrix * get_mdpi_state_space(void);

/**
* Get the state space size
*
* WARNING: Works for DTMC, CTMC, DMRM, CMRM and SCTMDPI modes!
*
* @return the number of states in the currently used model, i.e.
*		"get_state_space()->rows" or "get_mdpi_state_space()->n"
*	depending on the run-time mode.
*/
extern int get_state_space_size(void);

/**
* Set MDP state space for global access.
* @param mdp pointer to mdpi state space
*/
extern void set_mdpi_state_space(NDSparseMatrix *mdp);

/**
* Frees the MDP state space.
*/
extern void freeNDSparseMatrix(void);

 /**
 * Sets method to be used for CTMDP transient analysis.
 * The method may be one of CTMDPI_UNIFORM_METHOD or
 * CTMDPI_NONUNIFORM_METHOD.
 *
 * @param new_method method to be set
 */
extern void set_method_ctmdpi_transient(int);

extern int get_method_ctmdpi_transient(void);

/************************************************************************************/
/*****************************THE LABELLING ACCESS FUNCTIONS*************************/
/************************************************************************************/

/**
* Sets labelling function for global access using get_labeller.
* @param labelling the labelling function.
*/
extern void set_labeller(const labelling *);

/**
* Access the labelling function.
* @return returns a pointer to the structure labelling.
*/
extern const labelling * get_labeller(void);

/************************************************************************************/
/***************************THE STATE REWARDS ACCESS FUNCTIONS***********************/
/************************************************************************************/

/**
* Set State Rewards
*/
extern void setStateRewards(double * _pRewards);

/**
* Get State Rewards
*/
extern /*@dependent@*/ double * getStateRewards(void);

/**
* Free memory
*/
extern void freeStateRewards(void);

/************************************************************************************/
/**************************IMPULSE REWARDS ACCESS FUNCTIONS**************************/
/************************************************************************************/

/**
* Get the value of pImpulseRewards
* @return impulse rewards
*/
extern const sparse * getImpulseRewards(void);

/**
* Set the value of pImpulseRewards
* @param sparse* _impulse_rewards
*/
extern void setImpulseRewards(sparse *);

/**
* This method is used to free the impulse reward structure.
*/
extern void freeImpulseRewards(void);

/************************************************************************************/
/********************************SIMULATION MODE SETTINGS****************************/
/************************************************************************************/

/**
* This method should be used to set the simulation ON/OFF
* indicator in the runtime settings.
* @param BOOL val TRUE to turn simulation ON.
* Provide FALSE for switching it OFF.
*/
extern void setSimulationStatus(BOOL val);

/**
* This method should be used to get the simulation ON/OFF
* indicator in the runtime settings.
* @return TRUE if the simulation engine is ON, otherwise FALSE.
*/
extern BOOL isSimulationOn(void);

/**
* This method allows to set the initial state for the "one initial state" simulation.
* NOTE: We do not call the setSimInitialState(int) method of simulation.h directly
*	because if case of formula independent lumping the state-space is partitioned.
* @param initial_state the initial state index
*/
extern
inline void setSimInitialStateRuntime( const int user_initial_state );

/************************************************************************************/
/*********************************LUMPING MODE SETTINGS******************************/
/************************************************************************************/

/**
* This method is used for setting the partitioning
* for formula independent lumping
* @param _P the partition for formula independent lumping
*/
extern void setPartition(partition *_P);

/**
* This method is used for getting the partitioning
* for formula independent lumping
* @return The partition that was set before by the
*	  setPartition(partition *) method or NULL.
*/
extern const partition * getPartition(void);

/**
* This method is used for freeing the memory of partition
* for formula independent lumping
*/
extern void freePartition(void);

/**
* Lets us know whether we are working with the statespace after it was lumped.
* @return TRUE if we are working with the lumped markov chain
*/
extern BOOL isFormulaLumpingDone(void);

/**
* Is used to set Formula dependent lumping ON/OFF
* @param _formula_lumping_is_done TRUE if we just did the formula dependent lumping
*					  FALSE if we unlumped the state space etc. back.
*/
extern void setFormulaLumpingDone(BOOL _formula_lumping_is_done);

/**
* Here we get the original state-space size,
* because there could be lumping done before.
* @return the original state-space size
*/
extern int get_original_state_space_size(void);

/************************************************************************************/
/****************************MINOR SETTING ACCESS FUNCTIONS**************************/
/************************************************************************************/

/**
* Gets the value of w
*/
extern double get_w(void);

/**
* Sets the value of w
* @param double the value of w
*/
extern void set_w(double);

/*****************************************************************************
name		: get_d_factor
role		: get the value of d
@param		:
@return         : double: d
******************************************************************************/
extern double get_d_factor(void);

/*****************************************************************************
name		: set_d_factor
role		: set the value of d_factor
@param		: double: _d_factor
remark		:
******************************************************************************/
extern void set_d_factor(double _d_factor);

/*****************************************************************************
name		: get_max_iterations
role		: get the value of max_iterations
@param		:
@return         : int: max_iterations
******************************************************************************/
extern int get_max_iterations(void);

/*****************************************************************************
name		: set_max_iterations
role		: set the value of max_iterations
@param		: int: max_iterations
remark		:
******************************************************************************/
extern void set_max_iterations(int);

/*****************************************************************************
name		: get_underflow
role		: get the value of underflow
@param		:
@return         : double: underflow
******************************************************************************/
extern double get_underflow(void);

/*****************************************************************************
name		: set_underflow
role		: set the value of underflow
@param		: double: underflow
remark		:
******************************************************************************/
extern void set_underflow(double underflow);

/*****************************************************************************
name		: get_overflow
role		: get overflow
@param		:
@return         : double: overflow
******************************************************************************/
extern double get_overflow(void);

/*****************************************************************************
name		: set_overflow
role		: set the value of overflow
@param		: double: overflow
remark		:
******************************************************************************/
extern void set_overflow(double overflow);

/*****************************************************************************
name		: set_method_path
role		: set method_path for global access.
@param          : int: method_path
remark		: should be either GJ or GS
******************************************************************************/
extern void set_method_path(int);

/*****************************************************************************
name		: get_method_path
role		: get method_path
@param		:
@return         : int: method_path
******************************************************************************/
extern int get_method_path(void);

/*****************************************************************************
name		: set_method_steady
role		: set method_steady for global access.
@param          : int: method_steady
remark		: should be either GJ or GS
******************************************************************************/
extern void set_method_steady(int _method_steady);

/*****************************************************************************
name		: get_method_steady
role		: get method_steady
@param		:
@return         : int: method_steady
******************************************************************************/
extern int get_method_steady(void);

/*****************************************************************************
name		: set_method_until_rewards
role		: set the method for the evaluation of time-reward-bounded until.
@param		: int: the method.
remark		:
******************************************************************************/
extern void set_method_until_rewards(int);

/*****************************************************************************
name		: get_method_until_rewards
role		: globally access the result bitset
@return		: int: the method.
remark		:
******************************************************************************/
extern int get_method_until_rewards(void);

/*****************************************************************************
name		: set_error_bound
role		: set error_bound for global access.
@param		: double: the error bound
remark		:
******************************************************************************/
extern void set_error_bound(double);

/*****************************************************************************
name		: get_error_bound
role		: get error_bound.
@param		:
@return         : double: the error bound
******************************************************************************/
extern double get_error_bound(void);

/**
* Set method for the logic comparator.
* @param comp the comparator to be set
* NOTE: Should be one of {GREATER, GREATER_EQUAL, LESS, LESS_EQUAL}
*/
extern void set_comparator(int comp);

/**
* Get method for the logic comparator.
*/
extern int get_comparator(void);


/************************************************************************************/
/******************************THE BSCC SEARCH SETTINGS******************************/
/************************************************************************************/

/**
* Set method for the BSCC search
* @param the method to be set
* NOTE: The method should be either REC or NON_REC
*/
extern void set_method_bscc(int);

/**
* Get method for the BSCC search
* @param the method to be set
* NOTE: The method should be either REC or NON_REC
*/
extern int get_method_bscc(void);

/************************************************************************************/
/***************************THE STEADY-STATE DETECTION SETTINGS**********************/
/************************************************************************************/

/**
* This method is used to set the steady-state detection on and off
* @param _on_off TRUE for setting the steady-state detection on, otherwise FALSE.
*/
extern
void set_ssd(BOOL _on_off);

/**
* Enables steady-state detection for uniformization (CSL)
*/
extern void set_ssd_on(void);

/**
* Disables steady-state detection for uniformization (CSL)
*/
extern void set_ssd_off(void);

/**
* @return TRUE if steady-state detection is enabled, otherwise FALSE
*/
extern BOOL is_ssd_on(void);

/************************************************************************************/
/*********************************THE TECHNICAL FUNCTIONS****************************/
/************************************************************************************/

/**
* This method allows to globally access the row sums of the currently set state space.
* @return an array that contains the sum of ALL row elements of the current
*       state_space matrix. It means that including the diagonal values and
*       this it is not exactly what you have on the diagonal of the Generator
*	matrix! To obtain the i'th diagonal value of the generator matrix you
*	should do (state_space->val[i].diag - row_sums[i])
*/
extern const double * get_row_sums(void);

/**
* This method frees the row_sums array.
* WARNING: It is called automatically if the set_state_space(...) method is called with the NULL parameter
*/
extern void free_row_sums(void);

/************************************************************************************/
/*********************************PRINTING MODE SETTINGS*****************************/
/************************************************************************************/

/**
* This method should be used to set the printing ON/OFF
* indicator in the runtime settings.
* @param BOOL val TRUE to turn printing of the
*				  probabilities and states ON.
*				  Provide FALSE for switching
*				  it OFF.
*/
extern void setPrintingStatus(BOOL val);

/**
* This method should be used to get the printing ON/OFF
* indicator in the runtime settings.
* @return TRUE if the printing of the probabilities and
*	    states is ON, otherwise FALSE.
*/
extern BOOL isPrintingOn(void);

/************************************************************************************/
/**********************************THE PRINTING FUNCTIONS****************************/
/************************************************************************************/

/**
* This method is used for printing the logic name and also extra options such as Lumping mode.
* @param logic: if set to TRUE then the Logic info is printed
* @param extra: if set to TRUE then the Lumping info is printed
*/
extern inline void print_run_mode(BOOL logic, BOOL extra);

/**
* This method is used to print the runtime settings of error
* bounds, maximum number of iterations etc.
*/
extern void print_runtime_info(void);

/**
* This method is used for getting the printing pattern from the error bound.
* NOTE: If we have error 1.231e-10 then the precision
* is 10+1 digits after the decimal point, but not 13!
* @param index the state index, is needed if we want to know the pressision
*		defined by the error: pErrorBounds[index]
*	WARNING: If we are to print an array of values with the given pattern,
*		then use "index <= -1"
* @param error_bound the general error bound, is used if pErrorBounds == NULL
* @param size the size of pErrorBounds
* @param pErrorBounds the array of error bounds for each initial state
* @return the precision for the printing pattern by taking the position
*	    of the first significant digit after the decimal point + 1
* NOTE: if the error bound turns out to be zero, then we return the pattern
* for the error bound provided by the get_error_bound() function.
*/
extern int get_error_bound_precision(const int index, const double error_bound,
                const int size, /*@null@*/ const double * pErrorBounds);

/**
* This method is used for getting the printing pattern from the value.
* NOTE: If we have a value 1.231e-10 then the precision
* is 10+1 digits after the decimal point, but not 13!
* @param value: the value to evaluate
* @return the precision for the printing pattern by taking the position
*	    of the first significant digit after the decimal point + 1
*/
extern inline int get_precision(double value);

/**
* This method is used for printing a bitset which is an
* indicator set of states that satisfy the property
* @param b the bitset which has to be printed
* @param name the bitset name string
*/
extern void print_states(const bitset * b, const char * name);

/**
* This method is used for printing TRUE/FALSE result for a state index.
* TRUE if the state 'index' is in the set, otherwise FALSE.
* This method handles lumping.
* @param pBitset: the bitset which has to be printed
* @param index: the state index in the original state space.
* @param name the bitset name string
*/
extern void print_state( const bitset * pBitset, const int index, const char * name );

/**
* This method is used for printing the errors of state probabilities
* Depending whether there was lumping done or not to the
* state space it does different printing, Note that, since we are to
* print errors, we do nto use any template printing
* @param size the size of the probability-error vector
* @param pErrorValues the vector of probability errors
* @param pName the name string
*/
extern
void print_error_probs( const int size, const double * pErrorValues, const char * pName );

/**
* This method is used for printing the state probabilities
* Depending whether there was lumping done or not to the
* state space it does different printing
* @param size the size of the probabilities vector
* @param probs the vector of probabilities
* @param name the bitset name string
* @param error_bound the error bound for all probabilities
* @param pErrorBounds the error bounds for all probabilities
* NOTE: Either we have pErrorBounds != NULL, or the error bound is defined by error_bound
*/
extern void print_state_probs(int size, const double * probs, const char * name,
                        const double error_bound,
			const double * pErrorBounds );

/**
* This function prints the probability for the state with the given index.
* Possible lumping is taken into account.
* @param index: the original-state-space index (even if lumping was done)
* @param size: the size of 'probs' array
* @param probs: the array of probabilities
* @param name the bitset name string
* @param error_bound the error bound for the given state probability
* @param pErrorBounds the error bounds for all probabilities
* NOTE: Either we have pErrorBounds != NULL, or the error bound is defined by error_bound
*/
extern
inline void print_state_prob( const int index, const int size, const double * probs, const char * name,
				const double error_bound, const double * pErrorBounds );

/**
* This function prints the states, that have confidence intervals not tighter
* than the specified indiff_width.
*
* NOTE: Just in case we also check if the conf. int. is invalid i.e. ( ciRightBorder - ciLeftBorder ) < 0
*
* @param result_size the size of the probabilities vector
* @param probCILeft the array of left confidence interval borders
* @param probCIRight the array of right confidence interval borders
* @param indiff_width the specified width a confidence interval may not exceed
* @param bitset_name the bitset name string for the error states
*/
extern void print_indiff_err_states( const int result_size, const double *  probCILeft, const double * probCIRight,
                                        const double indiff_width,
                                        const BOOL isSimOneInitState_local,
					const int initial_state, const char * bitset_name);

#endif
