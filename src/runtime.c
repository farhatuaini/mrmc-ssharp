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
*	Source description: Store global variable: state_space, labelling,
*				result_bitset.
*	Uses: DEF: sparse.h, bitset.h, label.h, runtime.h
*		LIB: sparse.c, bitset.c, label.c, runtime.c
*		Definition of runtime - runtime.h
*/

#include "runtime.h"

#include "simulation_utils.h"
#include "simulation.h"
#include "rand_num_generator.h"

#include <float.h>
#include <math.h>

/**
* This boolean variable indicates the formula independent lumping is done already.
* It should allow to avoid cases when we have hested call of operators each of
* which wants to do lumping.
*/
static BOOL formula_lumping_is_done = FALSE;

/**
* This will be used as a bitset of 32 bits to hold the input
* and runtime flag settings, such as: CTMC_MODE, ... , F_DEP_LUMPING_MODE,
* TEST_VMV_MODE or combinations of them.
*/
static unsigned int mode = BLANK_MODE;

/**
* True if the steady-state detection is on
*/
static BOOL ssd_on = FALSE;

/**
* Holds the curent statespace, it is changed if lumping is used.
*/
static sparse *state_space = NULL;

/**
* The partitioning for formula independent lumping
*/
static partition *P = NULL;

/**
* Allowed to be NULL, if not specified, for MRMs (continuous).
*/
static sparse *pImpulseRewards = NULL;

static const labelling * labeller = NULL;

static NDSparseMatrix *mdpi_state_space = NULL;
static int method_ctmpdi_transient = CTMDPI_HD_AUTO_METHOD;

/**
* This variable stores, which logic comparator is given in the formula.
*/
static int comparator;

/**
* This array contains the sum of ALL row elements of the current
* state_space matrix. It means that including the diagonal values
* and this it is not exactly what you have on the diagonal of the
* Generator matrix! To obtain the i'th diagonal value of the
* generator matrix you should do (state_space->val[i].diag - row_sums[i])
*/
static double *row_sums = NULL;

/**
* Rewards for DTRM
*/
static double * pRewards = NULL;
static double error_bound=1e-6;
static int method_path=GS, method_steady=GS;
static int method_bscc=REC;
static int method_until_rewards=DTV;
static int max_iterations=1000000;
static double un=DBL_MIN, ov=DBL_MAX;
static double d_factor = (double)1/32.0;
static double w = 1e-11;

/************************************************************************************/
/******************************THE RUN-MODE ACCESS FUNCTIONS*************************/
/************************************************************************************/

/**
* This method sets the mode in which we run the mcc tool.
* It sets a flag to a bitset which can have the following flags
* on or off:  CTMC_MODE, DTMC_MODE, DMRM_MODE, CMRM_MODE etc.
*/
void addRunMode(unsigned int flag)
{
	set_flag(mode, flag);
}

/**
* This method clears part of the mode in which we run the mcc tool.
* It clears the flag from a bitset which can have the following flags
* on or off:  CTMC_MODE, DTMC_MODE, DMRM_MODE, CMRM_MODE etc.
*/
void clearRunMode(unsigned int flag)
{
	clear_flag(mode, flag);
}

/**
* This method tests the mode in which we runthe mcc tool
* It returns true(The logic ID) if the given flag is set
* in a mode bitset which can have the following flags on
* or off: CTMC_MODE, DTMC_MODE, DMRM_MODE, CMRM_MODE etc.
*
* Example: CTMC_MODE == isRunMode(CTMC_MODE) in case we are
* in CTMC_MODE otherwise BLANK_MODE == isRunMode(CTMC_MODE)
*/
unsigned int isRunMode(unsigned int flag)
{
	return test_flag(mode, flag);
}

/**
* This method using the mode local file variable and the ANY_MODEL_MODE
* macros value tests whether the logic has been already set. Returns
* True if yes, otherwise, False.
*/
unsigned int isRunModeSet(void)
{
	return test_flag(mode, ANY_MODEL_MODE);
}

/************************************************************************************/
/****************************THE FORMULA-TREE ACCESS FUNCTIONS***********************/
/************************************************************************************/

/* This variable stores the globally accessible result of the */
/* lately model-checked formula. It is in a form of a formulatree */
static void* pFormulaTreeRootNode = NULL;

/**
* Sets the formula tree after it has been modelchecked
* @param pFormulaTreeRootNode the formula tree with the results.
*/
void set_formula_tree_result(void * _pFormulaTreeRootNode){
	pFormulaTreeRootNode = _pFormulaTreeRootNode;
}

/**
* Gets the formula tree after it has been modelchecked
* @return the formula tree with the results.
*/
void * get_formula_tree_result(void) {
	return pFormulaTreeRootNode;
}

/************************************************************************************/
/****************************THE STATE-SPACE ACCESS FUNCTIONS************************/
/************************************************************************************/

/**
* This function sets the state space for global access.
* @param space the state space
* WARNING: If space == NULL then the row sums are freed using the free_row_sums() method.
*/
void set_state_space(sparse *space){
	state_space = space;
	/* WARNING: In principle set_state_space should not be called with a NULL parameter */
	/* This has to be done only if we want to reset the matrix*/
	if( state_space != NULL ){
		row_sums = get_mtx_row_sums(space);
                if ( NULL == row_sums ) {
                        exit(err_macro_3(err_CALLBY, "set_state_space(%p[%dx"
                                "%d])", (void *) space, mtx_rows(space),
                                mtx_cols(space), EXIT_FAILURE));
                }
	}else{
		free_row_sums();
	}
}


/**
* Get the currently set state space.
*
* WARNING: Can be used only in DTMC, CTMC, DMRM, and CMRM modes!
*
* @return the pointer to the sparse matrix containing the current state.
*/
sparse * get_state_space(void) {
	return state_space;
}

/**
* Globally access the MDP state space.
* @return the pointer to the sparse matrix of the MDP while im CTMDPI_MODE.
*/
NDSparseMatrix * get_mdpi_state_space(void)
{
	return mdpi_state_space;
}

/**
 * Sets method to be used for CTMDP transient analysis.
 * The method may be one of CTMDPI_UNIFORM_METHOD or
 * CTMDPI_NONUNIFORM_METHOD.
 *
 * @param new_method method to be set
 */
void set_method_ctmdpi_transient(int new_method) {
	method_ctmpdi_transient = new_method;
}

int get_method_ctmdpi_transient(void) {
	return method_ctmpdi_transient;
}

/**
* Get the state space size
*
* WARNING: Works for DTMC, CTMC, DMRM, CMRM and SCTMDPI modes!
*
* @return the number of states in the currently used model, i.e.
*               "mtx_rows(get_state_space())" or "get_mdpi_state_space()->n"
*	depending on the run-time mode.
*/
int get_state_space_size(void) {
	int state_space_size = 0;
	if( isRunMode(CTMDPI_MODE) ){
		IF_SAFETY( mdpi_state_space != NULL )
			state_space_size = mdpi_state_space->n;
		ELSE_SAFETY
                        printf("ERROR: trying to retrieve dimensions of a "
                                "NULL pointed CTMDPI.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	}else{
		IF_SAFETY( state_space != NULL )
                        state_space_size = mtx_rows(state_space);
		ELSE_SAFETY
                        printf("ERROR: trying to retrieve dimensions of a "
                                "NULL pointed DTMC/CTMC/DMRM/CMRM.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	}
	return state_space_size;
}

/**
* Set MDP state space for global access.
* @param mdp pointer to mdpi state space
*/
void set_mdpi_state_space(NDSparseMatrix *mdp)
{
	mdpi_state_space = mdp;
}

/**
* Frees the MDP state space.
*/
void freeNDSparseMatrix(void) {
	if (NULL != mdpi_state_space) {
		NDSparseMatrix_free(mdpi_state_space);
	}
}

/************************************************************************************/
/*****************************THE LABELLING ACCESS FUNCTIONS*************************/
/************************************************************************************/

/**
* Sets labelling function for global access using get_labeller.
* @param labelling the labelling function.
*/
void set_labeller(const labelling * labellin)
{
	labeller = labellin;
}

/**
* Access the labelling function.
* @return returns a pointer to the structure labelling.
*/
const labelling * get_labeller(void)
{
	return labeller;
}

/************************************************************************************/
/***************************THE STATE REWARDS ACCESS FUNCTIONS***********************/
/************************************************************************************/

/**
* Set State Rewards
*/
void setStateRewards(double * _pRewards)
{
	pRewards = _pRewards;
}

/**
* Get Rewards
*/
double * getStateRewards(void)
{
	return pRewards;
}

/**
* Free memory
*/
void freeStateRewards(void)
{
	if( pRewards ) free(pRewards);
	pRewards = NULL;
}

/************************************************************************************/
/**************************IMPULSE REWARDS ACCESS FUNCTIONS**************************/
/************************************************************************************/

/**
* Get the value of impulse_rewards
* @return impulse rewards
*/
const sparse * getImpulseRewards(void)
{
	return pImpulseRewards;
}

/**
* Set the value of pImpulseRewards
* @param sparse* _impulse_rewards
*/
void setImpulseRewards(sparse * _pImpulseRewards)
{
	pImpulseRewards=_pImpulseRewards;
	set_method_until_rewards(UQS);
}

/**
* This method is used to free the impulse reward structure.
*/
void freeImpulseRewards(void) {
	if(pImpulseRewards){
                if ( err_state_iserror(free_sparse_ncolse(pImpulseRewards)) ) {
                        exit(err_macro_0(err_CALLBY, "freeImpulseRewards()",
                                EXIT_FAILURE));
                }
		pImpulseRewards = NULL;
	}
}

/************************************************************************************/
/********************************SIMULATION MODE SETTINGS****************************/
/************************************************************************************/

/**
* This method should be used to set the simulation ON/OFF
* indicator in the runtime settings.
* @param BOOL val TRUE to turn simulation ON.
* Provide FALSE for switching it OFF.
* NOTE: in case the run mode is F_DEP_LUMP_MODE the simulation mode does not work
*/
void setSimulationStatus(BOOL val){
	if( val ){
		if( ! isRunMode(F_DEP_LUMP_MODE) && ! isRunMode(F_IND_LUMP_MODE) &&
			( isRunMode(CTMC_MODE) || isRunMode(CMRM_MODE) ) ){
			/* Update the runtime mode */
			addRunMode( SIMULATION_MODE );

			/* The default RNG methods are set beforehand, but GSL functions */
			/* need an explicit initialization. */
			/* Generating of the new seeds is done in the set methods */
			setRNGMethodDiscrete( getRNGMethodDiscrete() );
			setRNGMethodExp( getRNGMethodExp() );
		}else{
			printf("WARNING: The simulation engine is only available in CTMC (CMRM) mode without lumping.\n");
			printf("WARNING: The 'set' command is ignored.\n");
		}
	}else{
		/* Update the runtime mode */
		clearRunMode( SIMULATION_MODE );

		/* Free the random-number generator data, especially needed by GSL functions */
		freeRNGDiscrete();
		freeRNGExp();
	}
}

/**
* This method should be used to get the simulation ON/OFF
* indicator in the runtime settings.
* @return TRUE if the simulation engine is ON, otherwise FALSE.
*/
BOOL isSimulationOn(void) {
	return (isRunMode( SIMULATION_MODE ) ? TRUE : FALSE);
}

/**
* This method allows to set the initial state for the "one initial state" simulation.
* NOTE: We do not call the setSimInitialState(int) method of simulation.h directly
*	because if case of formula independent lumping the state-space is partitioned.
* @param initial_state the initial state index
*/
inline void setSimInitialStateRuntime( const int user_initial_state ){
	int true_initial_state = user_initial_state - 1;
	if( isRunMode(F_IND_LUMP_MODE) ){
		true_initial_state = getLumpedStateIndex( getPartition(), true_initial_state );
	}
        /* Check if the index is within the allowed range. Note that the number
           of states in the original state space is never smaller than in the */
	/* lumped one. Therefore this check does no harm in any case. */
	if( ( true_initial_state >= 0 ) && ( true_initial_state < get_state_space_size() ) ){
		setSimInitialState( true_initial_state );
	} else {
		printf("ERROR: The state index is outside the allowed state range.\n");
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/************************************************************************************/
/*********************************LUMPING MODE SETTINGS******************************/
/************************************************************************************/

/**
* This method is used for setting the partitioning
* for formula independent lumping
* @param _P the partition for formula independent lumping
*/
void setPartition(partition *_P){
	P = _P;
}

/**
* This method is used for getting the partitioning
* for formula independent lumping
* @return The partition that was set before by the
*	  setPartition(partition *) method or NULL.
*/
const partition * getPartition(void) {
	return P;
}

/**
* This method is used for freeying the memory of partition
* for formula independent lumping
*/
void freePartition(void) {
        if ( isRunMode(F_IND_LUMP_MODE) && NULL != P ) {
		free_partition(P);
		P = NULL;
	}
}

/**
* Lets us know whether we are working with the statespace after it was lumped.
* @return TRUE if we are working with the lumped markov chain
*/
BOOL isFormulaLumpingDone(void) {
	return formula_lumping_is_done;
}

/**
* Is used to set Formula dependent lumping ON/OFF
* @param _formula_lumping_is_done TRUE if we just did the formula dependent lumping
*					  FALSE if we unlumped the state space etc. back.
*/
void setFormulaLumpingDone(BOOL _formula_lumping_is_done){
	if( formula_lumping_is_done && _formula_lumping_is_done ){
		printf("ERROR: Trying to do second formula dependent lumping in a row!\n");
	}
	formula_lumping_is_done = _formula_lumping_is_done;
}

/**
* Here we get the original state-space size,
* because there could be lumping done before.
* @return the original state-space size
*/
int get_original_state_space_size(void) {
	int result = 0;
        const partition * P_local = getPartition();
        if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
                result = get_unlumped_state_space_size(P_local);
	}else{
		result = get_labeller()->ns;
	}
	return result;
}

/************************************************************************************/
/****************************MINOR SETTING ACCESS FUNCTIONS**************************/
/************************************************************************************/

/**
* Gets the value of w
*/
double get_w(void)
{
	return w;
}

/**
* Sets the value of w
* @param double the value of w
*/
void set_w(double _w)
{
	w=_w;
}


/*****************************************************************************
name		: get_d_factor
role		: get the value of d
@param		:
@return         : double: d
******************************************************************************/
double get_d_factor(void)
{
	return d_factor;
}

/*****************************************************************************
name		: set_d_factor
role		: set the value of d_factor
@param		: double: _d_factor
remark		:
******************************************************************************/
void set_d_factor(double _d_factor)
{
	d_factor=_d_factor;
}

/*****************************************************************************
name		: get_max_iterations
role		: get the value of max_iterations
@param		:
@return         : int: max_iterations
******************************************************************************/
int get_max_iterations(void)
{
	return max_iterations;
}

/*****************************************************************************
name		: set_max_iterations
role		: set the value of max_iterations
@param		: int: max_iterations
remark		:
******************************************************************************/
void set_max_iterations(int _max_iterations)
{
	max_iterations=_max_iterations;
}

/*****************************************************************************
name		: get_underflow
role		: get the value of underflow
@param		:
@return         : double: underflow
******************************************************************************/
double get_underflow(void)
{
	return un;
}

/*****************************************************************************
name		: set_underflow
role		: set the value of underflow
@param		: double: underflow
remark		:
******************************************************************************/
void set_underflow(double underflow)
{
	un=underflow;
}

/*****************************************************************************
name		: get_overflow
role		: get overflow
@param		:
@return         : double: overflow
******************************************************************************/
double get_overflow(void)
{
	return ov;
}

/*****************************************************************************
name		: set_overflow
role		: set the value of overflow
@param		: double: overflow
remark		:
******************************************************************************/
void set_overflow(double overflow)
{
	ov=overflow;
}

/*****************************************************************************
name		: set_method_path
role		: set method_path for global access.
@param          : int: method_path
remark		: should be either GJ or GS
******************************************************************************/
void set_method_path(int _method_path)
{
	method_path = _method_path;
}

/*****************************************************************************
name		: get_method_path
role		: get method_path
@param		:
@return         : int: method_path
******************************************************************************/
int get_method_path(void)
{
	return method_path;
}

/*****************************************************************************
name		: set_method_steady
role		: set method_steady for global access.
@param          : int: method_steady
remark		: should be either GJ or GS
******************************************************************************/
void set_method_steady(int _method_steady)
{
	method_steady = _method_steady;
}

/*****************************************************************************
name		: set_method_until_rewards
role		: set the method for the evaluation of time-reward-bounded until.
@param		: int: the method.
remark		:
******************************************************************************/
void set_method_until_rewards(int method_until_rew)
{
	method_until_rewards = method_until_rew;
}

/*****************************************************************************
name		: get_method_until_rewards
role		: globally access the result bitset
@return		: int: the method.
remark		:
******************************************************************************/
int get_method_until_rewards(void)
{
	return method_until_rewards;
}

/*****************************************************************************
name		: get_method_steady
role		: get method_steady
@param		:
@return         : int: method_steady
******************************************************************************/
int get_method_steady(void)
{
	return method_steady;
}

/*****************************************************************************
name		: set_error_bound
role		: set error_bound for global access.
@param		: double: the error bound
remark		:
******************************************************************************/
void set_error_bound(double _error_bound)
{
	error_bound = _error_bound;
}

/*****************************************************************************
name		: get_error_bound
role		: get error_bound.
@param		:
@return         : double: the error bound
******************************************************************************/
double get_error_bound(void)
{
	return error_bound;
}

/**
* Set method for the logic comparator.
* @param comp the comparator to be set
* NOTE: Should be one of {C_GREATER, C_GREATER_EQUAL, C_LESS, C_LESS_EQUAL}
*/
void set_comparator(int comp)
{
	comparator = comp;
}

/**
* Get method for the comparator sign.
*/
int get_comparator(void)
{
	return comparator;
}


/************************************************************************************/
/******************************THE BSCC SEARCH SETTINGS******************************/
/************************************************************************************/

/**
* Set method for the BSCC search
* @param the method to be set
* NOTE: The method should be either REC or NON_REC
*/
void set_method_bscc(int _method_bscc)
{
	method_bscc = _method_bscc;
}

/**
* Get method for the BSCC search
* @param the method to be set
* NOTE: The method should be either REC or NON_REC
*/
int get_method_bscc(void)
{
	return method_bscc;
}

/************************************************************************************/
/***************************THE STEADY-STATE DETECTION SETTINGS**********************/
/************************************************************************************/

/**
* This method is used to set the steady-state detection on and off
* @param _on_off TRUE for setting the steady-state detection on, otherwise FALSE.
*/
void set_ssd(BOOL _on_off)
{
  ssd_on = _on_off;
}

/**
* Enables steady-state detection for uniformization (CSL)
*/
void set_ssd_on(void)
{
  set_ssd(TRUE);
}

/**
* Disables steady-state detection for uniformization (CSL)
*/
void set_ssd_off(void)
{
  set_ssd(FALSE);
}

/**
* @return TRUE if steady-state detection is enabled, otherwise FALSE
*/
BOOL is_ssd_on(void)
{
   return ssd_on;
}

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
const double * get_row_sums(void) {
	return row_sums;
}

/**
* This method frees the row_sums array.
* WARNING: It is called automatically if the set_state_space(...) method is called with the NULL parameter
*/
void free_row_sums(void) {
	if( row_sums != NULL){
		free(row_sums);
		row_sums = NULL;
	}
}

/************************************************************************************/
/*********************************PRINTING MODE SETTINGS*****************************/
/************************************************************************************/

static BOOL printing_status = TRUE;

/**
* This method should be used to set the printing ON/OFF
* indicator in the runtime settings.
* Prints a warning message if the status is set to FALSE
* @param BOOL val TRUE to turn printing of the
*				  probabilities and states ON.
*				  Provide FALSE for switching
*				  it OFF.
*/
void setPrintingStatus(BOOL val){
	printing_status = val;
	if( ! val ){
		printf("WARNING: Printing of results is now mostly disabled, use 'help' for more options.\n");
	}
}

/**
* This method should be used to get the printing ON/OFF
* indicator in the runtime settings.
* @return TRUE if the printing of the probabilities and
*	    states is ON, otherwise FALSE.
*/
BOOL isPrintingOn(void) {
	return printing_status;
}

/************************************************************************************/
/**********************************THE PRINTING FUNCTIONS****************************/
/************************************************************************************/

/**
* This method is used for printing the logic name and also extra options such as Lumping mode.
* @param logic: if set to TRUE then the Logic info is printed
* @param extra: if set to TRUE then the Lumping info is printed
*/
inline void print_run_mode(BOOL logic, BOOL extra){
	if( logic ){
		/* Logic mode */
		if( isRunMode(TEST_VMV_MODE) ){
			printf(" Mode\t\t\t = SELF TEST\n");
		}else{
			printf(" Logic\t\t\t = ");
			if( isRunMode(CTMC_MODE) ){
				printf("CSL\n");
			}else{
				if( isRunMode(DTMC_MODE) ){
					printf("PCTL\n");
				}else{
					if( isRunMode(DMRM_MODE) ){
						printf("PRCTL\n");
					}else{
						if( isRunMode(CMRM_MODE) ){
							printf("CSRL\n");
						}else{
							if( isRunMode(CTMDPI_MODE) ){
								printf("CSL\n");
							}else{
								printf("UNKNOWN!\n");
							}
						}
					}
				}
			}
		}
	}
	if( extra ){
		printf(" Formula ind. lumping\t = %s\n",(isRunMode(F_IND_LUMP_MODE) ? "ON":"OFF"));
		printf(" Formula dep. lumping\t = %s\n",(isRunMode(F_DEP_LUMP_MODE) ? "ON":"OFF"));
	}
}

/**
* This method is used to print the runtime settings of error
* bounds, maximum number of iterations etc.
*/
void print_runtime_info(void) {
	int mp  = get_method_path();
	int ms  = get_method_steady();
	int mur = get_method_until_rewards();
	int mb  = get_method_bscc();

	printf(" ---General settings:\n");
	print_run_mode( TRUE, FALSE );
	print_run_mode( FALSE, TRUE );
	printf(" M. C. simulation\t = %s\n", (isSimulationOn() ? "ON":"OFF" ) );
	if( isRunMode(CTMC_MODE) || isRunMode(CMRM_MODE) ){
		printf(" Steady-state detection\t = %s\n", (is_ssd_on() ? "ON":"OFF"));
	}
	printf(" Method Path\t\t = ");
	switch(mp){
		case GJ:
			printf("Gauss-Jacobi\n");
			break;
		case GS:
			printf("Gauss-Seidel\n");
			break;
                default:
                        fprintf(stderr,
                                "print_runtime_info: illegal Method Path\n");
                        exit(EXIT_FAILURE);
	}
	printf(" Method Steady\t\t = ");
	switch(ms){
		case GJ:
			printf("Gauss-Jacobi\n");
			break;
		case GS:
			printf("Gauss-Seidel\n");
			break;
                default:
                        fprintf(stderr,
                                "print_runtime_info: illegal Method Steady\n");
                        exit(EXIT_FAILURE);
	}
	if( isRunMode(CMRM_MODE) ){
		printf(" Method Until Rewards\t = ");
		switch(mur){
			case US:
				printf("Uniformization-Sericola\n");
				break;
			case UQS:
				printf("Uniformization-Qureshi-Sanders\n");
				break;
			case DTV:
				printf("Discretization-Tijms-Veldman\n");
				break;
                        default:
                                fprintf(stderr, "print_runtime_info: "
                                        "illegal Method Until Rewards\n");
                                exit(EXIT_FAILURE);
		}
	}
	printf(" Method BSCC\t\t = ");
	switch(mb){
		case REC:
			printf("Recursive\n");
			break;
		case NON_REC:
			printf("Non-Recursive\n");
			break;
                default:
                        fprintf(stderr,
                                "print_runtime_info: illegal Method BSCC\n");
                        exit(EXIT_FAILURE);
	}
	printf(" Results printing\t = %s\n", (isPrintingOn()? "ON":"OFF") );
	printf("\n");

	/* Print the simulation runtime parameters */
	printRuntimeSimInfo( isRunMode( CTMC_MODE ) || isRunMode( CMRM_MODE ) );

	printf(" ---Numerical methods:\n");
	printf(" -Iterative solvers:\n");
	printf("   Error Bound\t\t\t = %e\n", get_error_bound());
	printf("   Max Iterations\t\t = %ld\n", (long) get_max_iterations());
	if( isRunMode(CTMC_MODE) || isRunMode(CMRM_MODE) ){
		printf(" -Fox-Glynn algorithm:\n");
		printf("   Overflow\t\t\t = %e\n", get_overflow());
		printf("   Underflow\t\t\t = %e\n", get_underflow());
	}
	if( isRunMode(CMRM_MODE) ){
		printf(" -Uniformization Qureshi-Sanders:\n");
		printf("   Probability threshold\t = %e\n", get_w());
		printf(" -Discretization Tijms-Veldman:\n");
		printf("   Discretization factor\t = %e\n", get_d_factor());
	}
}

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
int get_error_bound_precision(const int index, const double error_bound_local,
                const int size, const double * pErrorBounds)
{
	double temp_error_bound;
	int i;

	if( pErrorBounds == NULL ){
		/* There are no separate error bounds for states, thus use the global one */
                temp_error_bound = error_bound_local;
	} else {
		if( index >= 0){
			/* There is a state index provided, and pErrorBounds != NULL */
			IF_SAFETY( index < size )
				/* Take the error of the given state */
				temp_error_bound = pErrorBounds[index];
			ELSE_SAFETY
				printf("ERROR: Trying to get an error bound for the state index %d, whereas we only have %d states.\n", index+1, size);
                                exit(EXIT_FAILURE);
			ENDIF_SAFETY
		}else{
			IF_SAFETY( size > 0 )
				/* We are to print an array of values, for each of which we have */
				/* a specific error. In this case compute the largest error and */
				/* use it for generating the precision pattern */
				temp_error_bound = pErrorBounds[0];
				for( i = 1; i < size; i++){
					if( temp_error_bound < pErrorBounds[i] ) {
						temp_error_bound = pErrorBounds[i];
					}
				}
			ELSE_SAFETY
				printf("ERROR: Trying to get error bounds but the error-array size is %d.\n", size);
                                exit(EXIT_FAILURE);
			ENDIF_SAFETY
		}
	}

	/* If the error bound computed so far seems to be zero, then return */
        /* the pattern for the default one, i.e. given by get_error_bound() */
	return get_precision( ( temp_error_bound == 0.0 ) ? get_error_bound() : temp_error_bound ) ;
}

/**
* This method is used for getting the printing pattern from the value.
* NOTE: If we have a value 1.231e-10 then the precision
* is 10+1 digits after the decimal point, but not 13!
* @param value: the value to evaluate
* @return the precision for the printing pattern by taking the position
*	    of the first significant digit after the decimal point + 1
*/
inline int get_precision(double value){
	int precision = 0;
        double integer = 0.0, fraction = 0.0;

        value = fabs(value);
        /* NOTE: If we have error 1.1e-10 then the precision */
	/* is 10 digits after the decimal point, but not 11. */
        do {
		value *= 10;
		fraction = modf(value, &integer);
		precision++;
        } while ( integer <= 0.0 );

	return precision + 1;
}

/**
* This method is used for printing a bitset which is an
* indicator set of states that satisfy the property
* @param b the bitset which has to be printed
* @param name the bitset name string
*/
void print_states(const bitset * b, const char * name) {
        const partition * P_local = getPartition();

	printf("%s: ", name);
        if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
                print_original_states(P_local, b);
	}else{
		IF_SAFETY( b != NULL )
		print_bitset_states(b);
		ELSE_SAFETY
			printf("??\n");
			printf("ERROR: Trying to print a non-existing bitset.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	}
	printf("\n");
}

/**
* This method is used for printing TRUE/FALSE result for a state index.
* TRUE if the state 'index' is in the set 'pBitset', otherwise FALSE.
* This method handles lumping.
* @param pBitset: the bitset which has to be printed
* @param index: the state index in the original state space.
* @param name the bitset name string
*/
void print_state( const bitset * pBitset, const int index, const char * name ){
        const partition * P_local = getPartition();

	printf("%s[%d] = ", name, index+1);
        if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
                print_original_state(P_local, pBitset, index);
	}else{
		if( pBitset != NULL ){
                        if ( 0 <= index && index < bitset_size(pBitset) ) {
				if( get_bit_val( pBitset, index ) ){
					printf("TRUE");
				}else{
					printf("FALSE");
				}
			}else{
				printf("??\n");
                                printf("WARNING: Invalid index %d, required to "
                                       "be in the [1, %d] interval.",
                                       index + 1, bitset_size(pBitset));
			}
		}else{
			printf("??\n");
			printf("WARNING: Trying to print an element of a non-existing bitset.");
		}
	}
	printf("\n");
}

/**
* This method is used for printing the errors of state probabilities
* Depending whether there was lumping done or not to the
* state space it does different printing, Note that, since we are to
* print errors, we do not use any template printing
* @param size the size of the probability-error vector
* @param pErrorValues the vector of probability errors
* @param pName the name string
*/
void print_error_probs( const int size, const double * pErrorValues, const char * pName ){
        const partition * P_local = getPartition();

	printf("%s: ", pName);
        if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
                print_error_probs_partition(P_local, pErrorValues);
	}else{
                if ( err_state_iserror(print_pattern_vec_double("%1.2e", size,
                                pErrorValues)) )
                        exit(err_macro_3(err_CALLBY, "print_error_probs(%d,%p,"
                                "\"%s\")", size, (const void *) pErrorValues,
                                pName, EXIT_FAILURE));
	}
	printf("\n");
}

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
void print_state_probs(int size, const double * probs, const char * name,
                        const double error_bound_local,
			const double * pErrorBounds ){
        const partition * P_local = getPartition();

	printf("%s: ", name);
        if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
                print_state_probs_partition(P_local, probs, error_bound_local,
                                pErrorBounds);
	}else{
		/* Calculate the right pattern according to the precision */
		char buffer[255];
                sprintf(buffer, "%%1.%df", get_error_bound_precision(-1,
                                error_bound_local, size, pErrorBounds));
                if ( err_state_iserror(print_pattern_vec_double(buffer, size,
                                                probs)) )
                        exit(err_macro_5(err_CALLBY, "print_state_probs(%d,%p,"
                                "\"%s\",%g,%p)", size, (const void *) probs,
                                name, error_bound_local,
                                (const void *) pErrorBounds, EXIT_FAILURE));
	}
	printf("\n");
}

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
inline void print_state_prob( const int index, const int size, const double * probs, const char * name,
                                const double error_bound_local,
                                const double * pErrorBounds)
{
        const partition * P_local = getPartition();

	printf( "%s[%d] = ", name, index + 1 );
        if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
                print_state_prob_partition(P_local, index, probs,
                                error_bound_local, pErrorBounds);
	}else{
		if( probs != NULL ){
			if( 0 <= index && index < size ){
				/* Calculate the right pattern according to the precision */
                                printf("%1.*f", get_error_bound_precision(
                                                index, error_bound_local,
                                                size, pErrorBounds ),
                                        probs[index]);
			}else{
				printf("??\n");
				printf("WARNING: Invalid index %d, required to be in the [1, %d] interval.", index + 1, size);
			}
		}else{
			printf("??\n");
			printf("WARNING: Trying to print an element of a non-existing array.");
		}
	}
	printf("\n");
}

/****************************************************************************/
/**********************PRINT THE ERROR SIMULATION STATES*********************/
/****************************************************************************/

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
void print_indiff_err_states( const int result_size, const double *  probCILeft, const double * probCIRight,
                                const double indiff_width,
                                const BOOL isSimOneInitState_local,
				const int initial_state, const char * bitset_name){
	int i;
	/* Indicated whether there are states found for which the conf. int. is invalid or not tight enough*/
	BOOL isIndiffErrStates = FALSE;
	/*This bitset will be allocated later and will store invalid states*/
	bitset * pErrStatesBitset = NULL;

	IF_SAFETY( ( probCILeft != NULL ) && ( probCIRight != NULL ) )
	    /*Check if there are bad conf. int., i.e. an invalid one or the one with the width >= indiff_width*/
            if ( isSimOneInitState_local ) {
		    /*To make things uniform we create an entire bitset which can store the initial state index*/
		    pErrStatesBitset = get_new_bitset(initial_state + 1);
		    /*In case of one initial state there should be only one element in each of the conf. int. arrays*/
		    checkErrorConfIntState(probCILeft[0], probCIRight[0], initial_state, indiff_width, &isIndiffErrStates, pErrStatesBitset);
	    } else {
		    /*Allocate a bitset for all considered initial states*/
		    pErrStatesBitset = get_new_bitset(result_size);
		    /*Check  the conf. int. for all initial states*/
		    for(i = 0; i < result_size; i++){
			    checkErrorConfIntState(probCILeft[i], probCIRight[i], i, indiff_width, &isIndiffErrStates, pErrStatesBitset);
		    }
	    }
	ELSE_SAFETY
		printf("ERROR: At least one of the confidence-interval bound arrays is NULL.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY

	/*If there are bad confidence intervals, then we report on them*/
	if( isIndiffErrStates ){
		print_states(pErrStatesBitset, bitset_name);
		printf("WARNING: Increase max_sample_size for obtaining the conf. int. of the desired width.\n");
	}

	free_bitset(pErrStatesBitset);
}
