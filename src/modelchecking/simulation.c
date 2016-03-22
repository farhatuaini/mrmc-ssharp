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
*	Source description: This is a header file for the simulation engine
*	here we intend to define the common functions and data structures.
*/

#include "simulation.h"

#include "rand_num_generator.h"

#include "runtime.h"

/****************************************************************************/
/********************MANAGE GENERAL SIMULATION PARAMETERS********************/
/****************************************************************************/

/* Defines the mode for model checking the steady-state operator. */
/* Should be FALSE for the Hybrid mode, otherwise TRUE */
static BOOL isSteadyStatePure = FALSE;

/**
* This method should be used to turn on the pure simulation
* for the steady-state (long-run) operator
* @param isOnOff TRUE to turn on the pure simulation mode, FALSE for the hybrid mode
*/
inline void setSimSteadyStateModePure( BOOL isOnOff ){
	isSteadyStatePure = isOnOff;
}

/**
* This method should be used to get the pure-simulation ON/OFF
* (for the steady-state (long-run) operator) indicator in the runtime settings.
* @return TRUE if the pure-simulation is ON, otherwise FALSE.
*/
BOOL isSimSteadyStateModePure(void) {
	return isSteadyStatePure;
}

/**
* This method should be used to get the hybrid-simulation ON/OFF
* (for the steady-state (long-run) operator) indicator in the runtime settings.
* @return TRUE if the hybrid-simulation is ON, otherwise FALSE.
*/
BOOL isSimSteadyStateModeHybrid(void) {
	return ! isSimSteadyStateModePure();
}

/* Defines the regeneration state initializing for the regenerative */
/* method used when model checking the steady-state operator. */
static BOOL isRegMethodPure = FALSE;

/**
* This method should be used to turn on the pure regenerative method
* for the steady-state operator.
* @param ifOnOff TRUE to turn the pure regenerative method on, FALSE for the
*			regenerative method with heuristic regeneration state
*			initialization
*/
inline void setRegMethodModePure( BOOL isOnOff ){
	isRegMethodPure = isOnOff;
}

/**
* This method should be used to get the pure regenerative method ON/OFF indicator
* for the steady-state operator.
* @return TRUE if the pure regenerative method is ON, otherwise FALSE.
*/
BOOL isRegMethodModePure(void) {
	return isRegMethodPure;
}

/**
* This method should be used to get the heuristic regenerative method ON/OFF indicator
* for the steady-state operator.
* @return TRUE if the heuristic regenerative method is ON, otherwise FALSE.
*/
BOOL isRegMethodModeHeuristic(void) {
	return !isRegMethodModePure();
}

/* The initial state for the "one initial state" simulation */
static int sim_initial_state = 0;

/**
* Sets the initial state that is used for the "one initial state" simulation.
* @param the initial state for doing simulations
*/
inline void setSimInitialState( int initial_state ){
	sim_initial_state = initial_state;
}

/**
* Gets the initial state that is used for the "one initial state" simulation.
* @return the initial state for doing simulations
*/
int getSimInitialState(void) {
	return sim_initial_state;
}

/* Defines if one or all states should be considered as initial */
/* TRUE for one initial state, otherwise FALSE; */
static BOOL isOneInitState = FALSE;

/**
* This method sets the type of simulation we have in mind,
* for one initial state only or for all stetes
* @param isOnOff TRUE for just one initial state, FALSE for all states
*/
void setSimOneInitState(BOOL isOnOff){
	isOneInitState = isOnOff;
}

/**
* This method indicates the type of simulation we have in mind,
* for one initial state only or for all stetes.
* @return TRUE for just one initial state, FALSE for all states
*/
BOOL isSimOneInitState(void) {
	return isOneInitState;
}

/* The dewised conf. int. width should be smaller than indiff_width */
static double indiff_width = DEF_INDIFFERENCE_REGION;

/**
 * This method sets the desired supremum of the indefference-region width.
 * @param _indiff_width the dewised conf. int. width should be smaller than indiff_width
 */
inline void setSupIndifferenceWidth( const double _indiff_width ){
	if( ( _indiff_width > MIN_INDIFFERENCE_REGION ) && ( _indiff_width <= MAX_INDIFFERENCE_REGION ) ){
		indiff_width = _indiff_width;
	}else{
		printf("WARNING: The width of the indifference region should be > %e and <= %e\n",
				MIN_INDIFFERENCE_REGION, MAX_INDIFFERENCE_REGION);
		printf("WARNING: The set command is ignored.\n");
	}
}

/**
 * This method gets the desired supremum of the indefference-region width.
 * @return the dewised conf. int. width should be smaller than this value
 */
double getSupIndifferenceWidth(void) {
	return indiff_width;
}

/* Stores the global confidence level for model */
/* checking the formula with simulation */
/* WARNING: Should have DEF_GENERAL_CONFIDENCE >= MIN_GENERAL_CONFIDENCE */
static double confidence = DEF_GENERAL_CONFIDENCE;

/**
* Sets the general confidence level in case of simulations.
* @param _confidence the confidence level for the formula
*/
void setSimGeneralConfidence(double _confidence){
	if( ( _confidence >= MIN_GENERAL_CONFIDENCE ) && ( _confidence < MAX_GENERAL_CONFIDENCE ) ){
		confidence = _confidence;
	}else{
		printf("WARNING: The confidence level should be >= %e and < %e.\n",
				MIN_GENERAL_CONFIDENCE, MAX_GENERAL_CONFIDENCE);
		printf("WARNING: The set command is ignored.\n");
	}
}

/**
* Gets the general confidence level in case of simulations.
* @return the confidence level for the formula
*/
double getSimGeneralConfidence(void) {
	return confidence;
}

/****************************************************************************/
/**********************MANAGE THE SAMPLE PARAMETERS**************************/
/****************************************************************************/

/* Stores the maximum sample size */
static int max_sample_size = DEF_MAX_SAMPLE_SIZE;

/* Stores the minimum sample size */
static int min_sample_size = DEF_MIN_SAMPLE_SIZE;

/* Stores the minimum sample size */
static int sample_size_step = DEF_SAMPLE_SIZE_STEP;

/**
* Sets the maximum sample size
* @param _max_sample_size the maximum sample size
*/
inline void setSimMaxSampleSize( int _max_sample_size ){
	if( _max_sample_size >= min_sample_size ){
		max_sample_size = _max_sample_size;
	} else {
		printf("WARNING: The maximum sample size should be >= %d.\n", min_sample_size );
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/**
* Gets the maximum sample size
* @return the maximum sample size
*/
int getSimMaxSampleSize(void) {
	return max_sample_size;
}

/**
* Sets the minimum sample size
* @param _min_sample_size the minimum sample size
*/
inline void setSimMinSampleSize( int _min_sample_size ){
	if( _min_sample_size >= MIN_MIN_SAMPLE_SIZE && _min_sample_size <= max_sample_size ){
		min_sample_size = _min_sample_size;
	} else {
		printf("WARNING: The minimum sample size should be >= %d and <= %d.\n", MIN_MIN_SAMPLE_SIZE, max_sample_size);
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/**
* Gets the minimum sample size
* @return the minimum sample size
*/
int getSimMinSampleSize(void) {
	return min_sample_size;
}

/**
* Sets the sample size step
* @param _sample_size_step the sample size step
*/
inline void setSimSampleSizeStep( int _sample_size_step ){
	if( _sample_size_step >= MIN_SAMPLE_SIZE_STEP ){
		sample_size_step = _sample_size_step;
	} else {
		printf("WARNING: The sample-size step should be >= %d.\n", MIN_SAMPLE_SIZE_STEP );
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/**
* Gets the sample size step
* @return the sample size step
*/
int getSimSampleSizeStep(void) {
	return sample_size_step;
}

/* Defines if the sample size step should be computed automatically*/
/* or set manually */
/* TRUE for automatically, otherwise FALSE; */
static BOOL isSampleSizeStepAuto = TRUE;

/**
* This method sets the type of sample size step generation
* we have in mind, automatically or manually with fixed step size.
* @param isOnOff TRUE for automatically, FALSE for manually
*/
void setSimSampleSizeStepAuto(BOOL isOnOff){
	isSampleSizeStepAuto = isOnOff;
}

/**
* This method indicates the type of sample size step generation
* we have in mind, automatically or manually with fixed step size.
* @return TRUE for automatically, FALSE for manually
*/
BOOL isSimSampleSizeStepAuto(void) {
	return isSampleSizeStepAuto;
}

/****************************************************************************/
/*******************MANAGE THE SIMULATION DEPTH PARAMETERS*******************/
/*******************IS NEEDED FOR M.C. UNBOUNDED UNTL************************/
/***********************USING PURE SIMULATION APPROACH***********************/
/****************************************************************************/

/* Stores the maximum simulation depth */
static int max_simulation_depth = DEF_MAX_SIMULATION_DEPTH;

/* Stores the minimum simulation depth */
static int min_simulation_depth = DEF_MIN_SIMULATION_DEPTH;

/* Stores the minimum simulation depth */
static int simulation_depth_step = DEF_SIMULATION_DEPTH_STEP;

/**
* Sets the maximum simulation depth
* @param _max_simulation_depth the maximum simulation depth
*/
inline void setSimMaxSimulationDepth( int _max_simulation_depth ){
	if( _max_simulation_depth >= min_simulation_depth ){
		max_simulation_depth = _max_simulation_depth;
	} else {
		printf("WARNING: The maximum simulation depth should be >= %d.\n", min_simulation_depth );
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/**
* Gets the maximum simulation depth
* @return the maximum simulation depth
*/
int getSimMaxSimulationDepth(void) {
	return max_simulation_depth;
}

/**
* Sets the minimum simulation depth
* @param _min_simulation_depth the minimum simulation depth
*/
inline void setSimMinSimulationDepth( int _min_simulation_depth ){
	if( _min_simulation_depth >= MIN_MIN_SIMULATION_DEPTH && _min_simulation_depth <= max_simulation_depth ){
		min_simulation_depth = _min_simulation_depth;
	} else {
		printf("WARNING: The minimum simulation depth should be >= %d and <= %d.\n", MIN_MIN_SIMULATION_DEPTH, max_simulation_depth );
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/**
* Gets the minimum simulation depth
* @return the minimum simulation depth
*/
int getSimMinSimulationDepth(void) {
	return min_simulation_depth;
}

/**
* Sets the simulation depth step
* @param _simulation_depth_step the simulation depth step
*/
inline void setSimSimulationDepthStep( int _simulation_depth_step ){
	if( _simulation_depth_step >= MIN_SIMULATION_DEPTH_STEP ){
		simulation_depth_step = _simulation_depth_step;
	} else {
		printf("WARNING: The simulation-depth step should be >= %d.\n", MIN_SIMULATION_DEPTH_STEP );
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/**
* Gets the simulation depth step
* @return the simulation depth step
*/
int getSimSimulationDepthStep(void) {
	return simulation_depth_step;
}

/****************************************************************************/
/*******************MANAGE THE BSCC MULTIPLIER PARAMETER*********************/
/*******************IS NEEDED FOR M.C. STEADY STATE**************************/
/*******************USING HYBRID SIMULATION APPROACH*************************/
/****************************************************************************/

/* Stores the multiplier for choosing the sample-based regeneration point */
static int sampleBSCCDimensionMultiplier = DEF_SAMPLE_BSCC_MULTIPLIER;

/**
* Sets the BSCC dimension multiplier for choosing the sample-based regeneration point.
* @param _sampleBSCCDimensionMultiplier the BSCC dimension multiplier
*/
inline void setSampleBSCCDimensionMultiplier(int _sampleBSCCDimensionMultiplier){
	if( (_sampleBSCCDimensionMultiplier >= MIN_MULTIPLIER_SIZE) && (_sampleBSCCDimensionMultiplier <= MAX_MULTIPLIER_SIZE) ){
		sampleBSCCDimensionMultiplier = _sampleBSCCDimensionMultiplier;
	} else {
		printf("WARNING: The multiplier should be >= %d and <= %d.\n", MIN_MULTIPLIER_SIZE, MAX_MULTIPLIER_SIZE );
		printf("WARNING: The 'set' command is ignored.\n");
	}
}

/**
* Gets the BSCC dimension multiplier for choosing the sample-based regeneration point.
* @return the BSCC dimension multiplier
*/
int getSampleBSCCDimensionMultiplier(void) {
	return sampleBSCCDimensionMultiplier;
}


/****************************************************************************/
/*******************PRINT THE SIMULATION RUNTIME PARAMETERS******************/
/****************************************************************************/

/**
* Print the simulation runtime parameters if the simulation is on.
* Prints nothing otherwise. Also: does not print the RNG for exp.
* if isExpNeeded == FALSE.
* @param isExpNeeded if it is needed to print the info about the exp distrib.
*			random number generator.
*/
inline void printRuntimeSimInfo(BOOL isExpNeeded){
	if( isSimulationOn() ){
		printf(" ---Monte Carlo simulation:\n");
		printf(" Simulation type\t = %s\n", (isSimOneInitState()? "ONE":"ALL") );
		if( isSimOneInitState() ){
			printf(" Sim. initial state\t = %d\n", getSimInitialState()+1 );
		}
		printf(" Sim. steady state\t = %s\n", (isSimSteadyStateModePure() ? "PURE":"HYBRID" ) );
		printf(" Reg. method steady\t = %s\n", (isRegMethodModePure() ? "PURE":"HEURISTIC" ) );
		printf(" Confidence level\t = %e\n", getSimGeneralConfidence() );
		printf(" Indiff. reg. width\t = %e\n", getSupIndifferenceWidth() );
		printf(" Max sample size\t = %d\n", getSimMaxSampleSize() );
		printf(" Min sample size\t = %d\n", getSimMinSampleSize() );
		printf(" Sample-size step type\t = %s\n", (isSimSampleSizeStepAuto() ? "AUTO":"MANUAL" ) );
		printf(" Sample-size step\t = %d\n", getSimSampleSizeStep() );

		/* Print the runtime parameters of the random number generator */
		printRuntimeRNGInfoDiscrete();
		if( isExpNeeded ){
			printRuntimeRNGInfoExp();
		}

		/* Print the options for the unbounded-until simulation */
		/* The should be also used in the steady-state operator simulation */
		/* Then the pure simulation mode is on. */
		printf(" Max simulation depth\t = %d\n", getSimMaxSimulationDepth() );
		printf(" Min simulation depth\t = %d\n", getSimMinSimulationDepth() );
		printf(" Simulation-depth step\t = %d\n", getSimSimulationDepthStep() );

		/* Print the options for the steady-state simulation */
		printf(" BSCC dim. multiplier\t = %d\n", getSampleBSCCDimensionMultiplier() );

		printf("\n");
	}
}
