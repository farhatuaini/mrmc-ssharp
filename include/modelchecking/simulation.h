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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "macro.h"

/****************************************************************************/
/********************MANAGE GENERAL SIMULATION PARAMETERS********************/
/****************************************************************************/

/* The default value of the the dewised conf. int. width */
#define DEF_INDIFFERENCE_REGION 0.02

/* The excluded minimum value for the indifference region */
#define MIN_INDIFFERENCE_REGION 0.0
#define MAX_INDIFFERENCE_REGION 1.0

/* The min/max general confidence allowed */
/* See the theory in the Thesis of Ivan S. Zapreev */
#define MIN_GENERAL_CONFIDENCE 0.25
#define MAX_GENERAL_CONFIDENCE 1.0

/* The default confidence level for model checking any formula */
#define DEF_GENERAL_CONFIDENCE 0.95

/****************************************************************************/
/********************MANAGE GENERAL SIMULATION PARAMETERS********************/
/****************************************************************************/

/**
* This method should be used to turn on the pure simulation
* for the steady-state (long-run) operator
* @return TRUE if the pure-simulation is ON, otherwise FALSE.
*/
extern inline void setSimSteadyStateModePure( BOOL isOnOff );

/**
* This method should be used to get the pure-simulation ON/OFF
* (for the steady-state (long-run) operator) indicator in the runtime settings.
* @return TRUE if the pure-simulation is ON, otherwise FALSE.
*/
extern BOOL isSimSteadyStateModePure(void);

/**
* This method should be used to get the hybrid-simulation ON/OFF
* (for the steady-state (long-run) operator) indicator in the runtime settings.
* @return TRUE if the hybrid-simulation is ON, otherwise FALSE.
*/
extern BOOL isSimSteadyStateModeHybrid(void);

/**
* This method should be used to turn on the pure regenerative method
* for the steady-state operator.
* @param ifOnOff TRUE to turn the pure regenerative method on, FALSE for the
*			regenerative method with heuristic regeneration state
*			initialization
*/
extern inline void setRegMethodModePure( BOOL isONOFF );

/**
* This method should be used to get the pure regenerative method ON/OFF indicator
* for the steady-state operator.
* @return TRUE if the pure regenerative method is ON, otherwise FALSE.
*/
extern BOOL isRegMethodModePure(void);

/**
* This method should be used to get the heuristic regenerative method ON/OFF indicator
* for the steady-state operator.
* @return TRUE if the heuristic regenerative method is ON, otherwise FALSE.
*/
extern BOOL isRegMethodModeHeuristic(void);

/**
* Sets the initial state that is used for the "one initial state" simulation.
* @param the initial state for doing simulations
*/
extern inline void setSimInitialState( int initial_state );

/**
* Gets the initial state that is used for the "one initial state" simulation.
* @return the initial state for doing simulations
*/
extern int getSimInitialState(void);

/**
* This method sets the type of simulation we have in mind,
* for one initial state only or for all stetes
* @param isOnOff TRUE for just one initial state, FALSE for all states
*/
extern inline void setSimOneInitState(BOOL isOnOff);

/**
* This method indicates the type of simulation we have in mind,
* for one initial state only or for all stetes
* @return TRUE for just one initial state, FALSE for all states
*/
extern BOOL isSimOneInitState(void);

/**
 * This method sets the desired supremum of the indefference-region width.
 * @param _indiff_width the dewised conf. int. width should be smaller than sup_width
 */
extern inline void setSupIndifferenceWidth( const double _indiff_width );

/**
 * This method gets the desired supremum of the indefference-region width.
 * @return the dewised conf. int. width should be smaller than this value
 */
extern double getSupIndifferenceWidth(void);

/**
* Sets the general confidence level in case of simulations.
* @param _confidence the confidence level for the formula
*/
extern inline void setSimGeneralConfidence(double _confidence);

/**
* Gets the general confidence level in case of simulations.
* @return the confidence level for the formula
*/
extern double getSimGeneralConfidence(void);

/****************************************************************************/
/**********************MANAGE THE SAMPLE PARAMETERS**************************/
/****************************************************************************/

/* This is the default maximum sample size */
#define DEF_MAX_SAMPLE_SIZE 100000
/* This is the smallest possible minimum sample size */
/* WARNING: Should be at least >=3 ! */
#define MIN_MIN_SAMPLE_SIZE 10
/* WARNING: Should be DEF_MAX_SAMPLE_SIZE >= DEF_MIN_SAMPLE_SIZE >= MIN_MIN_SAMPLE_SIZE */
#define DEF_MIN_SAMPLE_SIZE 10000
/* This is the smallest possible sample size step */
#define MIN_SAMPLE_SIZE_STEP 1
/* WARNING: Should be DEF_SAMPLE_SIZE_STEP >= MIN_SAMPLE_SIZE_STEP */
#define DEF_SAMPLE_SIZE_STEP 100

/**
* Sets the maximum sample size
* @param _max_sample_size the maximum sample size
*/
extern inline void setSimMaxSampleSize( int _max_sample_size );

/**
* Gets the maximum sample size
* @return the maximum sample size
*/
extern int getSimMaxSampleSize(void);

/**
* Sets the minimum sample size
* @param _min_sample_size the minimum sample size
*/
extern inline void setSimMinSampleSize( int _min_sample_size );

/**
* Gets the minimum sample size
* @return the minimum sample size
*/
extern int getSimMinSampleSize(void);

/**
* Sets the sample size step
* @param _sample_size_step the sample size step
*/
extern inline void setSimSampleSizeStep( int _sample_size_step );

/**
* Gets the sample size step
* @return the sample size step
*/
extern int getSimSampleSizeStep(void);

/**
* This method sets the type of sample size step generation
* we have in mind, automatically or manually with fixed step size.
* @param isOnOff TRUE for automatically, FALSE for manually
*/
extern inline void setSimSampleSizeStepAuto(BOOL isOnOff);

/**
* This method indicates the type of sample size step generation
* we have in mind, automatically or manually with fixed step size.
* @return TRUE for automatically, FALSE for manually
*/
extern BOOL isSimSampleSizeStepAuto(void);

/****************************************************************************/
/*******************MANAGE THE SIMULATION DEPTH PARAMETERS*******************/
/*******************IS NEEDED FOR M.C. UNBOUNDED UNTL************************/
/***********************USING PURE SIMULATION APPROACH***********************/
/****************************************************************************/

/* This is the default maximum simulation depth */
#define DEF_MAX_SIMULATION_DEPTH 100000
/* This is the smallest possible minimum simulation depth */
#define MIN_MIN_SIMULATION_DEPTH 1
/* WARNING: Should be DEF_MAX_SIMULATION_DEPTH >= DEF_MIN_SIMULATION_DEPTH >= MIN_MIN_SIMULATION_DEPTH */
#define DEF_MIN_SIMULATION_DEPTH 10000
/* This is the smallest possible simulation depth step */
#define MIN_SIMULATION_DEPTH_STEP 1
/* WARNING: Should be DEF_SIMULATION_DEPTH_STEP >= MIN_SIMULATION_DEPTH_STEP */
#define DEF_SIMULATION_DEPTH_STEP 1000

/**
* Sets the maximum simulation depth
* @param _max_sample_size the maximum simulation depth
*/
extern inline void setSimMaxSimulationDepth( int _max_sample_size );

/**
* Gets the maximum simulation depth
* @return the maximum simulation depth
*/
extern int getSimMaxSimulationDepth(void);

/**
* Sets the minimum simulation depth
* @param _min_sample_size the minimum simulation depth
*/
extern inline void setSimMinSimulationDepth( int _min_sample_size );

/**
* Gets the minimum simulation depth
* @return the minimum simulation depth
*/
extern int getSimMinSimulationDepth(void);

/**
* Sets the simulation depth step
* @param _sample_size_step the simulation depth step
*/
extern inline void setSimSimulationDepthStep( int _sample_size_step );

/**
* Gets the simulation depth step
* @return the simulation depth step
*/
extern int getSimSimulationDepthStep(void);

/****************************************************************************/
/*******************MANAGE THE BSCC MULTIPLIER PARAMETER*********************/
/*******************IS NEEDED FOR M.C. STEADY STATE**************************/
/*******************USING HYBRID SIMULATION APPROACH*************************/
/****************************************************************************/

/* This is the maximum BSCC dimension multiplier for choosing the sample-based regeneration point */
#define MAX_MULTIPLIER_SIZE 15

/* This is the minimum BSCC dimension multiplier for choosing the sample-based regeneration point */
#define MIN_MULTIPLIER_SIZE 2

/* This is the default BSCC dimension multiplier for choosing the sample-based regeneration point */
#define DEF_SAMPLE_BSCC_MULTIPLIER 3

/**
* Sets the BSCC dimension multiplier for choosing the sample-based regeneration point.
* @param _sampleBSCCDimensionMultiplier the BSCC dimension multiplier
*/
extern inline void setSampleBSCCDimensionMultiplier(int _sampleBSCCDimensionMultiplier);

/**
* Gets the BSCC dimension multiplier for choosing the sample-based regeneration point.
* @return the BSCC dimension multiplier
*/
extern int getSampleBSCCDimensionMultiplier(void);

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
extern inline void printRuntimeSimInfo(BOOL isExpNeeded);

#endif
