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
*	Source description: This is a source file for the sample-vector
*	structure and operations on it. It is used in the simulation engine.
*/


#ifndef SAMPLE_VEC_H
#define SAMPLE_VEC_H

#include "bitset.h"

	/****************************************************************************/
	/**************THE DEFINITION OF THE SAMPLE OF OBSERVATIONS TYPE*************/
	/****************************************************************************/

	/**
	* This is the common base of all sample vectors presented here
	* initial_state - the initial observation state index
	* curr_sample_size - the current sample size, i.e. the length of the pObservationsVec array
	* num_visited_states - the number of states visited during simulation of this sample
	*			we should also take into account the initial states (initial observations)
	*/
	typedef struct SSampleVec{
		int initial_state;
		int curr_sample_size;
		unsigned int num_visited_states;
	} TSampleVec;
	typedef TSampleVec * PTSampleVec;

	/**
	* This is a sample-vector structure, this one is used for the unbounded-until simulation.
	* sampleVecBase - the common part of all sample vectors
	* curr_simulation_depth - the current simulation depth needed or the unbounded until
	* sum_good - the current number of good absorbing states in the sample
	* sum_trans - the current number of transient states in the sample
	* pObservationsVec - the vector of observations, of size sampleVecBase.curr_sample_size
	* pTransStateInd - the bitset that contains the indexes of the transient states
	*			in the pObservationsVec vector.
	*/
	typedef struct SSampleVecUntil{
		TSampleVec sampleVecBase;	/* WARNING: Has to be the first field! */
		int curr_simulation_depth;
		int sum_good;
		int sum_trans;
		int * pObservationsVec;
		bitset * pTransStateInd;
	} TSampleVecUntil;
	typedef TSampleVecUntil * PTSampleVecUntil;

	/**
	* This is a sample-vector structure, this one is used for the interval-until simulation "Phi U[tl, tr] Psi"
	* sampleVecBase - the common part of all sample vectors
	* sum_good - the current number of observations that reached a Psi state within time interval [tl, tr]
	*/
	typedef struct SSampleVecIntUntil{
		TSampleVec sampleVecBase;	/* WARNING: Has to be the first field! */
		int sum_good;
	} TSampleVecIntUntil;
	typedef TSampleVecIntUntil * PTSampleVecIntUntil;

	/**
	* This is a sample-vector structure, this one is used for steady-state simulation.
	* NOTE: this structure is needed for every reachable BSCC!
	* sampleVecBase - the common part of all sample vectors
	*			Note that, sampleVecBase.curr_sample_size stores the number of regeneration cycles
	* pTimeObservationsVec - The vector of the expected time spent in the i'th regeneration cycle
	* pGoodObservationsVec - The vector of the accumulated number of "good" states divided by exit rates
	*			 in the i'th regeneration cycle.
	*/
	typedef struct SSampleVecSteady{
		TSampleVec sampleVecBase;	/* WARNING: Has to be the first field! */
		double * pTimeObservationsVec;
		double * pGoodObservationsVec;
	} TSampleVecSteady;
	typedef TSampleVecSteady * PTSampleVecSteady;

	/****************************************************************************/
	/************THE ALLOCATION AND DEALLOCATION OF THE SAMPLE VECTOR************/
	/****************************************************************************/

	/**
	* This method is used to allocate a new sample-vector structure for the
	* unbounded-until simulations CTMC/DTMC:
	* We set "curr_sample_size" to "sample_size".
	* We set "curr_simulation_depth" to zero.
	* We set "sum_good" to zero.
	* We set "sum_trans" to "sample_size".
	* We allocate "pTransStateInd" and "pObservationsVec" of size "sample_size".
	* We set "pTransStateInd" to all ones.
	* We fill "pObservationsVec" with "initial_state".
	* @param sample_size the initial sample size
	* @param initial_state the initial observation.
	* @return the newly-created sample vector.
	*/
	extern inline PTSampleVecUntil allocateSampleVectorUntil( const int sample_size, const int initial_state );

	/**
	* This method is used to allocate a new sample-vector structure for the
	* interval-until simulations CTMC/DTMC:
	* We set "curr_sample_size" to "sample_size".
	* We set "sum_good" to zero.
	* We allocate "pObsCurrStateVec" and "pObsExitTimeVec" of size "sample_size".
	* We fill "pObsCurrStateVec" with "initial_state".
	* We fill "pObsExitTimeVec" with zeroes.
	* @param sample_size the initial sample size
	* @param initial_state the initial observation.
	* @return the newly-created sample vector.
	*/
	extern inline PTSampleVecIntUntil allocateSampleVectorIntUntil( const int sample_size, const int initial_state );

	/**
	* This method is used to allocate a new sample-vector structure for the
	* steady-state simulations CTMC:
	* We set "initial_state" to "initial_state".
	* We set "curr_reg_cycles" to reg_cycles.
	* We set "pRegenerationPoints" to all zeros.
	* We set "pTimeObservationsVec" to all zeros.
	* We set "pGoodObservationsVec" to all zeros.
	* @param reg_cycles the initial number of regeneration cycles.
	* @param initial_state the regeneration point (the state index).
	* @return the newly-created sample vector.
	*/
	extern inline PTSampleVecSteady allocateSampleVectorSteady(const int reg_cycles, const int initial_state);

	/**
	* This method is used to extend a sample-vector structure for the
	* unbounded-until simulations CTMC/DTMC:
	* We can only increase the sample size, i.e. we expect
	*	"new_sample_size >= pSampleVecUntil->curr_sample_size"
	* @param pSampleVecUntil the sample-vector scheduled for reallocation.
	* @param sample_size the new sample size
	* @return the newly-created sample vector.
	*/
	extern inline void extendSampleVectorUntil( PTSampleVecUntil pSampleVecUntil, const int new_sample_size );

	/**
	* This method is used to extend a sample-vector structure for the
	* interval-until simulations CTMC/DTMC:
	* We can only increase the sample size, i.e. we expect
	*	"new_sample_size >= pSampleVecUntil->curr_sample_size"
	* @param pSampleVecIntUntil the sample-vector scheduled for reallocation.
	* @param sample_size the new sample size
	* @return the newly-created sample vector.
	*/
	extern inline void extendSampleVectorIntUntil( PTSampleVecIntUntil pSampleVecIntUntil, const int new_sample_size );

	/**
	* This method is used to extend a sample-vector structure for the
	* steady-state simulations CTMC:
	* We can only increase the regeneration cycles, i.e. we expect
	*	"new_reg_cycles >= pSampleVecSteady->curr_reg_cycles"
	* @param pSampleVecSteady the sample-vector scheduled for reallocation.
	* @param new_reg_cycles the new number of regeneration cycles.
	* @return the newly-created sample vector.
	*/
	extern inline void extendSampleVectorSteady( PTSampleVecSteady pSampleVecSteady, const int new_reg_cycles );

	/****************************************************************************/
	/***************FREE THE MEMORY ALLOCATED FOR THE SAMPLE VECTOR**************/
	/****************************************************************************/

	/**
	* This method is used to free the sample-vector structure for the
	* unbounded-until simulations CTMC/DTMC:
	* @param pSampleVecUntil the sample vector that has to be freed.
	*/
	extern inline void freeSampleVectorUntil( PTSampleVecUntil pSampleVecUntil );

	/**
	* This method is used to free the sample-vector structure for the
	* interval-until simulations CTMC/DTMC:
	* @param pSampleVecIntUntil the sample vector that has to be freed.
	*/
	extern inline void freeSampleVectorIntUntil( PTSampleVecIntUntil pSampleVecIntUntil );

	/**
	* This method is used to free the sample-vector structure for the steady-state
	* simulations CTMC:
	* @param pSampleVecSteady the sample vector that has to be freed.
	*/
	extern inline void freeSampleVectorSteady ( PTSampleVecSteady pSampleVecSteady );

	/****************************************************************************/
	/***************************PRINT THE SAMPLE VECTOR**************************/
	/****************************************************************************/

	/**
	* This method prints the sample-vector structure for unbounded-until simulation
	* @param pSampleVecUntil the sample-vector to be printed
	*/
	extern inline void printSampleVectorUntil( PTSampleVecUntil pSampleVecUntil );

	/**
	* This method prints the sample-vector structure for interval-until simulation CTMC/DTMC
	* @param pSampleVecIntUntil the sample-vector to be printed
	*/
	extern inline void printSampleVectorIntUntil( PTSampleVecIntUntil pSampleVecIntUntil );

	/**
	* This method prints the sample-vector structure for steady-state simulation CTMC
	* @param pSampleVecSteady the sample-vector to be printed
	*/
	extern inline void printSampleVectorSteady( PTSampleVecSteady pSampleVecSteady );

#endif
