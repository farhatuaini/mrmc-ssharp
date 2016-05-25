/**
*	WARNING: Do Not Remove This Section
*
*       $LastChangedRevision: 415 $
*       $LastChangedDate: 2016-03-23  $
*       $LastChangedBy: joleuger $
*
*	MRMC is a model checker for discrete-time and continuous-time Markov
*	reward models. It supports reward extensions of PCTL and CSL (PRCTL
*	and CSRL), and allows for the automated verification of properties
*	concerning long-run and instantaneous rewards as well as cumulative
*	rewards.
*
*	Copyright (C) The University of Twente, 2004-2008.
*	Copyright (C) RWTH Aachen, 2008-2009.
*	Copyright (C) University of Augsburg, 2016.
*	Authors: Johannes Leupolz
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
*	Source description: Writes a .res file containing the results of the last calculation.
*/

#include <stdlib.h>

#include "write_res_file.h"
#include "parser_to_tree.h"
#include "runtime.h"


#define DYNAMIC_ARRAY_INITIAL_CAPACITY 8
/* see https://en.wikipedia.org/wiki/Dynamic_array */
const char * res_file  = NULL;

struct write_res_file_list_of_states
{
	int *statesToWrite;
	size_t length;
	size_t capacity;
};

struct write_res_file_list_of_states write_res_file_statesToWrite;

void write_res_file_initialize(void)
{
	write_res_file_statesToWrite.statesToWrite = malloc(DYNAMIC_ARRAY_INITIAL_CAPACITY*sizeof(int));
	write_res_file_statesToWrite.length = 0;
	write_res_file_statesToWrite.capacity = DYNAMIC_ARRAY_INITIAL_CAPACITY;
}

void write_res_file_reset(void)
{
	write_res_file_statesToWrite.length = 0;
}

void write_res_file_add(int stateToWrite)
{
	int double_capacity;
	int nextIndex;
	
	nextIndex=write_res_file_statesToWrite.length;
	if(nextIndex == write_res_file_statesToWrite.capacity)
	{
		/* resize */
		double_capacity = write_res_file_statesToWrite.capacity * 2;
		write_res_file_statesToWrite.statesToWrite = (int*) realloc(write_res_file_statesToWrite.statesToWrite, double_capacity * sizeof(int));
	    write_res_file_statesToWrite.capacity = double_capacity;
	}
	write_res_file_statesToWrite.statesToWrite[nextIndex] = stateToWrite;
	write_res_file_statesToWrite.length++;
}

/**
* Derived from print_state_prob in runtime.c
*/
static void print_state_prob_to_file( FILE *p, const int index, const int size, const double * probs,
                                const double error_bound_local,
                                const double * pErrorBounds)
{
    const partition * P_local = getPartition();

	fprintf( p, "%d ", index + 1 );
	if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
			/* TODO
			print_state_prob_partition(P_local, index, probs,
							error_bound_local, pErrorBounds);
			*/
	}else{
		if( probs != NULL ){
			if( 0 <= index && index < size ){
				/* Calculate the right pattern according to the precision */
                                fprintf( p, "%1.*f\n", get_error_bound_precision(
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
}

/**
* Derived from print_state in runtime.c
*/
static void print_state_sat_to_file( FILE *p, const bitset * pBitset, const int index ){
	const partition * P_local = getPartition();

	fprintf(p, "%d ", index + 1);
       
	if ( isRunMode(F_IND_LUMP_MODE) && NULL != P_local ) {
		/*
		TODO
			print_original_state(P_local, pBitset, index);
		*/
	}else{
		if( pBitset != NULL ){
			if ( 0 <= index && index < bitset_size(pBitset) ) {
				if( get_bit_val( pBitset, index ) ){
					fprintf(p,"TRUE\n");
				}else{
					fprintf(p,"FALSE\n");
				}
			}else{
				fprintf(p,"??\n");
                                printf("WARNING: Invalid index %d, required to "
                                       "be in the [1, %d] interval.",
                                       index + 1, bitset_size(pBitset));
			}
		}else{
			fprintf(p,"??\n");
			printf("WARNING: Trying to print an element of a non-existing bitset.");
		}
	}
}

/**
* Derived from printResultingStateProbability in parser_to_core.c
*/
static void write_result_of_state_to_res_file(PTFTypeRes pFTypeRes, FILE *p, const int index){
	PTCompStateF pCompStateF; PTFTypeRes pFTypeResSubForm;
	
	/* Convert the user-level state index into the internal index */
	const int internal_state_index = index - 1;

	if( ( pFTypeRes != NULL ) && ( pFTypeRes->formula_type == COMPARATOR_SF ) ){
		pCompStateF = (PTCompStateF) pFTypeRes;
		pFTypeResSubForm = (PTFTypeRes) pCompStateF->unary_op.pSubForm;
		/* Check if the simulations are done here */
		if( pFTypeResSubForm->doSimHere ){
			/* WARNING: Here we assume that the conf. int. borders are computed exactly */
			/* TODO
			print_state_prob_to_file(p, internal_state_index, pFTypeResSubForm->prob_result_size,
						pFTypeResSubForm->pProbCILeftBorder, 0.0, NULL );
			print_state_prob_to_file(p, internal_state_index, pFTypeResSubForm->prob_result_size,
						pFTypeResSubForm->pProbCIRightBorder, 0.0, NULL );
			*/
		} else {
			print_state_prob_to_file(p, internal_state_index, pFTypeResSubForm->prob_result_size,
					pFTypeResSubForm->pProbRewardResult,
					pFTypeResSubForm->error_bound, pFTypeResSubForm->pErrorBound );
		}
	}else{
		printf("WARNING: There are NO results to print.\n");
	}
}

/**
* Derived from printResultingStateSatisfyability in parser_to_core.c
*/
static void write_satisfiability_of_state_to_res_file(PTFTypeRes pFTypeRes, FILE *p, const int index){
	PTCompStateF pCompStateF; PTFTypeRes pFTypeResSubForm;
	
	/* Convert the user-level state index into the internal index */
	const int internal_state_index = index - 1;

	if( pFTypeRes != NULL ){
		if( pFTypeRes->doSimHere || pFTypeRes->doSimBelow ){
			/* TODO:
			print_state_to_file(p, pFTypeRes->pYesBitsetResult, index - 1, YES_STATES_STR );
			print_state_to_file(p, pFTypeRes->pNoBitsetResult, index - 1, NO_STATES_STR );
			*/
		} else {
			print_state_sat_to_file(p, pFTypeRes->pYesBitsetResult, internal_state_index );
		}
	}else{
		printf("WARNING: There are NO results to print.\n");
	}
}


/*****************************************************************************
name		: write_res_file_state
role		: writes the res file with all requested states. Prints
              satisfiability.
@return		: void
remark		:
******************************************************************************/
void write_res_file_state() {	
	FILE *p;
	
	int i;
	int numberOfElements;
	
	numberOfElements = write_res_file_statesToWrite.length;
	
	PTFTypeRes pFTypeRes = (PTFTypeRes) get_formula_tree_result();
	
	printf("Writing results to file '%s'\n", res_file);
	
	p = fopen(res_file,"w+");
	
	for(i=0;i<numberOfElements;i++)
	{		
		int listElement = numberOfElements - i - 1; /* List is reverse */
		int state = write_res_file_statesToWrite.statesToWrite[listElement];
		/* printf("Write state %d\n",state); */
		write_satisfiability_of_state_to_res_file(pFTypeRes,p,state);
	}
	
	(void)fclose(p);
	write_res_file_reset();
}

/*****************************************************************************
name		: write_res_file_result
role		: writes the res file with all requested states. Prints calculated
              result (either probability or reward).
@return		: void
remark		:
******************************************************************************/
void write_res_file_result() {	
	FILE *p;
	
	int i;
	int numberOfElements;
	
	numberOfElements = write_res_file_statesToWrite.length;
	
	PTFTypeRes pFTypeRes = (PTFTypeRes) get_formula_tree_result();
	
	printf("Writing results to file '%s'\n", res_file);
	
	p = fopen(res_file,"w+");
	
	for(i=0;i<numberOfElements;i++)
	{		
		int listElement = numberOfElements - i - 1; /* List is reverse */
		int state = write_res_file_statesToWrite.statesToWrite[listElement];
		/* printf("Write state %d\n",state); */
		write_result_of_state_to_res_file(pFTypeRes,p,state);
	}
	
	(void)fclose(p);
	write_res_file_reset();
}
