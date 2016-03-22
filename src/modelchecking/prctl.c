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
*	Authors: Maneesh Khattri
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
*	Source description: Solve formulas in PRCTL.
*/

#include "prctl.h"

#include "transient.h"

#include "runtime.h"

#include <string.h>


/*****************************************************************************
name		: ef_bounded: E[n}[subJ, supJ][PHI]
role		: solve N-R bounded E formula in DMRM
@param		: int steps: n.
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of N-R bounded E formula in DMRM.
remark		:
******************************************************************************/
static double * ef_bounded(int steps, const bitset * phi) {
	int i, j;
	const int size = get_state_space_size();
        double * result = (double *) calloc((size_t) size, sizeof(double)),
                        * temp_result;
        const double * rewards = getStateRewards();
	bitset *tt=get_new_bitset(size), *psi=get_new_bitset(size);
	fill_bitset_one(tt);

        /* get_idx_next_non_zero() is more efficient than checking every bit in
           phi individually. David N. Jansen. */
        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(phi, i)) != state_index_NONE ) {
			set_bit_val(psi, i, BIT_ON);
                        temp_result = until(TIME_INTERVAL_FORM, tt, psi, 0.0,
                                        (double) steps, FALSE);
                        if ( NULL == temp_result ) {
                                err_msg_3(err_CALLBY, "ef_bounded(%d,%p[%d])",
                                        steps, (const void *) phi,
                                        bitset_size(phi), NULL);
                        }
			for(j=0;j<size;j++) {
				result[j]+=temp_result[j]*rewards[i];
			}
			free(temp_result);
			set_bit_val(psi, i, BIT_OFF);
	}
	for(i=0;i<size;i++) {
		result[i]/=(steps+1);
	}
	free_bitset(psi);
	free_bitset(tt);
	return result;
}

/*****************************************************************************
name		: ef_unbounded: E[subJ, supJ][PHI]
role		: solve unbounded E formula in DMRM
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of N-R bounded E formula in DMRM.
remark		:
******************************************************************************/
static double * ef_unbounded(const bitset * phi) {
	int i, j;
	const int size = get_state_space_size();
        double * result = (double *) calloc((size_t) size, sizeof(double)),
                        * temp_result;
        const double * rewards = getStateRewards();
	bitset *tt=get_new_bitset(size), *psi=get_new_bitset(size);
	fill_bitset_one(tt);

        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(phi, i)) != state_index_NONE ) {
			set_bit_val(psi, i, BIT_ON);
			temp_result=until( TIME_UNBOUNDED_FORM, tt, psi, 0.0, 0.0, FALSE);
                        if ( NULL == temp_result ) {
                                err_msg_2(err_CALLBY, "ef_unbounded(%p[%d])",
                                        (const void *) phi, bitset_size(phi),
                                        NULL);
                        }
			for(j=0;j<size;j++){
				result[j]+=temp_result[j]*rewards[i];
			}
			free(temp_result);
			set_bit_val(psi, i, BIT_OFF);
	}
	free_bitset(psi);
	free_bitset(tt);
	return result;
}

/*****************************************************************************
name		: ef
role		: solve E formula in DMRM
@param		: int steps: n.
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of N-R bounded E formula in DMRM.
remark		:
******************************************************************************/
double * ef(int steps, const bitset * phi) {
	double *result;
	if( isRunMode(DMRM_MODE) )
	{
		if(steps==0)
                {
			result=ef_unbounded(phi);
                }
		else if(steps!=0)
                {
			result=ef_bounded(steps, phi);
                }
                if ( NULL == result ) {
                        err_msg_3(err_CALLBY, "ef(%d,%p[%d])", steps,
                                (const void *) phi, bitset_size(phi), NULL);
                }
	}
	else
	{
		printf("ERROR: Unsupported mode for E formula\n");
                exit(EXIT_FAILURE);
	}
	return result;
}

/*****************************************************************************
name		: cf
role		: solve C formula in DMRM
@param		: int steps: n.
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of N-R bounded C formula in DMRM.
remark		:
******************************************************************************/
double * cf(int steps, const bitset * phi) {
	int i, j;
        const
	sparse * state_space = get_state_space();
	const int size = get_state_space_size();
        double * result = (double *) calloc((size_t) size, sizeof(double));
        const double * rewards = getStateRewards();
        double * result_1 = (double *) calloc((size_t) size, sizeof(double));
        double * result_2 = (double *) calloc((size_t) size, sizeof(double)),
                        * pTmp = NULL;

	if( isRunMode(DMRM_MODE) ) {
                i = state_index_NONE;
                while ( (i = get_idx_next_non_zero(phi,i)) != state_index_NONE )
                {
                                memset(result_1, '\0', size * sizeof(double));
				result_1[i] = 1.0;
				/*Compute P^supi*i_phi*/
				for(j = 0; j < steps ; j++) {
                                        if ( err_state_iserror(multiply_mtx_MV(
                                                        state_space, result_1,
                                                        result_2)) )
                                        {
                                                err_msg_3(err_CALLBY, "cf(%d,"
                                                        "%p[%d])", steps,
                                                        (const void *) phi,
                                                        bitset_size(phi), NULL);
                                        }
					pTmp = result_1;
					result_1 = result_2;
					result_2 = pTmp;
				}
				for(j=0;j<size;j++) {
					result[j]+=result_1[j]*rewards[i];
				}
		}
	} else {
		printf("ERROR: Unsupported mode for C formula\n");
                exit(EXIT_FAILURE);
	}
	free(result_1);
	free(result_2);
	return result;
}

/*****************************************************************************
name		: yf
role		: solve Y formula in DMRM
@param		: int steps: n.
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of N-R bounded C formula in DMRM.
remark		:
******************************************************************************/
double * yf(int steps, const bitset * phi) {
	int i, j;
        const
	sparse * state_space = get_state_space();
	const int size = get_state_space_size();
        const
	double *rewards = getStateRewards();
        double * result1 = (double *) calloc((size_t) size, sizeof(double));
        double * result2 = (double *) calloc((size_t) size, sizeof(double));
        double * rewards2 = (double *) calloc((size_t) size, sizeof(double));

        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(phi, i)) != state_index_NONE ) {
			rewards2[i] = rewards[i];
	}
	if( isRunMode(DMRM_MODE) ) {
		for(i=0;i<steps;i++) {
                        if ( err_state_iserror(multiply_mtx_MV(state_space,
                                                        result1, result2)) )
                        {
                                err_msg_3(err_CALLBY, "yf(%d,%p[%d])", steps,
                                        (const void*)phi,bitset_size(phi),NULL);
                        }
			for(j=0;j<size;j++) {
				result1[j] = rewards2[j]+result2[j];
			}
		}
	} else {
		printf("ERROR: Unsupported mode for Y formula\n");
                exit(EXIT_FAILURE);
	}

	free(result2);
	free(rewards2);
	return result1;
}
