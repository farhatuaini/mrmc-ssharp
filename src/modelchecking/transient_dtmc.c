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
*	Authors: Maneesh Khattri, Ivan Zapreev, Tim Kemna
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
*	Source description: Perform Transient Analysis for PCTL - X, U.
*/

#include "transient_dtmc.h"

#include "transient_common.h"
#include "lump.h"

#include "runtime.h"

/**
* Creates new R matrix and copies all the rows that belong to the pValidStates into it
* @param pStateSpace the initial matrix
* @param pP the new matrix to where to copy rows
*                     (the matrix is assumed to be empty)
* @param pValidStates the valid states i.e. ids of the rows to be copied
*                     this array contains the number of nodes as the first element
*/
static void makeAbsorbingDTMC(const sparse * pStateSpace, sparse * pP,
                const int * pValidStates)
{
	int i;
	int valid_state;

        /*Set diagonal values to 1 for all rows.*/
        const int size = mtx_rows(pStateSpace);
	for( i = 0; i < size  ; i++ )
	{
                if ( err_state_iserror(mtx_set_diag_val(pP, i, 1.0)) )
                        exit(err_macro_7(err_CALLBY,"makeAbsorbingDTMC(%p[%dx"
                                "%d],%p[%dx%d],%p)", (const void *) pStateSpace,
                                mtx_rows(pStateSpace), mtx_cols(pStateSpace),
                                (void *) pP, mtx_rows(pP), mtx_cols(pP),
                                (const void *) pValidStates, EXIT_FAILURE));
	}

	/*Form the pP matrix, here the 1 values of the pValidStates get correct values*/
	for( i = 1; i <= pValidStates[0] ; i++ )
	{
		valid_state = pValidStates[i];
                if ( err_state_iserror(mtx_copy_row(pP, valid_state,
                                                                pStateSpace)) )
                        exit(err_macro_7(err_CALLBY,"makeAbsorbingDTMC(%p[%dx"
                                "%d],%p[%dx%d],%p)", (const void *) pStateSpace,
                                mtx_rows(pStateSpace), mtx_cols(pStateSpace),
                                (void *) pP, mtx_rows(pP), mtx_cols(pP),
                                (const void *) pValidStates, EXIT_FAILURE));
	}
}

/**
* Solve the unbounded until operator for all states for DTMC.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @return	: double *: result of unbounded until operator for all states.
* NOTE: 1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
double * dtmc_unbounded_until(const bitset *phi, const bitset *psi)
{
	sparse *state_space = get_state_space();
	bitset * EU = get_exist_until(state_space, phi, psi);
        bitset * AU = NULL;
        int i, size = mtx_rows(state_space);
        double * initial = NULL, * rhs = NULL, * result = NULL;
        sparse * pP = NULL;
	bitset *dummy=NULL;
	int *pValidStates=NULL;

        const char * error_str = err_CALLBY;
        if ( NULL == EU
                        || (AU = get_always_until(state_space, phi, psi, EU))
                                        == NULL
                        || (error_str = err_MEMORY,
                            initial = (double *) calloc((size_t) size,
                                        sizeof(double))) == NULL
                        || (rhs = (double *) calloc((size_t) size,
                                        sizeof(double))) == NULL
                        || (error_str = err_CALLBY,
                            pP = allocate_sparse_matrix(size, size)) == NULL
        /* not(not(EU) or AU) = EU and not(AU). This is more efficient than the
           previous implementation. David N. Jansen */
                        || (dummy = not(AU)) == NULL
                        || err_state_iserror(and_result(EU, dummy))
                        || (pValidStates = count_set(dummy)) == NULL )
        {
                err_msg_4(error_str, "dtmc_unbounded_until(%p[%d],%p[%d])",
                        (const void *) phi, bitset_size(phi), (const void*) psi,
                        bitset_size(psi),
                        ((void) (NULL == EU || (free_bitset(EU),
                            NULL == AU || (free_bitset(AU),
                                NULL == initial || (free(initial),
                                    NULL == rhs || (free(rhs),
                                        NULL == pP || (free_mtx_sparse(pP),
                                            NULL==dummy || (free_bitset(dummy),
                                                free(pValidStates),FALSE))))))),
                        NULL));
        }
	/* make ~(Phi /\ ~Psi) states absorbing */
	makeAbsorbingDTMC(state_space, pP, pValidStates);

        i = state_index_NONE;
        /* get_idx_next_non_zero() is more efficient than checking every bit of
           AU individually. David N. Jansen */
        while ( (i = get_idx_next_non_zero(AU, i)) != state_index_NONE ) {
			rhs[i] = 1.0;
			initial[i] = 1.0;
	}
	/* solve (I-P)x = i_Psi */
	result = unbounded_until_universal(pP, pValidStates, initial, rhs, TRUE);

	free(rhs);
	free_bitset(EU);
	free_bitset(AU);
	free_bitset(dummy);
	free(pValidStates);
        if ( err_state_iserror(free_mtx_wrap(pP)) ) {
                err_msg_4(err_CALLBY, "dtmc_unbounded_until(%p[%d],%p[%d])",
                        (const void *) phi, bitset_size(phi),
                        (const void *) psi, bitset_size(psi), NULL);
        }

	return result;
}

/**
* Solve the unbounded until operator for all states for DTMC with lumping.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @return		: double *: result of unbounded until operator for all states.
*/
double * dtmc_unbounded_until_lumping(const bitset *phi, const bitset *psi)
{
        const
	sparse *state_space = get_state_space();
        int size = mtx_rows(state_space);
	bitset * EU = get_exist_until(state_space, phi, psi);
        bitset * AU = NULL;
        double * initial = NULL, * rhs = NULL, * result, * lumped_result = NULL;
        bitset * dummy = NULL, * lumped_dummy = NULL;
	int *pValidStates=NULL;
        partition * P = NULL;
        sparse * pP = NULL;

        const char * error_str = err_CALLBY;
        if ( NULL == EU
                        || (AU = get_always_until(state_space, phi, psi, EU))
                                        == NULL
	/* create initial partition for unbounded until formula and lump */
                        || (P = init_partition_formula(EU, AU, FALSE)) == NULL
                        || (pP = lump(P, state_space)) == NULL

        /* not(not(EU) or AU) = EU and not(AU) */
                        || (dummy = not(AU)) == NULL
                        || err_state_iserror(and_result(EU, dummy))
                        || (lumped_dummy = lump_bitset(P, dummy)) == NULL
                        || (pValidStates = count_set(lumped_dummy)) == NULL

	/* initialize for lumped DTMC: initial and rhs vectors */
                        || (error_str = err_MEMORY,
                            initial = (double *) calloc((size_t) mtx_rows(pP),
                                        sizeof(double))) == NULL
                        || (rhs = (double*)calloc((size_t) mtx_rows(pP),
                                        sizeof(double))) == NULL )
        {
                err_msg_4(error_str, "dtmc_unbounded_until_lumping(%p[%d],%p["
                        "%d])", (const void *) phi, bitset_size(phi),
                        (const void *) psi, bitset_size(psi),
                        ((void) (NULL == EU || (free_bitset(EU),
                           NULL == AU || (free_bitset(AU),
                              NULL == P || (free_partition(P),
                                 NULL == pP || (free_sparse_ncolse(pP),
                                    NULL == dummy || (free_bitset(dummy),
                                       NULL == lumped_dummy || (free_bitset(
                                                                lumped_dummy),
                                          NULL == pValidStates || (free(
                                                                pValidStates),
                                             NULL == initial || (free(initial),
                                                free(rhs), FALSE))))))))),
                        NULL));
        }
        part_walk_blocks(P, B) {
                state_index lumped_id = state_index_NONE;
                state_index orig_id = part_unlump_state_block(P, B);
                if ( state_index_ERROR == orig_id
                                || (get_bit_val(AU, orig_id)
                                    && (lumped_id = part_lump_state_block(P, B))
                                                == state_index_ERROR) )
                {
                        exit(err_macro_4(err_CALLBY,
                                "dtmc_unbounded_until_lumping(%p[%d],%p[%d])",
                                (const void *) phi, bitset_size(phi),
                                (const void *) psi, bitset_size(psi),
                                EXIT_FAILURE));
                }
                if ( state_index_NONE != lumped_id ) {
                        rhs[lumped_id] = 1.0;
                        initial[lumped_id] = 1.0;
		}
	}
        end_part_walk_blocks;

	/* solve (I-P)x = i_Psi */
	lumped_result = unbounded_until_universal(pP, pValidStates, initial, rhs, FALSE);

	/* transform result vector */
	result = unlump_vector(P, size, lumped_result);

	free(rhs);
	free(lumped_result);
	free_bitset(EU);
	free_bitset(AU);
	free_bitset(dummy);
	free_bitset(lumped_dummy);
	free_partition(P);
	free(pValidStates);
        if ( err_state_iserror(free_sparse_ncolse(pP)) ) {
                err_msg_4(err_CALLBY, "dtmc_unbounded_until_lumping(%p[%d],%p"
                        "[%d])", (const void *) phi, bitset_size(phi),
                        (const void *) psi, bitset_size(psi),
                        (free(result), NULL));
        }

	return result;
}

/**
* The common part of the dtmc_bounded_until method for
* the case with and without lumping.
* @param pP the transition matrix
* @param ppInOutData the pointer to the pointer of the array of values.
*			  We have to do it this way since then the *ppInOutData
*			  array is used later and we can free it as result_2 here.
*			  So the remaining result_1 array is assigned to *ppInOutData.
* @param supi the time bound - number of iterations.
*/
static void dtmc_bounded_until_universal(sparse *pP, double **ppInOutData, double supi)
{
	int i;
	double *result_2, *pTmp;
	double * result_1 = *ppInOutData; /* Access the data-array pointer value */

	/*Compute Q^supi*i_psi*/
        result_2 = (double *) calloc((size_t) mtx_rows(pP), sizeof(double));
	for(i = 1; i <= supi ; i++) {
                if ( err_state_iserror(multiply_mtx_MV(pP, result_1,result_2)) )
                {
                        exit(err_macro_5(err_CALLBY,
                                "dtmc_bounded_until_universal(%p[%dx%d],%p,%g)",
                                (void *) pP, mtx_rows(pP), mtx_cols(pP),
                                (void *) ppInOutData, supi, EXIT_FAILURE));
                }
		pTmp = result_1;
		result_1 = result_2;
		result_2 = pTmp;
	}

	/* Free allocated memory */
	free(result_2);
        /* IMPORTANT! WE NEED THE RIGHT POINTER BACK! */
	*ppInOutData = result_1; /* Store the remaining pointer */
}

/**
* Solve the bounded until operator for DTMC.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @param		: double supi: sup I should contain the Natural number
* @return	: double *: result of the unbounded until operator for all states.
* NOTE: 1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
double * dtmc_bounded_until(const bitset *phi, const bitset *psi, double supi)
{
        const
	sparse * state_space = get_state_space();
	sparse *pP;
	int i;
        const int size = mtx_rows(state_space);
        double * result_1 = (double *) calloc((size_t) mtx_rows(state_space),
                        sizeof (double));
	int * pValidStates;

	/* Do optimization, we exclude states from Phi which are Psi states and */
        /* also states from which you always go to bad i.e. not Phi ^ not Psi
           states */
	/* This is only possible if until has lower time bound 'subi' */
	/* (and reward bound 'subj' if any) equal to 0. */
	bitset *good_phi_states = get_good_phi_states( phi, psi, state_space);

	pValidStates = count_set(good_phi_states);

	pP = allocate_sparse_matrix(size,size);
        if ( NULL == pP ) {
                err_msg_5(err_CALLBY, "dtmc_bounded_until(%p[%d],%p[%d],%g)",
                        (const void *) phi, bitset_size(phi), (const void*) psi,
                        bitset_size(psi), supi, (free_bitset(good_phi_states),
                        free(pValidStates), free(result_1), NULL));
        }
	makeAbsorbingDTMC(state_space, pP, pValidStates);

	/*init result_1 with i_psi*/
        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(psi, i)) != state_index_NONE ) {
			result_1[i] = 1.0;
	}

	/*Compute Q^supi*i_psi*/
	dtmc_bounded_until_universal(pP, &result_1, supi);

	/* Free allocated memory */
	free_bitset(good_phi_states);
	free(pValidStates);
        if ( err_state_iserror(free_mtx_wrap(pP)) ) {
                err_msg_5(err_CALLBY, "dtmc_bounded_until(%p[%d],%p[%d],%g)",
                        (const void *) phi, bitset_size(phi),
                        (const void *) psi, bitset_size(psi), supi,
                        (free(result_1), NULL));
        }

	return result_1;
}

/**
* Solve the bounded until operator for DTMC with lumping.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @param		: double supi: sup I should contain the Natural number
* @return		: double *: result of the unbounded until operator for all states.
*/
double * dtmc_bounded_until_lumping(const bitset *phi, const bitset *psi, double supi)
{
        const
	sparse * state_space = get_state_space();

	/* Do optimization, we exclude states from Phi which are Psi states and */
        /* also states from which you always go to bad i.e. not Phi ^ not Psi
           states */
	/* This is only possible if until has lower time bound 'subi' */
	/* (and reward bound 'subj' if any) equal to 0. */
	bitset *good_phi_states = get_good_phi_states( phi, psi, state_space);

	/* create initial partition for bounded until formula and lump */
	partition *P = init_partition_formula(good_phi_states, psi, FALSE);
	sparse *Q = lump(P, state_space);

	/*init lumped_result with i_psi*/
        double * lumped_result = (double *) calloc((size_t) mtx_rows(Q),
                        sizeof(double));
	double *result;

        part_walk_blocks(P, B) {
                state_index lumped_id = state_index_NONE;
                state_index orig_id = part_unlump_state_block(P, B);
                if ( state_index_ERROR == orig_id
                                || (get_bit_val(psi, orig_id)
                                    && (lumped_id = part_lump_state_block(P, B))
                                                == state_index_ERROR) )
                {
                        err_msg_5(err_CALLBY,
                                "dtmc_bounded_until_lumping(%p[%d],%p[%d],%g)",
                                (const void *) phi, bitset_size(phi),
                                (const void *) psi, bitset_size(psi), supi,
                                (free_bitset(good_phi_states),free_partition(P),
                                free_sparse_ncolse(Q), free(lumped_result),
                                NULL));
                }
                if ( state_index_NONE != lumped_id ) {
                        lumped_result[lumped_id] = 1.0;
                                                /* if B is labeled with psi */
		}
	}
        end_part_walk_blocks;

	/*Compute Q^supi*i_psi*/
	dtmc_bounded_until_universal(Q, &lumped_result, supi);

	/* transform result vector */
        result = unlump_vector(P, mtx_rows(state_space), lumped_result);

	/* Free allocated memory */
	free(lumped_result);
	free_partition(P);
	free_bitset(good_phi_states);
        if ( err_state_iserror(free_sparse_ncolse(Q)) ) {
                err_msg_5(err_CALLBY, "dtmc_bounded_until_lumping(%p[%d],%p"
                        "[%d],%g)", (const void *) phi, bitset_size(phi),
                        (const void *) psi, bitset_size(psi), supi,
                        (free(result), NULL));
        }

	return result;
}

/**
* This method modelchecks the Xphi formula for the DTMC
* @param the phi states
* @return the probabilities
*/
double * dtmc_unbounded_next(const bitset * phi)
{
        const
	sparse *state_space = get_state_space();
        const int size = mtx_rows(state_space);
	int i;
        double * result = (double *) calloc((size_t) size, sizeof(double));
        double * i_phi = (double *) calloc((size_t) size, sizeof(double));

        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(phi, i)) != state_index_NONE ) {
			i_phi[i] = 1.0;
	}

        if ( err_state_iserror(multiply_mtx_MV(state_space, i_phi, result)) ) {
                err_msg_2(err_CALLBY, "dtmc_unbounded_next(%p[%d])",
                                (const void *) phi, bitset_size(phi),
                                (free(i_phi), free(result), NULL));
        }

	free(i_phi);

	return result;
}
