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
*	Source description: Contains some common methods for transient
*				 analysis of PCTL/CSL/PRCTL/CSRL - X, U.
*/

#include "transient_common.h"

#include "iterative_solvers.h"

#include "runtime.h"

/**
* Solve E(phi U psi) until formula.
* @param: sparse *state_space: the state space
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: bitset *psi: satisfaction relation for psi formula.
* @return: bitset *: result of E(SAT(phi) U SAT(psi)) for all states.
* NOTE: adapted for t-bounded-until from
*         1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
bitset * get_exist_until(const sparse *state_space, const bitset *phi, const bitset *psi)
{
        int i, size, reach;
	int * states = NULL;
        bitset * EU;

        if ( NULL == state_space || NULL == phi || NULL == psi
                        || mtx_rows(state_space) != bitset_size(psi)
                        || bitset_size(phi) != bitset_size(psi) )
        {
                err_msg_7(err_PARAM, "get_exist_until(%p[%dx%d],%p[%d],%p[%d])",
                        (const void *) state_space,
                        NULL != state_space ? mtx_rows(state_space) : 0,
                        NULL != state_space ? mtx_cols(state_space) : 0,
                        (const void *) phi, NULL != phi ? bitset_size(phi) : 0,
                        (const void *) psi, NULL != psi ? bitset_size(psi) : 0,
                        NULL);
        }

        size = bitset_size(phi);
        /* The following two statements are a more efficient way to copy psi to
           EU than the previously programmed detour via a third bitset. David N.
           Jansen */
        EU = get_new_bitset(size);
        copy_bitset(psi, EU);

	/* Find out and store all psi states first in the states structure */
        /* Calling count_set to get a list of states is more efficient than
           testing every bit of psi individually. David N. Jansen. */
        states = count_set(psi);
	/* The reach variable will contain the number of psi states */
        reach = states[0];
        /* Note that we will not update states[0], but only reach to reflect
           changes in size of states[]. */

	/*  */
        for ( i = 1 ; i <= reach ; i++ )
	{
                mtx_walk_column_nodiag_noval(state_space, back_set_j,
                                                        (const int) states[i])
                {
			/* If we can reach psi state from the phi state */
			/* and it's not yet been found then */
                        if ( get_bit_val(phi, back_set_j)
                                        && ! get_bit_val(EU, back_set_j) )
                        {
				/* Add this phi state to the set of states from which psi is reachable */
                                states = (int *) realloc(states,
                                                (reach + 2) * sizeof(int));
                                states[++reach] = back_set_j;
				/* Mark this state as the one from which psi is reachable. */
                                set_bit_val(EU, back_set_j, BIT_ON);
			}
		}
                end_mtx_walk_column_nodiag_noval;
	}
	free(states); states=NULL;
	return EU;
}

/**
* Solve A(phi U psi) until formula.
* @param: sparse *state_space: the state space
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: bitset *psi: satisfaction relation for psi formula.
* @param: bitset *e_phi_psi: The indicator set for the formula E(Phi U Psi)
* @return: bitset *: result of A(SAT(phi) U SAT(psi)) for all states.
* NOTE: adapted for t-bounded-until from
*         1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
bitset * get_always_until(const sparse *state_space, const bitset *phi, const bitset *psi, const bitset *e_phi_psi)
{
        int i, k, size = bitset_size(phi), * toremove = NULL, pres_size;
        BOOL flag;
	/* Copy the E(Phi U Psi) into the initial A(Phi U Psi) bitset */
        bitset * AU = get_new_bitset(size);
        copy_bitset(e_phi_psi, AU);

	/* find a state in E(Phi U Psi) that does not satisfy A(Phi U Psi) */
        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(AU, i)) != state_index_NONE ) {
                flag = FALSE;
		/* If it is a pure phi state */
                if ( ! get_bit_val(psi, i) ) {
			/* Check if we can go from the i'th state to some non psi state */
			/* If so then the i'th state and all of its predecessors will */
			/* have to be excluded from AU */
                        mtx_walk_row_nodiag(state_space, (const int) i, col_j,
                                        UNUSED(val))
                        {
                                if ( ! get_bit_val(AU, col_j) ) {
                                        flag = TRUE;
					break;
				}
			}
                        end_mtx_walk_row_nodiag;

			if( flag ) {
				/* remove i'th state and all of its predecessors */
				pres_size = 1;
                                toremove = (int *) calloc((size_t) pres_size,
                                                sizeof(int));
				toremove[0] = i;
				set_bit_val(AU, toremove[pres_size-1], BIT_OFF);
				for(k=0;k<pres_size;k++) {
					/* schedule predecessors for removal */
                                        mtx_walk_column_nodiag_noval(
                                                        state_space, back_set_l,
                                                        (const int) toremove[k])
                                        {
                                                if ( get_bit_val(AU, back_set_l)
                                                        && !get_bit_val(psi,
                                                                back_set_l) )
                                                {
							toremove=(int *)realloc(toremove, (pres_size+1)*sizeof(int));
                                                        toremove[pres_size++] =
                                                                back_set_l;
							set_bit_val(AU, toremove[pres_size-1], BIT_OFF);
						}
					}
                                        end_mtx_walk_column_nodiag_noval;
				}
				free(toremove);
				toremove = NULL;
			}
		}
	}
	return AU;
}

/**
* Universal part of PCTL and CSL unbounded until
* Solves the system of linear equations Ax=b
* @param: pM: the A matrix
* @param: pValidStates: this array contains the number of nodes as the
*	    first element all the other elements are the node ids, if it is
*	    NULL then all the nodes from the A matrix are valid
* @param: pX: the initial x vector.
*	NOTE: Is possibly freed inside by solveGaussJacobi method
* @param: pB: the b vector
* @param: BOOL revert_matrix: revert the matrix afterwards, if true
* @return: the solution of the system
*/
double * unbounded_until_universal(sparse * pM, int * pValidStates, double * pX,
                const double * pB, BOOL revert_matrix)
{
	double err = get_error_bound();
	int max_iterations = get_max_iterations();
	int method;
	double *pResult;

	/* get (I-P) */
        if ( err_state_iserror(mult_mtx_cer_const(pM, -1.0, pValidStates))
                        || err_state_iserror(add_cer_cons_diagonal_mtx(pM,
                                        pValidStates, 1)) )
        {
                err_msg_7(err_CALLBY, "unbounded_until_universal(%p[%dx%d],%p,"
                        "%p,%p,%d)", (void *) pM, mtx_rows(pM), mtx_cols(pM),
                        (void *) pValidStates, (void *) pX, (const void *) pB,
                        revert_matrix, NULL);
        }

	/* solve (I-P)x = i_Psi */
	if( get_method_path() == GJ) /* Retrieve the l.e. solution method from the runtime.c */
		method = GAUSS_JACOBI;
	else /* This should be "GS" otherwise */
		method = GAUSS_SEIDEL;

	pResult = solve_initial(method, pM, pX, pB, err, max_iterations, pValidStates);

	if(revert_matrix){
                if ( err_state_iserror(mult_mtx_cer_const(pM, -1.0,
                                                        pValidStates)) )
                {
                        err_msg_7(err_CALLBY, "unbounded_until_universal(%p[%d"
                                "x%d],%p,%p,%p,%d)", (void *) pM, mtx_rows(pM),
                                mtx_cols(pM), (void *) pValidStates, (void*) pX,
                                (const void *) pB, revert_matrix,
                                (free(pResult), NULL));
                }
	}

	return pResult;
}

/**
* Do optimization, we exclude states from Phi which are Psi states and
* also states from which you always go to bad i.e. not Phi ^ not Psi states
* NOTE: This optimization is applicable for until if only it has lower time bound 'subi'
* (and reward bound 'subj' if any) equal to 0. The latter is because Psi states should
* not be made absorbing before the lower time bound 'subi' is reached!
* So it should be used with time (and reward) bounded until, NOT interval until,
* though there it may be applied in a smart way as well.
* @param phi the Phi states
* @param psi the Psi states
* @param state_space this is the sparse matrix
* @return phi_and_not_psi without phi states from which you never reach psi states
*/
bitset * get_good_phi_states(const bitset * phi, const bitset * psi,
                const sparse * state_space)
{
        /* It is not necessary to allocate another bitset by calling and() and
           freeing the temporary bitset not(psi) afterwards. David N. Jansen */
        bitset * eu_reach_psi_through_phi = get_exist_until(state_space, phi,
                        psi);
        bitset * good_phi_states = NULL;

        if ( NULL == eu_reach_psi_through_phi
                        || (good_phi_states = not(psi)) == NULL
                        ||err_state_iserror(and_result(eu_reach_psi_through_phi,
                                        good_phi_states))
	/* Free memory */
                        || (err_state_iserror(free_bitset(
                                        eu_reach_psi_through_phi))
                                && (eu_reach_psi_through_phi = NULL, TRUE)) )
        {
                err_msg_7(err_CALLBY, "get_good_phi_states(%p[%d],%p[%d],%p[%dx"
                        "%d])", (const void *) phi, bitset_size(phi),
                        (const void *) psi, bitset_size(psi),
                        (const void *) state_space, mtx_rows(state_space),
                        mtx_cols(state_space),
                        (NULL == eu_reach_psi_through_phi
                                    || (free_bitset(eu_reach_psi_through_phi),
                                        FALSE),
                        NULL == good_phi_states
                                    || (free_bitset(good_phi_states), FALSE),
                        NULL));
        }

	return good_phi_states;
}

/**
* Make certain states (not in n_absorbing) absorbing.
* Creates a new empty sparse matrix, then assighs non absorbing rows via pointers
* and computes uniformization rates for them.
* @param: sparse *state_space: the state space
* @param: bitset *n_absorbing: not absorbing states.
* @param: double *pLambda(ref.): return value of uniformizing rate(also q).
* @param: double *abse: return row_sums
* @param: int    *pNonAbsorbing(ref.): return value for the number of
*                       non-absorbing states in the model
* @return: sparse *: new sparse matrix with some states made absorbing.
* NOTE: makes states absorbing (in ab_state_space) by assigning pointers
*	  to rows (in state_space) only for those rows which are not to be
*	  made absorbing.
* NOTE: initially *lambda should be equal to 0!!!!
******************************************************************************/
sparse * ab_state_space(const sparse * state_space, const bitset *n_absorbing, double *pLambda, double *abse, int *pNonAbsorbing)
{
	int i, n, non_absorbing = *pNonAbsorbing;
	const double *e;
    double lambda = *pLambda;
        sparse * ab_state_space_result;

	if ( lambda != 0.0 ) {
		printf("ERROR: lambda parameter in ab_state_space is not 0!");
                exit(EXIT_FAILURE);
	}
        n = bitset_size(n_absorbing);
	/*Get row sums of the initial matrix*/
        e = get_row_sums();     /* Note: DO NOT FREE the array you have
                                   acquired! */
						/* It is globally stored in runtime.c */
        ab_state_space_result = allocate_sparse_matrix(mtx_rows(state_space),
                        mtx_cols(state_space)); /*Init a new matrix structure*/
        if ( NULL == ab_state_space_result ) {
                err_msg_8(err_CALLBY, "ab_state_space(%p[%dx%d],%p[%d],%p,%p,"
                        "%p)", (const void*) state_space, mtx_rows(state_space),
                        mtx_cols(state_space), (const void *) n_absorbing,
                        bitset_size(n_absorbing), (void*) pLambda, (void*) abse,
                        (void *) pNonAbsorbing, NULL);
        }
        i = state_index_NONE;
        while( (i = get_idx_next_non_zero(n_absorbing, i)) != state_index_NONE )
        { /*For all non-absorbing states*/
			++non_absorbing;
			/*printf("lambda = %e, e[%d] = %e \n", lambda, i, e[i] );*/
			if( lambda < e[i] ){
				lambda = e[i]; /*Get max among available row sums*/
			}
			abse[i] = e[i];
                        /*Assign rows via pointers*/
                        if ( err_state_iserror(mtx_copy_row(
                                        ab_state_space_result, i,state_space)) )
                        {
                                err_msg_8(err_CALLBY,"ab_state_space(%p[%dx%d],"
                                        "%p[%d],%p,%p,%p)",
                                        (const void *) state_space,
                                        mtx_rows(state_space),
                                        mtx_cols(state_space),
                                        (const void *) n_absorbing,
                                        bitset_size(n_absorbing),
                                        (void *) pLambda, (void *) abse,
                                        (void *) pNonAbsorbing,
                                        (free_mtx_wrap(ab_state_space_result),
                                         NULL));
                        }
	}
        /* NOTE: I think we do not need aperiodicity here, since there are just
           absorbing and transient states! */
	/* It should not affect the steady state detection method since there the steady state is detected by the */
	/* amount of the probability mass left in the transient states. The steady state might be reached faster though, */
	/* because the number of self loops in the transient part of MRMC becomes smaller with this lambda. */
	/*  */
	/* if( lambda > 1.0 ) ++(lambda); //Make uniformization rate strictly > than max row sum */

    *pLambda = lambda;
    *pNonAbsorbing = non_absorbing;
        return ab_state_space_result;
}
