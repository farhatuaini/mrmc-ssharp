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
*	Authors: Maneesh Khattri, Ivan Zapreev
*
*       Copyright (C) (this file) Saarland University, 2007-2009
*       Author: Ernst Moritz Hahn
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
*       Our contacts (Saarland University).
*               Ernst Moritz Hahn, Dependable Systems and Software
*               Campus Saarbruecken
*               66123 Saarbruecken
*               Phone: +49 681 302 5607
*               Fax: +49 681 302 5636
*               Location:   Bldg. E1 3, Room 536
*               E-Mail: emh@cs.uni-sb.de
*
*       Source description: Contains the algorithm for bounded time
*       reachability for  CDMDPs
*/

#include "transient_ctmdpi_hd_uni.h"

#include "foxglynn.h"

#include "runtime.h"

#include <math.h>

/**
* Checks for uniformity, sets common rate and divides by it.
*
* @param sparse calculate E from here
* @return TRUE iff uniform
*/
static BOOL uniformize_ctmdpi(NDSparseMatrix * sparse_local)
{
	double E = 0.0;
	double E_primed;
	double delta;
	BOOL result = TRUE;
	unsigned i;
	int state_nr;
	unsigned choice_nr;
        unsigned * row_starts = (unsigned *) sparse_local->row_counts;
        unsigned * choice_starts = (unsigned *) sparse_local->choice_counts;
        double * non_zeros = sparse_local->non_zeros;
	const double epsilon = get_error_bound();

	/* calculate reference rate */
	/* Take the first non-zero rate we find.
	 * As soon as we have a situation where this is not the case, we
	 * leave the loops.
	 */
        for ( state_nr = 0 ; state_nr < sparse_local->n ; state_nr++ ) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			/* Add up all outgoing rates of the distribution */
			unsigned i_start = choice_starts[choice_nr];
			unsigned i_end = choice_starts[choice_nr + 1];
			E = 0.0;
			for (i = i_start; i < i_end; i++) {
				E = E + non_zeros[i];
			}
			if (E > 0.0) {
				break;
			}
		}
		if (E > 0.0) {
			break;
		}
	}

        for ( state_nr = 0 ; state_nr < sparse_local->n ; state_nr++ ) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			/* Add up all outgoing rates of the distribution */
			unsigned i_start = choice_starts[choice_nr];
			unsigned i_end = choice_starts[choice_nr + 1];
			E_primed = 0.0;
			for (i = i_start; i < i_end; i++) {
				E_primed = E_primed + non_zeros[i];
			}
			/* if non-zero rate differs more than epsilon
			 * from reference rate this is incorrect */
			delta = fabs(E - E_primed);
			if ((E_primed > 0) && (delta > epsilon)) {
				result = FALSE;
				break;
			}
		}
		if (FALSE == result) {
			break;
		}
	}

	if (FALSE == result) {
		printf("MDP is not uniform, cannot model check with "
		       "current method.\n");
                return FALSE;
	}

	/* set rate and divide by it */
        sparse_local->rate = E;
        for ( state_nr = 0 ; state_nr < sparse_local->n ; state_nr++ ) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			unsigned i_start = choice_starts[choice_nr];
			unsigned i_end = choice_starts[choice_nr + 1];
			for (i = i_start; i < i_end; i++) {
				non_zeros[i] /= E;
			}
		}
	}

	return TRUE;
}

/**
* Transforms uniformized DTMDP to CTMDPI again.
*
* @param sparse calculate E from here
* @return TRUE iff uniform
*/
static void deuniformize_ctmdpi(NDSparseMatrix * sparse_local)
{
	unsigned i;
	int state_nr;
	unsigned choice_nr;
        unsigned * row_starts = (unsigned *) sparse_local->row_counts;
        unsigned * choice_starts = (unsigned *) sparse_local->choice_counts;
        double * non_zeros = sparse_local->non_zeros;

	/* multiply by rate */
        for ( state_nr = 0 ; state_nr < sparse_local->n ; state_nr++ ) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			unsigned i_start = choice_starts[choice_nr];
			unsigned i_end = choice_starts[choice_nr + 1];
			for (i = i_start; i < i_end; i++) {
                                non_zeros[i] *= sparse_local->rate;
			}
		}
	}
}

/**
 * Calculates matrix P(s,alpha,B).
 *
 * The result will be a one-dimensional array. Entries will be in the
 * following order: The first entry is for the first state, first
 * non-deterministic selection, then the one for the second
 * non-deterministic decision follows and so forth. Then the ones for
 * the second state follow, etc.
 *
 * @param sparse sparse matrix to calculate for
 * @param B set of states in B
 * @return P(s,alpha,B)
 */
static double * calculate_P_s_alpha_B(const NDSparseMatrix * sparse_local,
                const bitset * B)
{
        double *P_s_alpha_B;
	unsigned ps_index;

        unsigned * row_starts = (unsigned *) sparse_local->row_counts;
        unsigned * choice_starts = (unsigned *) sparse_local->choice_counts;
        double * non_zeros = sparse_local->non_zeros;
        unsigned * cols = sparse_local->cols;

	/* calculate final size of P(s,alpha,B) */
	unsigned ps_size = 0;
        state_index state_nr;
        bitset * not_B = not(B);

        /* get_idx_next_non_zero() is more efficient than testing every bit in
           B individually. David N. Jansen. */
        state_nr = state_index_NONE;
        while ( (state_nr = get_idx_next_non_zero(not_B, state_nr))
                                != state_index_NONE )
        {
			ps_size += row_starts[state_nr + 1] - row_starts[state_nr];
	}

	/* allocate and fill P(s,alpha,B). We precompute it for later
	 * usage in main part of the algorithm. */
        P_s_alpha_B = (double *) malloc(ps_size * sizeof(double));
	ps_index = 0;

        state_nr = state_index_NONE;
        while ( (state_nr = get_idx_next_non_zero(not_B, state_nr))
                                != state_index_NONE )
        {
		BITSET_BLOCK_TYPE s_in_B;
			unsigned state_start = row_starts[state_nr];
			unsigned state_end = row_starts[state_nr + 1];
			unsigned choice_nr;
			for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
				unsigned i;
				unsigned i_start = choice_starts[choice_nr];
				unsigned i_end = choice_starts[choice_nr + 1];
				P_s_alpha_B[ps_index] = 0.0;
				for (i = i_start; i < i_end; i++) {
					s_in_B = get_bit_val(B, cols[i]);
					if ( !(s_in_B == BIT_OFF) ) {
						P_s_alpha_B[ps_index] += non_zeros[i];
					}
				}
				ps_index++;
			}
	}

        free_bitset(not_B);

	return P_s_alpha_B;
}

/**
* Returns new optimal value for m.
* If @a min is true, it will return the smaller number of @a m and @a
* m_primed, if false, the larger one is returned.
*
* @param min true for minimum, false for maximum
* @param m first number
* @param b second number
* @return optimum of numbers
*/
static double optimum(BOOL min, double a, double b) {
	double result = a;

	if (min) {
		if (b < a) result = b;
	} else {
		if (b > a) result = b;
	}

	return result;
}

/**
* Get probability to reach target states with given choice.
* (for current number of jumps k considered in outer loop,
* P_lambda(X=k))
*
* @param j choice number
* @param choice_starts see CTMDPI documentation
* @param non_zeros see CTMDPI documentation
* @param cols see CTMDPI documentation
* @param q q (see paper)
* @param P_s_alpha_B P(s,alpha,B) (see paper)
* @param ps_index index to P_s_alpha_B
* @param psi Poisson probability for given nr. step nr
*/
inline static double get_choice_probability_psi( const unsigned j, const int *choice_starts, const double *non_zeros,
						const unsigned *cols, const double *q, const double *P_s_alpha_B,
						unsigned ps_index, const double psi) {
	double result = 0.0;
	unsigned l2 = choice_starts[j];
	unsigned h2 = choice_starts[j+1];
	unsigned k;
	for (k = l2; k < h2; k++) {
		result += non_zeros[k] * q[cols[k]];
	}
	result += psi * P_s_alpha_B[ps_index];

	return result;
}

/**
* Get probability to reach target states with given choice.
* (for current number of jumps k considered in outer loop,
* P_lambda(X=k))
*
* @param j choice number
* @param choice_starts see CTMDPI documentation
* @param non_zeros see CTMDPI documentation
* @param cols see CTMDPI documentation
* @param q q (see paper)
*/
inline static double get_choice_probability( const unsigned j, const int *choice_starts, const double *non_zeros,
						const unsigned *cols, const double *q) {
	double result = 0.0;
	unsigned l2 = choice_starts[j];
	unsigned h2 = choice_starts[j+1];
	unsigned k;
	for (k = l2; k < h2; k++) {
		result += non_zeros[k] * q[cols[k]];
	}

	return result;
}

/**
* Swaps two double vectors.
*
* @param a first double vector
* @param b second double vector
*/
inline static void swap(double **a, double **b) {
	double *swap_var = *a;
	*a = *b;
	*b = swap_var;
}

/**
* Actual CTMDPI bounded reachability algorithm.
*
* @param unsigned num_states number of states of CTMDPI
* @param fg Poisson values by Fox-Glynn algorithm
* @param left_end use Poisson probabilities up to left limit
* @param min true for minimal probability, false for maximal one
* @param row_starts see CTMDPI documentation
* @param choice_starts see CTMDPI documentation
* @param non_zeros see CTMDPI documentation
* @param q current probability vector
* @param q_primed next probabilities vector
* @param cols see CTMDPI documentation
* @param P_s_alpha_B P(s,alpha,B) (see paper)
* @param B target states
*/
static double *ctmdpi_iter( const unsigned num_states, const FoxGlynn *fg, const unsigned left_end,
                                const BOOL min, const int * row_starts,
                                const int * choice_starts,
				const double *non_zeros, double *q, double *q_primed, const unsigned *cols,
				const double *P_s_alpha_B, const bitset *B) {
       /* in the algorithm the main iteration is in lines 3-12. For each
	* jump probability k starting from the right bound in the Fox-Glynn
	* algorithm down to 1 the loop body is to be executed once. In the
	* implementation here we split this loop into two parts: One for
	* the non-negligible Poisson probabilities between left and right
	* bound of the Fox-Glynn algorithm and one for the jump
	* probabilities below left bound, if any.
	*/
	unsigned i;

        state_index row;
        bitset * not_B;

	/* first part - use poisson probabilities */

	/* line 3 in paper */
	for (i = fg->right; i >= left_end; i--) {
		unsigned ps_index = 0;
		double psi = fg->weights[i-fg->left] / fg->total_weight;
		/* lines 4-12 in paper */
                for ( row = 0 ; (unsigned) row < num_states ; row++ ) {
			BITSET_BLOCK_TYPE s_in_B;
			s_in_B = get_bit_val(B, row);
			if ( s_in_B == BIT_OFF ) {
				/* get maximizing/minimizing decision */
				/* lines 4-10 in paper (non-target state) */
                                double m = min ? 2.0 : -1.0;
				unsigned l1 = row_starts[row];
				unsigned h1 = row_starts[row+1];
				unsigned j;
				for (j = l1; j < h1; j++) {
					double m_primed = get_choice_probability_psi
					(j, choice_starts, non_zeros, cols,q, P_s_alpha_B, ps_index, psi);
					m = optimum(min, m, m_primed);
					ps_index++;
				}
				/* in case some non-target state did
				 * not have any possible decisions we
				 * set m to zero afterwards */
				m = ((-1.0 == m) || (2.0 == m)) ? 0.0 : m;
				/* line 9 in paper*/
				q_primed[row] = m;
			} else {
				/* line 11 in paper (target state) */
				q_primed[row] = psi + q[row];
			}
		}
		/* swapping is done instead using a new q vector for
		 * each k */
		swap(&q, &q_primed);
	}

	/* we have to set the q_primed values for target states once
	 * more now; they won't change afterwards, as psi =
	 * 0.0. We do here what would be done in the first iteration
	 * of the next loop. As the values for q_primed do not
	 * influence the first iteration of the next loop, we can do
	 * this beforehand to avoid having to do this in each loop. */
        row = state_index_NONE;
        while ( (row = get_idx_next_non_zero(B, row)) != state_index_NONE )
        { /* target states only */
	    q_primed[row] = q[row];
	}

        not_B = not(B);
	/* now do part where we are below left bound of poisson probabilities */
	for (; i > 0; i--) {
		unsigned ps_index = 0;
		/* lines 4-12 in paper */
                row = state_index_NONE;
                while ( (row = get_idx_next_non_zero(not_B, row))
                                        != state_index_NONE )
                { /* non-target state? */
				/* get maximizing/minimizing decision */
				/* lines 4-10 in paper (non-target state) */
				double m = min ? 2 : -1.0;
				unsigned l1 = row_starts[row];
				unsigned h1 = row_starts[row+1];
				unsigned j;
				for (j = l1; j < h1; j++) {
					double m_primed = get_choice_probability
					(j, choice_starts, non_zeros, cols, q);
					m = optimum(min, m, m_primed);
					ps_index++;
				}
				/* in case some non-target state did
				 * not have any possible decisions we
				 * set m to zero afterwards */
				m = ((-1.0 == m) || (2.0 == m)) ? 0.0 : m;
				/* line 9 in paper*/
				q_primed[row] = m;
		}
		/* swapping is done instead using a new q vector for
		 * each k */
		swap(&q, &q_primed);
	}

        free_bitset(not_B);

	return q;
}

/**
 * Uniform history dependent CTMDPI time bounded reachability.
 * Algorithm taken from
 *
 * "Efficient computation of time-bounded reachability probabilities in uniform
 * continuous-time Markov decision processes"
 *
 * by C. Baier, B. Haverkort, H. Hermanns, J.-P. Katoen.
 *
 * @author Ernst Moritz Hahn (emh@cs.uni-sb.de)
 * @param t time bound
 * @param B reachability set
 * @return the calculated probability vector
 */
double *transient_ctmdpi_hd_uni(const double t, const bitset *B) {
        state_index row;
	FoxGlynn *fg;
	unsigned left_end;
	double *result;
	unsigned state;
	double E;
	const double u = get_underflow();
	const double o = get_overflow();
	const double epsilon = get_error_bound();
	double *P_s_alpha_B;
	/* if user typed a "<" or "<=" as comparator in the fomula, we
	 * compute the maximal probability, for ">" or ">=" the
	 * minimal one. These are the worst-case probabilities that
	 * is, if it is possible for some state not to fulfill the
	 * formula, it won't fulfill it with these probabilities. */
	int comparator = get_comparator();
        BOOL min = (C_GREATER == comparator) || (C_GREATER_EQUAL == comparator);
        const
	double *non_zeros;
        const
	int *row_starts;
        const
	int *choice_starts;
        const
	unsigned int *cols;
	double *q;
	double *q_primed;

	/* some shortcuts for several values (precision, # of states,
	 * etc.) */
	NDSparseMatrix *ctmdp = get_mdpi_state_space();
	unsigned num_states = ctmdp->n;
	if (!uniformize_ctmdpi(ctmdp)) {
	  /* in case of failure, return dummy result */
	  result = (double *)calloc(num_states, sizeof(double));
	  return result;
	}
	E = ctmdp->rate;

	P_s_alpha_B = calculate_P_s_alpha_B(ctmdp, B);

	/* store local pointers of stuff */
	non_zeros = ctmdp->non_zeros;
	row_starts = (int *) ctmdp->row_counts;
	choice_starts = (int *) ctmdp->choice_counts;
	cols = ctmdp->cols;

	/* In case Fox-Glynn algorithm terminated successfully, we proceed */
	if ( fox_glynn(E * t, u, o, epsilon, &fg) ) {
		printf("Fox-Glynn: ltp = %d, rtp = %d, w = %1.15e\n",
		 fg->left, fg->right, fg->total_weight);

		/* we have to fill the entries for all states with zero, even
		 * for target states.
		 * Corresponds to line 2 or algorithm in paper. */
                q = (double *) calloc((size_t) num_states, sizeof(double));
                q_primed = (double*)calloc((size_t) num_states, sizeof(double));

		/* we have to do this, because for very small rates the FG-algo also
		 * gives the possibility to not doing a jump at all. */
		left_end = (0 == fg->left) ? 1 : fg->left;

		/* start main iteration of algorithm */
		result = ctmdpi_iter( num_states, fg, left_end, min, row_starts, choice_starts,
					non_zeros, q, q_primed, cols, P_s_alpha_B, B );

		/* line 15 in paper */
                row = state_index_NONE;
                while ( (row=get_idx_next_non_zero(B,row)) != state_index_NONE )
                {
				result[row] = 1.0;
		}
                /* cleanup stuff */
                if (result != q){
                        free(q); q = NULL;
                }
                if (result != q_primed){
                        free(q_primed); q_primed = NULL;
                }
	} else {
		/* Otherwise return a dummy result */
                result = (double *) calloc((size_t) num_states, sizeof(double));
		for (state = 0; state < (unsigned) ctmdp->n; state++) {
		  result[state] = get_bit_val(B, state);
		}
	}

	/* cleanup stuff */
	free(P_s_alpha_B);
	freeFG(fg);

	deuniformize_ctmdpi(ctmdp);

	/* line 16 in paper */
	return result;
}
