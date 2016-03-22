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
*       Copyright (C) (this file) Saarland University, 2009
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
*       reachability for non-uniform CDMDPs
*/

#include "transient_ctmdpi_hd_non_uni.h"

#include "foxglynn.h"
#include "runtime.h"

#include <string.h>
#include <math.h>

/*
 * TODO
 * - efficient computation of hypoexpontial distributions for smaller k
 *   (space / time tradeoff)
 * - write test cases?
 * - precision guarantees?
 */

inline static void swap_double(double **, double **);
inline static double compute_max_rate(const double *, const unsigned);

/*****************Tools for hypoexponential distribution *****************/

/**
 * Data structure for computation of hypoexponential distributions.
 */
typedef struct hypoexp
{
	const double *rates;
	unsigned num_rates;
	FoxGlynn *fg;
	double max_rate;
	double *selfloop_probs;
	double *forward_probs;
	double *reach_probs;
	double *weighted_probs;
} hypoexp_t;

/**
 * prepare chain-structured markov chain to compute values
 *
 * @param rates_mult rates vector (corresponds to i)
 * @param num_rates size of rates vector
 * @param hypoexp hypoexp data structure
 */
static void build_hypoexp_chain
(unsigned *rates_mult, unsigned num_stages, hypoexp_t *hypoexp)
{
	unsigned multiplicity;
	unsigned stage;
	unsigned rate_nr;
	double selfloop_prob;
	double forward_prob;
	unsigned mult_nr;
	double rate;
	const double max_rate = hypoexp->max_rate;

	memset(hypoexp->reach_probs, 0, (num_stages+1) * sizeof(double));
	memset(hypoexp->weighted_probs, 0, (num_stages+1) * sizeof(double));

	hypoexp->reach_probs[0] = 1.0;
	stage = 0;
	for (rate_nr = 0; rate_nr < hypoexp->num_rates; rate_nr++) {
		multiplicity = rates_mult[rate_nr];
		if (multiplicity > 0) {
			rate = hypoexp->rates[rate_nr];
			selfloop_prob = (max_rate - rate) / max_rate;
			forward_prob = rate / max_rate;
			for (mult_nr = 0; mult_nr < multiplicity; mult_nr++) {
				hypoexp->selfloop_probs[stage] = selfloop_prob;
				hypoexp->forward_probs[stage] = forward_prob;
				stage++;
			}
		}
	}
	hypoexp->selfloop_probs[num_stages] = 1.0;
}

/**
 * Computes value of cumulative distribution function from given vector.
 * Corresponds to the l_i(t) from the paper.
 *
 * @param rates_mult rates vector (corresponds to i)
 * @param num_rates size of rates vector
 * @param hypoexp hypoexp data structure
 */
static double hypoexp_cumulate
(unsigned *rates_mult, unsigned num_stages, hypoexp_t *hypoexp)
{
	unsigned iter_nr;
	double weight;
	double result;
	unsigned stage;
	unsigned left_start;
	const FoxGlynn *fg = hypoexp->fg;

	if (0 == num_stages) {
		return 1.0;
	}

	build_hypoexp_chain(rates_mult, num_stages, hypoexp);

	for (iter_nr = 1; iter_nr < (unsigned) fg->left; iter_nr++) {
		for (stage = num_stages; stage > 0; stage--) {
			hypoexp->reach_probs[stage] =
				hypoexp->reach_probs[stage]
				* hypoexp->selfloop_probs[stage]
				+ hypoexp->reach_probs[stage - 1]
				* hypoexp->forward_probs[stage - 1];
		}
		hypoexp->reach_probs[0] =
			hypoexp->reach_probs[0] * hypoexp->selfloop_probs[0];
	}

	left_start = (0 == fg->left ? 1 : fg->left);
	if (0 == fg->left) {
		hypoexp->weighted_probs[0] += fg->weights[0];
	}
	for (iter_nr = left_start; iter_nr <= (unsigned) fg->right; iter_nr++) {
		weight = fg->weights[iter_nr-fg->left];
		for (stage = num_stages; stage > 0; stage--) {
			hypoexp->reach_probs[stage] =
				hypoexp->reach_probs[stage]
				* hypoexp->selfloop_probs[stage]
				+ hypoexp->reach_probs[stage - 1]
				* hypoexp->forward_probs[stage - 1];
			hypoexp->weighted_probs[stage] +=
				hypoexp->reach_probs[stage] * weight;
		}
		hypoexp->reach_probs[0] =
			hypoexp->reach_probs[0] * hypoexp->selfloop_probs[0];
		hypoexp->weighted_probs[0] +=
			hypoexp->reach_probs[0] * weight;
	}

	result = hypoexp->weighted_probs[num_stages] / fg->total_weight;

	return result;
}

/**
 * Poisson distribution from max of @a rates times @a time_bound.
 *
 * @param rates vector of rates
 * @param num_rates size of @a rates
 * @param time_bound time bound for reachability
 * @param fg is set to Poisson values
 */
static void compute_fg_probs
(const double *rates, unsigned num_rates, double time_bound, FoxGlynn **fg)
{
	const double max_rate = compute_max_rate(rates, num_rates);
	const double u = get_underflow();
	const double o = get_overflow();
	const double precision = get_error_bound();

	if (!fox_glynn(max_rate * time_bound, u, o, precision, fg)) {
		printf("Fox-Glynn algorithm failed\n");
		exit(1);
	}
}

/**
 * Initialize data structure for computation of hypoexp distributions.
 * Notice that @a rates is not "consumed" by this function. That is,
 * we only set a pointer to the data structure, but freeing it is left
 * to the creator of @a rates. Also, please don't free @a rates when still
 * using the result of this function.
 *
 * @param max_stages maximal number of stages ever used
 * @param rates vector of different rates
 * @param num_rates size of @a rates
 * @param time_bound time bound for following analyses
 */
static hypoexp_t *init_hypoexp
(unsigned max_stages, const double *rates, unsigned num_rates, double time_bound)
{
	hypoexp_t *hypoexp = malloc(sizeof(hypoexp_t));
	FoxGlynn *fg;

	hypoexp->rates = rates;
	hypoexp->num_rates = num_rates;
	compute_fg_probs(rates, num_rates, time_bound, &fg);
	hypoexp->fg = fg;
	hypoexp->max_rate = compute_max_rate(rates, num_rates);
	hypoexp->selfloop_probs = (double*)malloc((max_stages + 1) * sizeof(double));
	hypoexp->forward_probs = (double*)malloc(max_stages * sizeof(double));
	hypoexp->reach_probs = (double *)calloc((max_stages + 1), sizeof(double));
	hypoexp->weighted_probs = (double *)calloc(max_stages + 1, sizeof(double));

	return hypoexp;
}

/**
 * Free hypoexp data structure.
 * Notice that you have to free @a hypoexp->rates yourself.
 *
 * @param hypoexp data structure to be freed
 */
static void free_hypoexp(hypoexp_t *hypoexp)
{
	freeFG(hypoexp->fg);
	free(hypoexp->selfloop_probs);
	free(hypoexp->forward_probs);
	free(hypoexp->reach_probs);
	free(hypoexp->weighted_probs);
	free(hypoexp);
}

/***************** main part *****************/

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
static double optimum(BOOL min, double a, double b)
{
	double result = a;

	if (min) {
		if (b < a) result = b;
	} else {
		if (b > a) result = b;
	}

	return result;
}

/**
 * Swaps two double vectors.
 *
 * @param a first double vector
 * @param b second double vector
 */
inline static void swap_double(double **a, double **b)
{
	double *swap_var = *a;
	*a = *b;
	*b = swap_var;
}

/**
 * Swaps two unsigned integer vectors.
 *
 * @param a first unsigned integer vector
 * @param b second unsigned integer vector
 */
inline static void swap_unsigned(unsigned **a, unsigned **b)
{
	unsigned *swap_var = *a;
	*a = *b;
	*b = swap_var;
}

/**
 * Computes maximum rate in given rate vector
 *
 * @param rates vector of rates
 * @param num_rates size of rates vector @a rates
 * @return maximum of rates in given vector
 */
static inline double compute_max_rate
(const double *rates, const unsigned num_rates)
{
	double max_rate;
	unsigned rate_nr;

	max_rate = 0;
	for (rate_nr = 0; rate_nr < num_rates; rate_nr++) {
		max_rate = rates[rate_nr] > max_rate ? rates[rate_nr] : max_rate;
	}

	return max_rate;
}

/**
 * Computes the k-bound needed to guarantee prespecified precision.
 *
 * @param rates vector of rates occuring
 * @param num_rates size of rates vector @a rates
 * @param time_bound time bound to consider
 * @return k-bound needed for given precision
 */
static unsigned compute_k_bound
(double *rates, unsigned num_rates, double time_bound)
{
	const double precision = get_error_bound();
	const double max_rate = compute_max_rate(rates, num_rates);

	return (unsigned) ceil(max_rate*time_bound*exp(2.0) - log(precision/2.0));
}

/**
 * Computes numbers of choices in CTMDP.
 *
 * @param ctmdp CTMDP to compute number of choices of
 * @return number of choices
 */
static unsigned compute_num_choices(const NDSparseMatrix *ctmdp)
{
	unsigned state_nr;
	const unsigned *row_starts = (unsigned *) ctmdp->row_counts;

	unsigned num_choices = 0;
	for (state_nr = 0; state_nr < (unsigned) ctmdp->n; state_nr++) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		num_choices += state_end - state_start;
	}

	return num_choices;
}

/**
 * Prepare CTMDP for algorithm execution.
 * The steps done by this function are
 * <ul>
 *   <li>construct vector of all occuring different rates</li>
 *   <li>assign choices in CTMDP to corresponding vector entries</li>
 *   <li>embed choices, so we have probabilities instead of rates</li>
 * </ul>
 *
 * @param ctmdp CTMDP to prepare for transient algorithm
 * @param rates rates vector to be constructed
 * @param num_rates resulting size of rates vector
 * @param choices_to_rates assigns choices to rate vector entries
 */
static void prepare_ctmdpi
(NDSparseMatrix *ctmdp, double **rates, unsigned *num_rates,
 unsigned **choices_to_rates)
{
	unsigned choice_nr;
	unsigned rate_nr;
	unsigned state_nr;
	unsigned i;
	const unsigned num_choices = compute_num_choices(ctmdp);
	const double precision = get_error_bound();
	const unsigned *row_starts = (unsigned *) ctmdp->row_counts;
	const unsigned *choice_starts = (unsigned *) ctmdp->choice_counts;
	double *non_zeros = ctmdp->non_zeros;

	*rates = (double*) malloc(num_choices * sizeof(double));
	*choices_to_rates = (unsigned*) malloc(num_choices * sizeof(double));
	*num_rates = 0;

	for (state_nr = 0; state_nr < (unsigned) ctmdp->n; state_nr++) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			/* Add up all rates of a choice */
			unsigned choice_start = choice_starts[choice_nr];
			unsigned choice_end = choice_starts[choice_nr + 1];
			double choice_rate = 0.0;
			for (i = choice_start; i < choice_end; i++) {
				choice_rate += non_zeros[i];
			}
			for (rate_nr = 0; rate_nr < *num_rates; rate_nr++) {
				double delta = fabs(choice_rate - (*rates)[rate_nr]);
				if (delta < precision) {
					break;
				}
			}
			if (rate_nr == *num_rates) {
				/* rate not in rate vector, so add it */
				(*num_rates)++;
				(*rates)[rate_nr] = choice_rate;
			}
			(*choices_to_rates)[choice_nr] = rate_nr;
			for (i = choice_start; i < choice_end; i++) {
				non_zeros[i] /= choice_rate;
			}
		}
	}
}

/**
 * Undoes preparations for CTMDP for algorithm execution.
 * The steps done by this function are
 * <ul>
 *   <li>un-embed choices, so we have rates again</li>
 *   <li>free @a rates and @a choices_to_rates</li>
 * </ul>
 *
 * @param ctmdp CTMDP to prepare for transient algorithm
 * @param rates rates vector
 * @param choices_to_rates assigns choices to rate vector entries
 */
static void unprepare_ctmdpi
(NDSparseMatrix *ctmdp, double *rates, unsigned *choices_to_rates)
{
	double choice_rate;
	unsigned choice_nr;
	unsigned state_nr;
	unsigned i;
	const unsigned *row_starts = (unsigned *) ctmdp->row_counts;
	const unsigned *choice_starts = (unsigned *) ctmdp->choice_counts;
	double *non_zeros = ctmdp->non_zeros;

	choice_nr = 0;
	for (state_nr = 0; state_nr < (unsigned) ctmdp->n; state_nr++) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			unsigned i_start = choice_starts[choice_nr];
			unsigned i_end = choice_starts[choice_nr + 1];
			choice_rate = rates[choices_to_rates[choice_nr]];
			for (i = i_start; i < i_end; i++) {
				non_zeros[i] *= choice_rate;
			}
		}
	}

	free(rates);
	free(choices_to_rates);
}

/**
 * Computes binomial coefficient @a n over @a k.
 * The implementation takes care that numbers don't become
 * unnecessarily large during the computation to avoid
 * overflow.
 *
 * @param n "over" parameter
 * @param k "under" parameter
 * @return @a n over @a k
 */
static unsigned binom(unsigned n, unsigned k)
{
	unsigned result;
	unsigned i;

	if (0 == k) {
		return 1;
	}

	if (2*k > n) {
		result = binom(n, n - k);
	} else {
		result = n;
		for (i = 2; i <= k; i++) {
			result *= (n + 1 - i);
			result /= i;
		}
	}

	return result;
}

/**
 * Return multiset coefficient.
 *
 * @param n "over" parameter
 * @param k "under" parameter
 * @return @a n over @a k
 */
static unsigned multiset_coeff(unsigned n, unsigned k)
{
	return binom(n + k - 1, k);
}

/**
 * Push back vector to given vector field.
 * Places @a vector at given position @a next_index of @a vectors and
 * increases @a next_index. This way, all vectors occuring can be
 * memorized by called the function one time for each of them.
 *
 * @param vectors memory field of vectors
 * @param vector vector to memorize
 * @param next_index where to place vector, is incremented
 */
static void push_back_vector
(unsigned *vectors, unsigned *vector,
 unsigned vec_size, unsigned *next_index)
{
	memcpy(vectors + (*next_index) * vec_size,
	       vector, vec_size * sizeof(unsigned));
	(*next_index)++;
}

/**
 * Helper function for enumerate_k_vectors.
 *
 * @param new_entry vector using during recursion
 * @param level at which column of given vector are we
 * @param rest rest of cardinality to be placed in vector
 * @param vec_size size of vector
 * @param next_vector_index where to place next computed vector
 * @param vectors memory where to place vectors
 */
static void enumerate_k_vectors_rec
(unsigned *new_entry, unsigned level, unsigned rest, unsigned vec_size,
 unsigned *next_index, unsigned *vectors)
{
	unsigned i;

	if ((level+1 < vec_size) && (rest > 0)) {
		for (i = 0; i <= rest; i++) {
			new_entry[level] = i;
			enumerate_k_vectors_rec(new_entry, level+1, rest-i, vec_size,
						next_index, vectors);
			new_entry[level] = 0;
		}
	} else {
		new_entry[level] = rest;
		push_back_vector(vectors, new_entry, vec_size, next_index);
		new_entry[level] = 0;
	}
}

/**
 * Enumerate nonnegative vectors of given length and cardinality.
 * The result will be placed in @a vectors, which must have been
 * previously allocated. For this, a memory amount of "num_rates
 * multichoose cardinality" * sizeof(unsigned) will be needed. The
 * vectors will be placed into @a vectors one after another. Vectors
 * will be placed into @a vectors in lexicographical order.
 *
 * @param num_rates number of rates
 * @param cardinality cardinality of vectors computed
 * @param vectors here vectors will be placed.
 * @param num_vectors number of vectors computed
 */
static void enumerate_k_vectors
(const unsigned num_rates, const unsigned cardinality, unsigned *vectors,
 unsigned *num_vectors)
{
	unsigned next_vector_index;
	unsigned *new_entry;

	next_vector_index = 0;
	new_entry = calloc(num_rates, sizeof(unsigned));
	enumerate_k_vectors_rec(new_entry, 0, cardinality, num_rates,
				&next_vector_index, vectors);
	*num_vectors = next_vector_index;
	free(new_entry);
}

/**
 * Compares two vectors.
 * Returns -1 if @a vec1 < @a vec2, 1 if  @a vec1 > @a vec2 and
 * 0 if  @a vec1 = @a vec2. For comparing vectors, the lexicographic
 * is used.
 *
 * @param vec1 first vector
 * @param vec2 second vector
 * @param size sizes of vectors
 * @return -1 if @a vec1 < @a vec2, 1 if  @a vec1 > @a vec2, 0 else
 */
static int compare_vectors(unsigned *vec1, unsigned *vec2, unsigned size)
{
	unsigned i;

	for (i = 0; i < size; i++) {
		if (vec1[i] < vec2[i]) {
			return -1;
		} else if (vec1[i] > vec2[i]) {
			return 1;
		}
	}

	return 0;
}

/**
 * Computes number of given @a vector by finding it in @a vectors.
 * Searches for @a vector in @a vectors. @a vectors should consist of
 * a bunch of vectors one after another in memory, in lexicographic
 * order. Returned is the entry number to which @a vector corresponds.
 *
 * @param vector vector to compute number of
 * @param vectors vectors to search @a vector in
 * @param vec_length length of vectors
 * @return number of vector
 */
static unsigned vector_to_number
(unsigned *vector, unsigned *vectors, unsigned vec_length,
 unsigned left, unsigned right)
{
	unsigned mid = left;
	int comp_res;
	unsigned *comp_vec;
        BOOL found;

        found = FALSE;
	while (!found && (left <= right)) {
		mid = left + ((right - left) / 2);
		comp_vec = vectors + mid * vec_length;
		comp_res = compare_vectors(vector, comp_vec, vec_length);
		if (0 == comp_res) {
                        found = TRUE;
		} else if (-1 == comp_res) {
			right = mid - 1;
		} else {
			left = mid + 1;
		}
	}

	return mid;
}

/**
 * Get vector of given number from vector array.
 *
 * @param vectors vector array
 * @param number number of vector to return
 * @param size size of vectors
 * @return vector with given @a number
 */
static inline unsigned *get_ith_vec
(unsigned *vectors, const unsigned number, const unsigned size) {
	return vectors + (number * size);
}

/**
 * Compute state probabilities for largest k to consider.
 *
 * @param all_k_rates_quant all vectors of k bound cardinality
 * @param num_k_vecs number of vectors in @a all_k_rates_quant
 * @param num_rates length of rate vectors
 * @param hypoexp hypoexp data structure
 * @param k_bound largest k to consider
 * @param num_states number of states of CTMDP
 * @param B set of target states
 * @param probs probabilities to be computed
 */
static void compute_k_bound_probs
(unsigned *all_k_rates_quant, const unsigned num_k_vecs,
 const unsigned num_rates, hypoexp_t *hypoexp, unsigned k_bound,
 unsigned num_states, const bitset *B, double *probs)
{
	unsigned prob_nr;
	unsigned vec_nr;
	unsigned *rates_quant;
	double li;
	unsigned state_nr;

	prob_nr = 0;
	fflush(stdout);
	for (vec_nr = 0; vec_nr < num_k_vecs; vec_nr++) {
		rates_quant = get_ith_vec(all_k_rates_quant, vec_nr, num_rates);
		li = hypoexp_cumulate(rates_quant, k_bound, hypoexp);

		for (state_nr = 0; state_nr < num_states; state_nr++) {
			if (get_bit_val(B, state_nr)) {
				probs[prob_nr] = li;
			} else {
				probs[prob_nr] = 0.0;
			}
			prob_nr++;
		}
	}
}

/**
 * Compute transient reachability probabilities in nonuniform CTMDPs.
 *
 * @param ctmdp CTMDP to compute probabilities for
 * @param rates vector of different rates occuring in CTMDP
 * @param num_rates size of rates vector
 * @param choices_to_rates maps choice numbers to corresponding rate index
 * @param k_bound maximal cardinality to consider
 * @
 */
static double *ctmdpi_iterate
(const NDSparseMatrix *ctmdp, const double *rates, const unsigned num_rates,
 const unsigned *choices_to_rates, const unsigned k_bound, const bitset *B,
 const double time_bound, const BOOL min)
{
	double *probs_k;
	double *probs_km1;
	unsigned state_nr;
	double li;
	unsigned prob_nr;
	unsigned vec_nr;
	unsigned num_k_vecs;
	unsigned num_km1_vecs;
	unsigned *rates_quant;
	unsigned *all_k_rates_quant;
	unsigned *all_km1_rates_quant;
	unsigned k;
	const unsigned num_states = (unsigned) ctmdp->n;
	const unsigned max_entries = multiset_coeff(num_rates, k_bound);
	unsigned succ_nr;
	unsigned choice_nr;
	double choice_prob;
	double succ_prob;
	unsigned succ_vec_nr;
	const unsigned *row_starts = (unsigned *) ctmdp->row_counts;
	const unsigned *choice_starts = (unsigned *) ctmdp->choice_counts;
	unsigned rate_nr;
	double *non_zeros = ctmdp->non_zeros;
	hypoexp_t *hypoexp;
	unsigned *rate_to_succ_vec_nr;

	probs_k = malloc(max_entries * num_states * sizeof(double));
	probs_km1 = malloc(max_entries * num_states * sizeof(double));
	all_k_rates_quant = malloc(max_entries * num_rates * sizeof(unsigned));
	all_km1_rates_quant = malloc(max_entries * num_rates * sizeof(unsigned));
	enumerate_k_vectors(num_rates, k_bound, all_k_rates_quant, &num_k_vecs);
	hypoexp = init_hypoexp(k_bound, rates, num_rates, time_bound);
	rate_to_succ_vec_nr = malloc(num_rates * sizeof(unsigned));

	compute_k_bound_probs(all_k_rates_quant, num_k_vecs, num_rates,
			      hypoexp, k_bound, num_states, B, probs_k);

	for (k = k_bound; k > 0; k--) {
		enumerate_k_vectors(num_rates, k - 1, all_km1_rates_quant, &num_km1_vecs);
		prob_nr = 0;
		for (vec_nr = 0; vec_nr < num_km1_vecs; vec_nr++) {
			rates_quant = get_ith_vec(all_km1_rates_quant, vec_nr, num_rates);
			li = hypoexp_cumulate(rates_quant, k - 1, hypoexp);
			for (rate_nr = 0; rate_nr < num_rates; rate_nr++) {
				rates_quant[rate_nr] += 1;
				succ_vec_nr = vector_to_number
						(rates_quant, all_k_rates_quant, num_rates, 0, num_k_vecs-1);
				rates_quant[rate_nr] -= 1;
				rate_to_succ_vec_nr[rate_nr] = succ_vec_nr;
			}

			for (state_nr = 0; state_nr < num_states; state_nr++) {
				if (get_bit_val(B, state_nr)) {
					probs_km1[prob_nr] = li;
				} else {
					unsigned state_start = row_starts[state_nr];
					unsigned state_end = row_starts[state_nr + 1];
					double opt_choice_prob = min ? 2.0 : -1.0;
					for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
						unsigned choice_start = choice_starts[choice_nr];
						unsigned choice_end = choice_starts[choice_nr + 1];
						choice_prob = 0.0;
						rate_nr = choices_to_rates[choice_nr];
						succ_vec_nr = rate_to_succ_vec_nr[rate_nr];
						for (succ_nr = choice_start; succ_nr < choice_end; succ_nr++) {
							succ_prob = probs_k[succ_vec_nr * num_states + (ctmdp->cols)[succ_nr]];
							choice_prob += non_zeros[succ_nr] * succ_prob;
						}
						opt_choice_prob = optimum(min, choice_prob, opt_choice_prob);
					}
					if (state_start == state_end) {
						opt_choice_prob = 0.0;
					}
					probs_km1[prob_nr] = opt_choice_prob;
				}
				prob_nr++;
			}
		}
		swap_double(&probs_k, &probs_km1);
		swap_unsigned(&all_k_rates_quant, &all_km1_rates_quant);
		num_k_vecs = num_km1_vecs;
	}

	free_hypoexp(hypoexp);
	free(probs_km1);
	free(all_k_rates_quant);
	free(all_km1_rates_quant);
	free(rate_to_succ_vec_nr);

	return probs_k;
}

/**
 * Compute vector corresponding to reachability within zero time.
 *
 * @param B set of target states
 * @return double vector corresponding to reachability in zero time
 */
static double *time_zero_result(const bitset *B) {
	NDSparseMatrix *ctmdp;
	double *result;
	unsigned state;

	ctmdp = get_mdpi_state_space();
	result = malloc(ctmdp->n * sizeof(double));
	for (state = 0; state < (unsigned) ctmdp->n; state++) {
		result[state] = get_bit_val(B, state);
	}

	return result;
}

/**
 * Non-uniform history dependent CTMDPI time bounded reachability.
 * Algorithm taken from
 *
 * "Continuous-Time Stochastic Games with Time-Bounded Reachability"
 *
 * by Tomas Brazdil, Vojtech Forejt, Jan Krcal, Jan Kretinsky, Antonin
 * Kucera.
 *
 * @author Ernst Moritz Hahn (emh@cs.uni-sb.de)
 * @param time_bound time bound
 * @param B reachability set
 * @return computed probability vector
 */
double *transient_ctmdpi_hd_non_uni
(const double time_bound, const bitset *B)
{
	NDSparseMatrix *ctmdp;
	unsigned num_rates;
	double *rates;
	unsigned *choices_to_rates;
	unsigned k_bound;
	double *result;
	int comparator;
        BOOL min;

	if (0.0 == time_bound) {
		result = time_zero_result(B);
	} else {
		comparator = get_comparator();
		min = (C_GREATER == comparator) || (C_GREATER_EQUAL == comparator);

		ctmdp = get_mdpi_state_space();
		prepare_ctmdpi(ctmdp, &rates, &num_rates, &choices_to_rates);
		k_bound = compute_k_bound(rates, num_rates, time_bound);
		result = ctmdpi_iterate(ctmdp, rates, num_rates, choices_to_rates,
							k_bound, B, time_bound, min);
		unprepare_ctmdpi(ctmdp, rates, choices_to_rates);
	}

	return result;
}
