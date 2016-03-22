/**
*    WARNING: Do Not Remove This Section
*
*       $LastChangedRevision: 415 $
*       $LastChangedDate: 2010-12-18 17:21:05 +0100 (Sa, 18. Dez 2010) $
*       $LastChangedBy: davidjansen $
*
*    MRMC is a model checker for discrete-time and continuous-time Markov
*    reward models. It supports reward extensions of PCTL and CSL (PRCTL
*    and CSRL), and allows for the automated verification of properties
*    concerning long-run and instantaneous rewards as well as cumulative
*    rewards.
*
*    Copyright (C) The University of Twente, 2004-2008.
*    Copyright (C) RWTH Aachen, 2008-2009.
*    Authors: Maneesh Khattri, Ivan Zapreev, Tim Kemna
*
*    This program is free software; you can redistribute it and/or
*    modify it under the terms of the GNU General Public License
*    as published by the Free Software Foundation; either version 2
*    of the License, or (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*    MA  02110-1301, USA.
*
*    Main contact:
*        Lehrstuhl für Informatik 2, RWTH Aachen University
*        Ahornstrasse 55, 52074 Aachen, Germany
*        E-mail: info@mrmc-tool.org
*
*    Old contact:
*        Formal Methods and Tools Group, University of Twente,
*        P.O. Box 217, 7500 AE Enschede, The Netherlands,
*        Phone: +31 53 4893767, Fax: +31 53 4893247,
*        E-mail: mrmc@cs.utwente.nl
*
*    Source description: Perform Transient Analysis for CSRL - X, U.
*    I. Discretization for CSRL - U.
*        Uses: DEF: bitset.h, sparse.h, label.h, runtime.h
*            LIB: bitset.c, sparse.c, label.c, runtime.c
*        Remarks: Discretization for CSRL - U.
*            SLOW setting FTSKN
*        For discussion see:
*          1. B.R. Haverkort, L. Cloth, H.Hermanns, J.-P. Katoen,
*             C. Baier. Model Checking Performability Properties.
*             DSN-2002, IEEE CS Press, pp. 103-112, 2002.
*          2. M. Khattri and R. Pulungan. Model Checking Markov
*             Reward Models with Impulse Rewards. M.S. Thesis.
*             University of Twente. 2004.
*          3. L. Cloth, J.-P. Katoen, M. Khattri, R. Pulungan,
*             Model Checking Markov Reward Models with Impulse
*             Rewards, DSN-PDS-05.
*
*        NOTE: This is too slow though. haven't had time to think
*            of something better.
*    II. Uniformization based algorithm by Bruno Sericola for CSRL - U.
*        Uses: DEF: bitset.h, sparse.h, label.h, runtime.h
*            LIB: bitset.c, sparse.c, label.c, runtime.c
*        Remarks: Uniformization based algorithm by Bruno Sericola
*            for CSRL - U.
*        For discussion see:
*          1. Bruno Sericola. Occupation Times in Markov Processes.
*                    Stochastic Models 16(5), 2000. pp. 479-510.
*          2. B.R. Haverkort, L. Cloth, H.Hermanns, J.-P. Katoen,
*             C. Baier. Model Checking Performability Properties.
*             DSN-2002, IEEE CS Press, pp. 103-112, 2002.
*          Only valid for MRMs without impulse rewards.
*        NOTE: INCOMPLETE - Needs to be implemented
*    III. Uniformization based algorithm by Qureshi&Sanders for CSRL
*        - U.
*        Uses: DEF: bitset.h, sparse.h, label.h, runtime.h
*            LIB: bitset.c, sparse.c, label.c, runtime.c
*        Remarks: Uniformization based algorithm by
*            Qureshi&Sanders for CSRL - U.
*
*        For discussion see:
*                    1. M.A. Qureshi and W.H.Sanders. Reward Model
*            Solution Methods with Impulse and Rate Rewards:
*            An Algorithm and Numerical Results.
*            Performance Evaluation, 20(4), 1994.
*                    2. M.A. Qureshi and W.H. Sanders. A new
*            methodology for calculating distributions of
*            reward accumulated during a finite interval.
*            Fault-Tolerant Computing Symposium, IEEE CS Press,
*            pp. 116-125, 1996.
*
*        FOR OMEGA SEE
*                    3. M.C. Diniz, E. de Souza e Silva and H.R. Gail.
*            Calculating the distribution of a linear
*            combination of uniform of order statistics.
*            INFORMS Journal on Computing: 14(2), pp. 124-131,
*            2002.
*
*        FOR Model Checking Markov Reward Models with Impulse
*        Rewards see:
*                    4. M. Khattri, R. Pulungan, Model Checking Markov
*            Reward Models with Impulse Rewards, M.Sc. Thesis,
*            UT-FMG-2004.
*                    5. L. Cloth, J.-P. Katoen, M. Khattri, R. Pulungan,
*            Model Checking Markov Reward Models with Impulse
*            Rewards, DSN-PDS-05.
*        NOTE: THIS DOES NOT SEEM TO WORK CORRECTLY OR MAY BE IS
*            NUMERICALLY UNSTABLE
*/

#include "transient_ctmrm.h"

#include "transient_common.h"
#include "kjstorage.h"
#include "lump.h"

#include "runtime.h"

#include <string.h>
#include <math.h>
#ifndef isnan
#       define isnan(x) FALSE
#endif
#ifndef NAN
#       define NAN (0.0/0.0)
#endif


/***********************************************************************************/
/********************************INTERVAL NEXT OPERATOR*****************************/
/***********************************************************************************/

/* Just some predeclaration: */
static double * interval_next_k(const bitset * phi, const double * subk,
                const double * supk);

/**
* Solve next_rewards formula.
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: double subi: lower rime bound
* @param: double supi: upper time bound
* @param: double subj: lower reward bound
* @param: double supj: upper reward bound
* @return: double *: result of the next formula for all states.
* NOTE: ********without impulse rewards at present***********
*          ********to be tested***********
*/
double * ctmrm_interval_next(const bitset * phi, double subi, double supi,
                double subj, double supj)
{
    int i, size = mtx_rows(get_state_space());
    double *subk = NULL, *supk = NULL;
        const
    double *reward = NULL;
    double *pResult = NULL;

    subk = (double *) calloc((size_t) size, sizeof(double));
    supk = (double *) calloc((size_t) size, sizeof(double));
    reward = getStateRewards();
    for( i=0; i < size; i++ ){
        if( reward[i] == 0 ){
            subk[i] = 0.0;
            supk[i] = 0.0;
        }else{
            subk[i] = subj / reward[i];
            supk[i] = supj / reward[i];
            if( subk[i] > supi || supk[i] < subi ){
                subk[i] = 0;
                supk[i] = 0;
            }else{
                if( subk[i] <= subi && supk[i] >= supi ){
                    subk[i] = subi;
                    supk[i] = supi;
                }else{
                    if( subk[i] >= subi && subk[i] <= supi && supk[i] >= supi ){
                        supk[i] = supi;
                    }else{
                        if( subk[i] <= subi && supk[i] <= supi && supk[i] >= subi ){
                            subk[i] = subi;
                        }else{
                            /* DOTO: WHAT IS THIS CASE FOR?! */
                            if( subk[i] >= subi && subk[i] <= supi && supk[i] <= supi && supk[i] >= subi){
                            }
                        }
                    }
                }
            }
        }
    }

    pResult = interval_next_k(phi, subk, supk);

    /* Free allocated memory */
    free( subk ); free( supk );

    return pResult;
}

/**
* Solve the interval next operator for CTMRMs (with K(s)).
* @param: bitset *phi: SAT(phi).
* @param: double subi: lower time bound
* @param: double supi: upper time bound
* @return: double *: result of the interval next operator for all states.
* NOTE: ********without impulse rewards at present***********
*/
static double * interval_next_k(const bitset * phi, const double * subk,
                const double * supk)
{
    const sparse *state_space = get_state_space();
    int size = mtx_rows(state_space), i;
    double * vec = (double *) calloc((size_t) size, sizeof(double)), vec_t;
    double * res = (double *) calloc((size_t) size, sizeof(double));
    const double *row_sums;

    for( i = 0; i < size; i++ ) {
        if( get_bit_val( phi, i ) ) {
            vec[i] = 1.0;
        }
    }
    if ( err_state_iserror( multiply_mtx_MV( state_space, vec, res ) )
         || (free(vec), vec = NULL, (row_sums = get_row_sums()) == NULL) )
    /* Note: DO NOT FREE the array you have acquired! */
    /* It is globally stored in runtime.c */
    {
        err_msg_4(err_CALLBY, "interval_next_k(%p[%d],%p,%p)",
                  (const void *) phi, bitset_size(phi),
                  (const void *) subk, (const void *) supk,
                  (free(res), free(vec), NULL));
    }
    for(i=0;i<size;i++) {
        if( (vec_t = row_sums[i]) != 0.0 ) {
            res[i]=(res[i]*(exp(-vec_t*subk[i])-exp(-vec_t*supk[i])))/vec_t;
        }
    }

    return res;
}

/***********************************************************************************/
/**************************TIME AND REWARD BOUNDED UNTIL****************************/
/***********************************************************************************/

/* Just some predeclarations: */
static double * uniformization_qureshi_sanders( const bitset * phi,
                                                const bitset * psi, 
                                                double supi, double supj,
                                                double ** ppResultError);
static double * discretization_wrapper( const bitset * phi, const bitset * psi,
                                        double supi, double supj);
static double discretization( int state_local, const bitset * phi,
                              const bitset * psi, const double m_time,
                              const double m_reward, const double d_factor);
static double* ctmrm_bounded_until_universal( const bitset * phi,
                                              const bitset * psi, 
                                              double supi, double supj,
                                              double ** ppResultError);

/**
* Solve bounded until formula with rewards.
* Does optimization: excludes states from phi_and_not_psi
* from which you always go to bad i.e. not Phi and not Psi states
* Also removes Phi & Psi states
* NOTE: The lower bound for time and reward is zero that is
* why the optimization is possible!
* @param phi the phi states.
* @param psi the psi states.
* @param supi the upper time bound
* @param supj the upper reward bound
* @param ppResultError  the pointer to the array of doubles. This array will
*                       store the errors
*            for the resulting probabilities. (This is the return variable)
* @return the result of the until formula.
*/
double * ctmrm_bounded_until( const bitset * phi, const bitset * psi,
                              double supi, double supj, double ** ppResultError)
{
    double *result = NULL;

    /* Do optimization, we exclude states from Phi which are Psi states and */
    /* also states from which you always go to bad i.e. not Phi ^ not Psi states */
    /* This is only possible if until has lower time bound 'subi' */
    /* (and reward bound 'subj' if any) equal to 0. */
    bitset *good_phi_states = get_good_phi_states( phi, psi, get_state_space() );

    result = ctmrm_bounded_until_universal( good_phi_states, psi, supi, supj, ppResultError );

    free_bitset(good_phi_states);

    return result;
}

/**
* Solve bounded until formula with rewards with lumping.
* NOTE: The lower bound for time and reward is zero
* @param phi the phi states.
* @param psi the psi states.
* @param supi the upper time bound
* @param supj the upper reward bound
* @param ppResultError  the pointer to the array of doubles. This array will
*                       store the errors
*            for the resulting probabilities. (This is the return variable)
* @return the result of the until formula for all states.
*/
double * ctmrm_bounded_until_lumping(const bitset * phi, const bitset * psi,
                double supi, double supj, double ** ppResultError)
{
    sparse *original_state_space = get_state_space();
    double *result, *lumped_result, * pLumpedResultError = NULL;

    /* Do optimization, we exclude states from Phi which are Psi states and */
    /* also states from which you always go to bad i.e. not Phi ^ not Psi states */
    /* This is only possible if until has lower time bound 'subi' */
    /* (and reward bound 'subj' if any) equal to 0. */
    bitset * good_phi_states = get_good_phi_states(phi, psi, original_state_space);

    /* Create initial partition for bounded until formula and lump */
    partition *P = init_partition_formula(good_phi_states, psi, FALSE);
    sparse *Q = lump(P, original_state_space);
    bitset *lumped_good_phi_states, *lumped_psi;
    double * p_old_rew, * p_new_rew;

    /* Free the row sums vector because otherwise it will be lost */
    free_row_sums();
    /* Set lumped CTMC to be the state space, NOTE: resets the row_sums! */
    set_state_space(Q);

    /* Lump the bitsets */
    lumped_good_phi_states = lump_bitset(P, good_phi_states);
    lumped_psi = lump_bitset(P, psi);

    /* Change state rewards */
    p_old_rew = getStateRewards();
    p_new_rew = change_state_rewards(p_old_rew, P);
    setStateRewards(p_new_rew);

    lumped_result = ctmrm_bounded_until(lumped_good_phi_states, lumped_psi, supi, supj, &pLumpedResultError);

    /* Unlump the resulting vector */
    result = unlump_vector(P, mtx_rows(original_state_space), lumped_result );

    /* The vector of errors is only allocated when the Uniformization Qureshi-Sanders is used */
    if( pLumpedResultError != NULL ){
        /* Unlump the result errors vector */
        *ppResultError = unlump_vector(P,mtx_rows(original_state_space),
                                       pLumpedResultError);
        /* Free the unneeded memory */
        free( pLumpedResultError ); pLumpedResultError = NULL;
    }

    /* Restore state rewards */
    freeStateRewards();
    setStateRewards(p_old_rew);

    /* Free allocated memory */
    free_bitset(good_phi_states);
    free_bitset(lumped_good_phi_states);
    free_bitset(lumped_psi);
    free( lumped_result ); lumped_result = NULL;
    free_partition(P);
    /* You have to do that not to waste the allocated memory */
    free_row_sums();

    /* Restore the original state space, NOTE: resets the row_sums! */
    set_state_space(original_state_space);

    /* Free the lumped state space. */
    if ( err_state_iserror(free_sparse_ncolse(Q)) ) {
        err_msg_7(err_CALLBY, "ctmrm_bounded_until_lumping(%p[%d],%p"
                  "[%d],%g,%g,%p)", (const void *) phi, bitset_size(phi),
                  (const void *) psi, bitset_size(psi), supi, supj,
                  (void *) ppResultError, (free(result), NULL));
    }

    return result;
}

/**
* Solve bounded until formula with rewards.
* @param phi the phi states.
* @param psi the psi states.
* @param supi the upper time bound
* @param supj the upper reward bound
* @param ppResultError  the pointer to the array of doubles. This array will
*                       store the errors
*            for the resulting probabilities. (This is the return variable)
* @return the result of the until formula.
*/
static double * ctmrm_bounded_until_universal( const bitset *phi,
                                               const bitset * psi, 
                                               double supi, double supj,
                                               double ** ppResultError)
{
    const int mur = get_method_until_rewards();
    double * result = NULL;

    if( mur == UQS ){
        printf("Performing Uniformization Qureshi-Sanders\n");
        result = uniformization_qureshi_sanders(phi, psi, supi, supj, ppResultError);
    } else {
        if( mur == US ){
            result = (double *) calloc((size_t) bitset_size(phi), sizeof(double));
            printf("ERROR: The uniformization based algorithm by Bruno Sericola for CSRL is not supported");
        } else if( mur == DTV ){
            printf("Performing Discretization Tijms-Veldman\n");
            result = discretization_wrapper( phi, psi, supi, supj);
        }
    }
    return result;
}

/***********************************************************************************/
/**********************************DISCRETIZATION***********************************/
/***********************************************************************************/

/**
* Test for 'd > 1/E(s)', as in this case we should terminate the algorithm
* I.e. 'd_factor > row_sums[i]' in the program turms
* @return: TRUE if d <= 1/E(s) for all states s otherwise FALSE
*/
static BOOL isDiscretizationFactorOk(void) {
    const double d_factor = get_d_factor();
    const double * row_sums = get_row_sums();
    const sparse * pStateSpace = get_state_space();
    int i = 0;

    /* Find the maximum and minimum row sum */
    double max_row_sum = -HUGE_VAL;
    BOOL isDiscretizationFactorOk_result;
    BOOL isFatal;
    double new_d_factor;

    for ( i = 0 ; i < mtx_rows(pStateSpace) ; i++ ) {
        if ( max_row_sum < row_sums[i] ){
            max_row_sum = row_sums[i];
        }
    }

    /* Test for 'd <= 1/E(s)', if this is NOT the case then we should terminate */
    isDiscretizationFactorOk_result = TRUE;
    isFatal = FALSE; /* Allows to indicate the exceptional case of all states being absorbing! */
    new_d_factor = 0.0;
    if( max_row_sum >= 1.0 ){
        if( max_row_sum != 0.0 ){
            /* Then we should get 'd <= 1 / max_row_sum' otherwise it is BAD */
            new_d_factor = 1.0 / max_row_sum ;
            if( d_factor > new_d_factor ){
                isDiscretizationFactorOk_result = FALSE;
            }
        }else{
            /* This is dead code, probably a bug. David N. Jansen.*/

            /* In case all row sums are equal to zero! */
            isFatal = TRUE;
            isDiscretizationFactorOk_result = FALSE;
        }
    }

    /* Print error message */
    if ( ! isDiscretizationFactorOk_result ) {
        if( isFatal ){
            printf("ERROR: All states turned out to be absorbing with ZERO exit rates, do not know how to proceed!\n" );
        }else{
            /* NOTE: we have '<=' below because the value gets truncated while printing */
            /* It is important that we print with %g here otherwise the value gets rounded */
                        /* and a user gets into trouble resetting the
                           discretization factor. */
            printf("ERROR: Reduce the discretization factor, it should be <= %g. See the 'set d R' command.\n", new_d_factor );
        }
    }

    return isDiscretizationFactorOk_result;
}

/**
* Solve t-r-bounded until formula in MRM using discretization.
* This method does error checks and calls for discretization procedure.
* The error checks are related for checking the validity of the discretization rate
* NOTE: if for some s 'd > 1/E(s)' then computations are not possible
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: bitset *psi: satisfaction relation for psi formula.
* @param: double supi: upper time bound
* @param: double supj: upper reward bound
* @return: double *: result of the until formula for all states.
*/
static double * discretization_wrapper(const bitset * phi, const bitset * psi,
                double supi, double supj)
{
    /* Allocate the resulting vector */
    double * result = (double *) calloc((size_t) bitset_size(phi), sizeof(double));

    if ( isDiscretizationFactorOk() ) {
        const int size = bitset_size(phi);
        BOOL result_error = FALSE;
        int i = 0;

        /* Select the 's|= not Phi ^ Psi' states */
        bitset * not_phi_and_psi;
        and_result( psi, not_phi_and_psi = not( phi ) );

        /* Compute probability for each state */
        for( i = 0; i < size; i++ ){
            if( get_bit_val( not_phi_and_psi , i ) ){
                /* NOTE: The case for 's|= not Phi ^ Psi and subi == 0 and subj == 0' */
                /* We do not support for subi>0, subj>0 case in MRMC so there is no check */
                result[i] = 1.0;
            }else{
                /* Perform discretization for the initial state i */
                result[i] = discretization(i, phi, psi, supi, supj, get_d_factor());

                /* Check the result for validity */
                if ( ( result[i] < 0.0 ) ||
                     ( result[i] > 1.0 + get_error_bound() ) ||
                     /*According to the specification NaN != NaN*/
                     ( result[i] != result[i] ) ) {
                    result[i] = NAN; /* Set to nan value */
                    result_error = TRUE;
                }
            }
        }

        /* Report on error if any */
        if( result_error ){
            printf("ERROR: Computations failure, try to reduce the discretization factor!\n");
            printf("ERROR: Beforehand, at least the probability values marked as 'nan' were invalid.\n");
        }

        free_bitset( not_phi_and_psi );
    }

    return result;
}

/**
* We have this procedure here just for convenience, to avoid duplications of the source code.
* 1. If it is the psi state then we just copy
*    the value into the F_t_next_s_k matrix
* 2. If it is the phi state then we increase
*    the reward by the idx1 state reward
* 3. If the discretized reward bound is not reached then
*    we store the updated value into the NEW reward element
*/
inline static void accumulateRewardForStaying( const BITSET_BLOCK_TYPE psi_state_bit, const int current_state,
                                               const int discrete_reward_state, const double prev_probability,
                                               const double *rew_array, const int R, const double *row_sums,
                                               const double d_factor, sparse *F_t_next_s_k){
    if( psi_state_bit ){
        /* In case this is the PSI state then we've reached the */
        /* desired state and we preserve the same probability */
        /* for this <state,reward> state as it is absorbing */
        if ( err_state_iserror(add_mtx_val(F_t_next_s_k, current_state,
             discrete_reward_state, prev_probability)) )
        {
            exit(err_macro_11(err_CALLBY,
                 "accumulateRewardForStaying(%d,%d,%d,%g,%p,%d,"
                 "%p,%g,%p[%dx%d])", psi_state_bit,current_state,
                 discrete_reward_state, prev_probability,
                 (const void *) rew_array, R,
                 (const void *) row_sums, d_factor,
                 (void *) F_t_next_s_k, mtx_rows(F_t_next_s_k),
                 mtx_cols(F_t_next_s_k), EXIT_FAILURE));
        }
    }else{
        int new_reward = discrete_reward_state + rew_array[current_state];
        /* Perhaps the following check should be changed to
           "new_reward <= R", similar to accumulateRewardForComming().
           David N. Jansen */
       if( new_reward < R ){
            if ( err_state_iserror(add_mtx_val(F_t_next_s_k,
                 current_state, new_reward,
                 prev_probability *
                 (1-row_sums[current_state]*d_factor))))
            {
                exit(err_macro_11(err_CALLBY,
                        "accumulateRewardForStaying(%d,%d,%d,"
                        "%g,%p,%d,%p,%g,%p[%dx%d])",
                        psi_state_bit, current_state,
                        discrete_reward_state, prev_probability,
                        (const void *) rew_array, R,
                        (const void *) row_sums, d_factor,
                        (void *) F_t_next_s_k,
                        mtx_rows(F_t_next_s_k),
                        mtx_cols(F_t_next_s_k), EXIT_FAILURE));
            }
        }
    }
}

/**
* We have this procedure here just for convenience, to avoid duplications of the source code.
* 1. Increase the reward by the prev_state state reward
*    and the impulse reward of going from prev_state to current_state
* 2. If the discretized reward bound is not reached
*    then we store the updated value into the
*    NEW reward element
*/
inline static void accumulateRewardForComming(const int current_state, const int prev_state, const int discrete_reward_state,
                        const double impulse_reward, const double prev_probability,
                        const double from_prev_to_current_state_rate, const double *rew_array,
                        const int R, const double d_factor, sparse *F_t_next_s_k){
    /* WARNING: Here according to the algorithm we should do impulse_reward / d_factor */
    /* which yields a non integer number, I am not sure that this is safe! */
    int new_reward = discrete_reward_state + rew_array[prev_state] + (int) ( impulse_reward / d_factor );
    if ( new_reward <= R ){
        if ( err_state_iserror(add_mtx_val(F_t_next_s_k, current_state,
             new_reward, prev_probability *
             from_prev_to_current_state_rate * d_factor)) )
        {
            exit(err_macro_12(err_CALLBY,
                    "accumulateRewardForComming(%d,%d,%d,%g,%g,%g,"
                    "%p,%d,%g,%p[%dx%d])", current_state,prev_state,
                    discrete_reward_state, impulse_reward,
                    prev_probability,
                    from_prev_to_current_state_rate,
                    (const void *) rew_array, R, d_factor,
                    (void *) F_t_next_s_k, mtx_rows(F_t_next_s_k),
                    mtx_cols(F_t_next_s_k), EXIT_FAILURE));
        }
    }
}

/**
* Solve t-r-bounded until formula in MRM using discretization.
* NOTE: The test for 'd > 1/E(s)' is done in the discretization_wrapper(...) method
* @param: int state: the starting state.
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: bitset *psi: satisfaction relation for psi formula.
* @param: double supi: upper time bound
* @param: double supj: upper reward bound
* @param: double d_factor: d factor.
* @param: BOOL testForError: TRUE if we are goinf to check for 'd > 1/E(s)'
* @param: BOOL *pError: Is TRUE if 'd > 1/E(s)', in this
*                 case we do not start computations at all.
* @return: double: result of the until formula for one state.
*/
static double discretization( int state_local, const bitset * phi,
                              const bitset * psi, const double m_time,
                              const double m_reward, const double d_factor)
{
    /* Shouldn't this be (int) (m_reward / d_factor) and
       (int) (m_time / d_factor), respectively? David N. Jansen. */
    const int R = (int) m_reward / d_factor, T = (int) m_time / d_factor,
              rows = bitset_size(phi);
    const sparse *pStateSpace = get_state_space();
    const sparse *rewi_matrix = getImpulseRewards();
    const double *row_sums = get_row_sums();
    const double *rew_array = getStateRewards();

    /* NOTE: We do not make states absording here in an explicit manner! */
    /* What is done is using the valid_states set and then checking */
    /* for PSI states. In PSI states all the rewards are ZERO!!!! */

    int current_state, prev_state;
    bitset *valid_states = or(phi, psi);

    sparse *F_t_now_s_k, *F_t_next_s_k;
    int t;
    double from_prev_to_current_state_rate, total_prob = 0.0,
           impulse_reward = 0.0;
    BITSET_BLOCK_TYPE psi_state_bit;

    /* Allocate the F_t_now_s_k matrix */
    /* correction of an MRMC bug: the matrix has to have columns 0...R, that
       makes R+1 columns. See also the second allocate_sparse_matrix() call
       a few lines below. David N. Jansen */
    F_t_now_s_k = allocate_sparse_matrix(rows, R + 1);  /* TODO: change into
            allocate_sparse_matrix_ncols, it will give some speedup */
    if ( NULL == F_t_now_s_k
        /* Iterate through time: */
        /* The initial value is for time t == 1 */
        || err_state_iserror(set_mtx_val(F_t_now_s_k,
                        state_local,(int)rew_array[state_local],
                        1 / d_factor)) )
    {
        err_msg_8(err_CALLBY, "discretization(%d,%p[%d],%p[%d],%g,%g,"
                "%g)", state_local, (const void*) phi, bitset_size(phi),
                (const void *) psi, bitset_size(psi), m_time, m_reward,
                d_factor,
                (NULL == F_t_now_s_k
                        || (free_mtx_sparse(F_t_now_s_k), FALSE),
                free_bitset(valid_states), NAN));
    }
    /* Compute F^{t}(s,k) */
    for( t = 2; t <= T; t++ ){
        /* Allocate the F_t_next_s_k matrix */
        F_t_next_s_k = allocate_sparse_matrix(rows, R + 1);
        if ( NULL == F_t_next_s_k ) {
                err_msg_8(err_CALLBY, "discretization(%d,%p[%d],%p[%d],"
                        "%g,%g,%g)", state_local, (const void*) phi,
                        bitset_size(phi), (const void *) psi,
                        bitset_size(psi),m_time,m_reward,d_factor,
                        (free_mtx_sparse(F_t_now_s_k),
                        free_bitset(valid_states), NAN));
        }

        /* Iterate through the valid states, i.row_sums. or(phi, psi) states */
        /* For F^{t}(current_state,k) */
        current_state = -1;
        while( ( current_state = get_idx_next_non_zero(valid_states, current_state) ) != -1 ){
            /* Iterate through the set rewards of the current_state row: */
            /* This is the first summand where we consider the reward */
            /* and probability of staying in the current state */
            psi_state_bit = get_bit_val(psi, current_state);
            mtx_walk_row(F_t_now_s_k, (const int) current_state,
                                    discrete_reward_state,
                                    prev_probability)
            {
                accumulateRewardForStaying(psi_state_bit, current_state, discrete_reward_state, prev_probability, rew_array, R, row_sums, d_factor, F_t_next_s_k);
            }
            end_mtx_walk_row;

            /* Iterate through phi and psi states again */
            /* This is the second summand, we go through the */
            /* states from which we can come to the current state */
            prev_state = -1;
            while( ( prev_state = get_idx_next_non_zero(valid_states, prev_state) ) != -1 ){
                psi_state_bit = get_bit_val( psi, prev_state );
                if ( err_state_iserror(get_mtx_val(pStateSpace,
                     prev_state, current_state,
                     &from_prev_to_current_state_rate
                    )) )
                {
                    err_msg_8(err_CALLBY, "discretization("
                            "%d,%p[%d],%p[%d],%g,%g,%g)",
                            state_local, (const void *) phi,
                            bitset_size(phi),
                            (const void *) psi,
                            bitset_size(psi), m_time,
                            m_reward, d_factor,
                            (free_mtx_sparse(F_t_next_s_k),
                            free_mtx_sparse(F_t_now_s_k),
                            free_bitset(valid_states),NAN));
                }

                if ( ! psi_state_bit && 0.0 != from_prev_to_current_state_rate)
                {
                    /* For the pure phi state prev_state from which there */
                    /* is a transition to state current_state */
                    if( rewi_matrix != NULL ){
                        /* Get the impulse reward for this transition */
                        if ( err_state_iserror( get_mtx_val( rewi_matrix,
                                                             prev_state,
                                                             current_state,
                                                             &impulse_reward ) ) )
                        {
                                err_msg_8( err_CALLBY, "discretization(%d,%p"
                                           "[%d],%p[%d],%g,%g,%g)",
                                           state_local, (const void *) phi,
                                           bitset_size(phi), (const void *) psi,
                                           bitset_size(psi), m_time, m_reward,
                                           d_factor, (free_mtx_sparse( F_t_next_s_k),
                                           free_mtx_sparse( F_t_now_s_k),free_bitset
                                           (valid_states), NAN));
                        }
                    }
                    /* Iterate through the prev_state row and: */
                    mtx_walk_row(F_t_now_s_k, (const int) prev_state, discrete_reward_state, prev_probability)
                    {
                        accumulateRewardForComming(current_state, prev_state, discrete_reward_state, impulse_reward, prev_probability,
                                        from_prev_to_current_state_rate, rew_array, R, d_factor, F_t_next_s_k);
                    }
                    end_mtx_walk_row;
                }
            }
        }

        /* Free the old matrix!!! */
        if ( err_state_iserror(free_mtx_sparse(F_t_now_s_k)) ) {
                err_msg_8(err_CALLBY, "discretization(%d,%p[%d],%p[%d],"
                        "%g,%g,%g)", state_local, (const void*) phi,
                        bitset_size(phi), (const void *) psi,
                        bitset_size(psi),m_time,m_reward,d_factor,
                        (free_mtx_sparse(F_t_next_s_k),
                        free_bitset(valid_states), NAN));
        }
        F_t_now_s_k = F_t_next_s_k;    /* SWAP the F_t_now_s_k and F_t_next_s_k matrices!!! */
        F_t_next_s_k = NULL;
    }

    /* Obtain the total probability by summing up resulting probabilities over all psi sates */
    current_state = -1;
    /* For all Psi states */
    while( ( current_state = get_idx_next_non_zero(psi, current_state) ) != -1 ){
        /* For all rewards from 1 to R */
        mtx_walk_row(F_t_now_s_k, (const int) current_state,
                        discrete_reward_state, prev_probability)
        {
            total_prob += prev_probability * d_factor;
        }
        end_mtx_walk_row;
    }
    /* Free allocated memory */
    free_bitset(valid_states);
    if ( err_state_iserror(free_mtx_sparse(F_t_now_s_k)) ) {
            err_msg_8(err_CALLBY, "discretization(%d,%p[%d],%p[%d],%g,%g,"
                    "%g)", state_local, (const void*) phi, bitset_size(phi),
                    (const void *) psi, bitset_size(psi), m_time, m_reward,
                    d_factor, NAN);
    }
    F_t_now_s_k = NULL;
    return total_prob;
}

/***********************************************************************************/
/*************************UNIFORMIZATION QURESHI SANDERS****************************/
/***********************************************************************************/

  /* GLOBAL */
static double maxPoi;
static int *index_;
static double *omIteration;
static double *Poisson;
static double *pathErr;
  /* GLOBAL */

/**
* Get the list of distinct state rewards from rew_array.(Descending order)
* @param: double *rew_array: vector containing the original state rewards.
* @param: bitset *valid_states: transient states after making states absorbing.
* @return: double *distinct_state_rewards_array: vector of distinct state rewards.
* NOTE: the output is a sorted list, so a translation of index(index_) is required.
*/
static double *get_dsr(const double *rew_array, const bitset *valid_states, const int rows, int *_p_number_of_distinct_state_rewards)
{
    int i, j;
    int number_of_distinct_state_rewards = 1; /* There is at least a zero reward */
        double * distinct_state_rewards_array = (double *) calloc((size_t) rows,
                        sizeof(double));
    double reward;
    int lower_index, upper_index, middle_index;

    /* Do the sorting of rewards for the valid states here */
    for( i=0; i < rows; i++ ){
        /* If we have a valid state */
        if( get_bit_val( valid_states, i ) ){
            BOOL found;

            reward = rew_array[i];
            lower_index = 0;
            upper_index = number_of_distinct_state_rewards - 1;
            middle_index = 0;
            found = FALSE;
            /* Search for the place where to add the reward (sorting) */
            while( lower_index <= upper_index ){
                middle_index = ( lower_index + upper_index ) / 2.0;
                if( reward == distinct_state_rewards_array[ middle_index ] ){
                    /* We found that this reward is in the array of considered rewards already */
                    found = TRUE;
                    index_[i] = middle_index;
                    break;
                }else {
                    /* If not then keep looking, the rewards are in a descending order */
                    if( reward > distinct_state_rewards_array[ middle_index ] ){
                        upper_index = middle_index - 1;
                    }else{
                        if( reward < distinct_state_rewards_array[ middle_index ] ){
                            lower_index = middle_index + 1;
                        }
                    }
                }
            }
            /* If the reward was not found */
            if( ! found ){
                if(  ( middle_index == ( number_of_distinct_state_rewards - 1)
                       && ( reward > distinct_state_rewards_array[ middle_index ] ) )
                    || middle_index < ( number_of_distinct_state_rewards - 1 ) ){
                    /* Add an extra reward - i.e. consider a new array element */
                    number_of_distinct_state_rewards ++;
                    /* Shift all rewards to get space for the new */
                    /* reward in distinct_state_rewards_array[ middle_index ] */
                    for( j = ( number_of_distinct_state_rewards - 1 ); j > middle_index; j-- ){
                        distinct_state_rewards_array[j] = distinct_state_rewards_array[j-1];
                    }
                                        /* Shift all the "state --> distinct
                                           reward array indices" mapping */
                    for( j = 0; j < rows; j++ ){
                        if( index_[j] >= middle_index ){
                            index_[j]++;
                        }
                    }
                    /* Put the new reward into the array and update mapping */
                    distinct_state_rewards_array[ middle_index ] = reward;
                    index_[i] = middle_index;
                } else {
                    /* We have to add the reward to the very end of the array. */
                    distinct_state_rewards_array[ number_of_distinct_state_rewards ++ ] = reward;
                    index_[i] = number_of_distinct_state_rewards;
                }
            }
        }
    }

    /* Reallocate the array, its size could be smaller than the number of states */
    distinct_state_rewards_array = (double *)realloc(distinct_state_rewards_array, number_of_distinct_state_rewards * sizeof(double));
    /* Update the output number of distinct rewards */
    * _p_number_of_distinct_state_rewards = number_of_distinct_state_rewards;

    return distinct_state_rewards_array;
}

/**
* All indices less than pivot are in G,
* pivot and greater indices are in L.
*/
static void omIt(int number_of_distinct_state_rewards, const int * tmp_k,
                double r_dash, int j, int pivot, int cols,
                const double * distinct_state_rewards_array)
{
    double prevValue = 1.0, val_i, val_j;
    int * tempG = (int *) calloc((size_t) number_of_distinct_state_rewards, sizeof(int)), count, i;

    for(count=0; count<cols;count++) {
        for(i=0;i<pivot;i++)
        {
            if(tempG[i]<tmp_k[i])
            {
                tempG[i]++;
                val_i = distinct_state_rewards_array[i];
                val_j = distinct_state_rewards_array[j];
                omIteration[count] = ((val_i-r_dash)/(val_i-val_j))*omIteration[count]+((r_dash-val_j)/(val_i-val_j))*prevValue;
                prevValue = omIteration[count];
                break;
            }
        }
    }
    free(tempG);
    tempG=NULL;
}

double omega(const int * tmp_k, double r_dash, int pivot,
                int number_of_distinct_state_rewards,
                const double * distinct_state_rewards_array)
{
    int cols = 0, rows=0, count, *tempL, j;
    for(count=pivot;count<number_of_distinct_state_rewards;count++) {
        rows+=tmp_k[count];
    }
    for(count=0;count<pivot;count++) {
        cols+=tmp_k[count];
    }
    if(rows!=0 && cols==0) {
        return 1.0;
    } else { 
        if (rows==0 && cols!=0) {
            return 0.0;
        }
    }
    tempL = (int *) calloc((size_t) number_of_distinct_state_rewards,
                    sizeof(int));
    omIteration = (double *) calloc((size_t) cols, sizeof(double));
    for(count=0; count<rows;count++)
    {
        for(j=pivot;j<number_of_distinct_state_rewards;j++)
        {
            if(tempL[j]<tmp_k[j])
            {
                tempL[j]++;
                omIt(number_of_distinct_state_rewards, tmp_k, r_dash, j, pivot, cols, distinct_state_rewards_array);
                break;
            }
        }
    }
    free(tempL); tempL=NULL;
    return omIteration[cols-1];
}

/**
* findPivot
* double r_dash, double * distinct_state_rewards_array, int number_of_distinct_state_rewards
* NOTE: all indices less than pivot are in G
*      pivot and all indices greater than pivot are in L
*      Ref: Solution for Order Statistics by Diniz, de Souza e Silva & H.R. Gail
*      Since di = (a(i)+a(i+1)...+a(n))
*      In our problem a(i) = r(i) - r(i+1)
*      therefore
*        di = ( (r(i)-r(i+1)) + (r(i+1)-r(i+2)) + ...+ (r(n)-r(n+1)) )
*      and since r(n+1) = 0; (Provided there is always absorbing state)
*      hence
*        di = r(i)
*      since all r(i)'s are distinct
*        ci = di = ri
*      all indices less than pivot are in G
*      pivot and all indices greater than pivot are in L
*/
static int findPivot(double r_dash, const double *distinct_state_rewards_array, int number_of_distinct_state_rewards)
{
    int i=0;
    for(;i<number_of_distinct_state_rewards;i++) {
        if(distinct_state_rewards_array[i]<=r_dash) {
            break;
        }
    }
    return i;
}

/**
* maxPoisson
* int n
* calculate poisson probability.
*/
static void maxPoisson(const int n, const double lambda, const double mTime, int *psize)
{
    int i;
    if(*psize==0)
    {
        maxPoi = exp(-(lambda*mTime));
        (*psize)++;
        Poisson[0]=maxPoi;
        pathErr[0]=maxPoi;
        /*printf("psize = %d, maxPoi=%lf\n", *psize, maxPoi);*/
    }
    if(*psize<=n)
    {
        Poisson = (double *)realloc(Poisson, (n+1)*sizeof(double));
        pathErr = (double *)realloc(pathErr, (n+1)*sizeof(double));
        for(i=*psize;i<=n;i++)
        {
            maxPoi *= (lambda*mTime)/(i);
            (*psize)++;
            Poisson[i] = maxPoi;
            pathErr[i] = pathErr[i-1]+maxPoi;
            /*printf("Poisson[%d]=%lf, psize = %d, maxPoi=%lf\n", i, Poisson[i], *psize, maxPoi);*/
        }
    }
}

/**
* errorBounds
* int n, double prob
* NOTE: calculate error bounds cumulatively given one discarded path at a time.
*/
static void errorBounds(const int n, const double prob, double *tot_err, const double lambda, const double mTime, int *psize)
{
    maxPoisson(n, lambda, mTime, psize);
    if(n>0) {
        *tot_err += prob * (1-pathErr[n-1]);
    } else {
        *tot_err += prob;
    }
}

/* static void testMaxPoisson(const double lambda, const double supi, int *psize)
{
    int n=100, i;
    maxPoisson(n, lambda, supi, psize);
    for(i=0;i<n;i++)
        printf("Poisson[%d] = %f, pathErr = %f\n", i, Poisson[i], pathErr[i]);
} */

/* static */
/* int recursion_level = 0; */

/**
* @param: w - Probability Threshold
* @param: double supi: upper time bound
* @param: double supj: upper reward bound
* @return: The algorithm computes the probability of
*        all paths in the uniformized DTMC (D). If we multiply the probability of
*        some path (in D) with the Poisson probability for a path of length 'n',
*        given that n is the length of the path in question, then we obtain the
*        probability of the path to be in the CTMC from which you get (D). But these
*        path probabilities are generated one-by-one given (D) which can be of any
*        length. So, it is necessary to stop the computation at some point. This is
*        where "w" - probability threshold - comes in, only paths whose probability
*        is more than 'w' are considered in the computation.
* NOTE: depth first path generation
*/
static void dfpg(int * global_k, const sparse * abs_local,
                const sparse * rewi_matrix, int n, int state_local,
        double p, double r_impulse, double w, double lambda, double supi, double supj,
        int *psize, double *tot_err, kjstruc *kjstruct, const bitset *psi)
{
    double pr;

    /* printf("[DEBUG] ---- RECURSION LEVEL: %d IN ----\n",
                    recursion_level); */
    maxPoisson(n, lambda, supi, psize);
    pr = p*Poisson[n];
    /*printf("pr=%lf\n",pr);*/
    /*TODO:
     * Where is the check for s |= not Phi /\ not Psi ?
     * Why do we check for supj < r_impulse here?
     */
    if( supj < r_impulse || pr < w )
    {
        errorBounds(n, pr, tot_err, lambda, supi, psize);
    }else{
        if ( get_bit_val(psi, state_local) )
        {
            add_new_kjnode(kjstruct, global_k, r_impulse, p);
            /*TODO: Why is it ELSE here? Looks like there sould be no else
            just the following part from the ELSE clause!!!*/
        } else {
            double val=0.0, value;
            mtx_walk_row_nodiag(abs_local, (const int) state_local, col, vall)
            {
                int kIndex;

                kIndex = index_[col];
                global_k[kIndex]++;
                value=0.0;
                if(rewi_matrix!=NULL) {
                    if ( err_state_iserror(get_mtx_val( rewi_matrix,state_local, col, &value)) )
                    {
                            exit(err_macro_20(err_CALLBY,
                                    "dfpg(%p,%p[%dx%d],%p["
                                    "%dx%d],%d,%d,%g,%g,%g,"
                                    "%g,%g,%g,%p,%p,%p,%p["
                                    "%d])", (void*)global_k,
                                    (const void*) abs_local,
                                    mtx_rows(abs_local),
                                    mtx_cols(abs_local),
                                    (const void*)rewi_matrix
                                    , mtx_rows(rewi_matrix),
                                    mtx_cols(rewi_matrix),n,
                                    state_local,p,r_impulse,
                                    w, lambda, supi, supj,
                                    (void *) psize,
                                    (void *) tot_err,
                                    (void *) kjstruct,
                                    (const void *) psi,
                                    bitset_size(psi),
                                    EXIT_FAILURE));
                    }
                }

                /* recursion_level++; */
                /*TODO: WHAT ABOUT THE STATE REWARDS HERE?!*/
                dfpg(global_k, abs_local, rewi_matrix, n+1, col,
                                p * vall, r_impulse + value, w,
                                lambda, supi, supj, psize,
                                tot_err, kjstruct, psi);
                /* recursion_level--; */

                global_k[kIndex]--;
            }
            end_mtx_walk_row_nodiag;
            if ( err_state_iserror(get_mtx_diag_val(abs_local, state_local, &val)))
            {
                exit(err_macro_20(err_CALLBY, "dfpg(%p,%p[%dx"
                        "%d],%p[%dx%d],%d,%d,%g,%g,%g,%g,%g,%g,"
                        "%p,%p,%p,%p[%d])", (void *) global_k,
                        (const void *) abs_local,
                        mtx_rows(abs_local),mtx_cols(abs_local),
                        (const void *) rewi_matrix,
                        mtx_rows(rewi_matrix),
                        mtx_cols(rewi_matrix), n, state_local,p,
                        r_impulse, w, lambda, supi, supj,
                        (void *) psize, (void *) tot_err,
                        (void *) kjstruct, (const void *) psi,
                        bitset_size(psi), EXIT_FAILURE));
            }
            if ( 0.0 != val )
            {
                int kIndex = index_[state_local];
                global_k[kIndex]++;
                value=0.0;
                
                /* recursion_level++; */
                dfpg(global_k, abs_local, rewi_matrix, n + 1,
                        state_local, p * val, r_impulse + value,
                        w, lambda, supi, supj, psize, tot_err,
                        kjstruct, psi);
                /* recursion_level--; */

                global_k[kIndex]--;
            }
        }
    }
    /* printf("[DEBUG] ---- RECURSION LEVEL: %d OUT ----\n",
                    recursion_level); */
}

/**
* Solve t-r-bounded until formula (phi U[0, supi][0, supj] psi) in MRM without impulse rewards
* using uniformization by qureshi & sanders.
* @param phi the phi states.
* @param psi the psi states.
* @param supi the upper time bound.
* @param supj the upper reward bound.
* @param ppResultError  the pointer to the array of doubles. This array will
*                       store the errors
*            for the resulting probabilities. It is going to be
*            allocated here and will be returned using ppResultError.
* @return: the resulting probabilities for the until formula.
*/
static double * uniformization_qureshi_sanders(const bitset * phi,
                const bitset * psi, double supi, double supj,
                double ** ppResultError)
{
    const int rows = bitset_size(phi);
    const sparse *rewi_matrix = getImpulseRewards();
    const double *rew_array = getStateRewards();
    const double w=get_w();

    bitset * valid_states;
    int number_of_distinct_state_rewards;
    double *distinct_state_rewards_array;
    int ki, non_absorbing, psize;
    double lambda, total_prob, tot_err;
    kjstruc *kjstruct;
    double * diag;
    sparse * abs_local;
    /* Stores the resulting probabilities */
    double *result;
    /* Stores the errors of the resulting probabilities */
    double *pResultError;
    int * global_k;

    valid_states = not(psi);
    and_result(phi, valid_states);  /* Remains Constant after this */

    /* index_ is modified in get_dsr(...) */
    index_ = (int *) calloc((size_t) rows, sizeof(int));
    distinct_state_rewards_array = get_dsr(rew_array, valid_states, rows, &number_of_distinct_state_rewards);    /* Remains Constant after this */
    /* print_vec_double(rows, distinct_state_rewards_array); */

    ki=0; non_absorbing=0; psize=0;
    lambda = 0.0; total_prob = 0.0; tot_err = 0.0;
    kjstruct=NULL;

    diag = (double *) calloc((size_t) rows, sizeof(double));
    abs_local = ab_state_space(get_state_space(), valid_states, &lambda,
                    diag, &non_absorbing);
    /* print_mtx_sparse(abs_local); */

    result = (double *) calloc((size_t) rows, sizeof(double));
    pResultError = (double *) calloc((size_t) rows, sizeof(double));
    global_k = (int *) calloc((size_t) number_of_distinct_state_rewards,
                    sizeof(int));
    Poisson = (double *) calloc((size_t) 1, sizeof(double));
    pathErr = (double *) calloc((size_t) 1, sizeof(double));

/*    printf("[DEBUG] Starting uniformization Qureshi-sanders.\n"); */

    /* R-E(s) */
    if ( err_state_iserror(sub_mtx_diagonal(abs_local, diag)) ) {
            err_msg_7(err_CALLBY, "uniformization_qureshi_sanders(%p[%d],"
                    "%p[%d],%g,%g,%p)", (const void*) phi, bitset_size(phi),
                    (const void *) psi, bitset_size(psi), supi, supj,
                    (void *) ppResultError,
                    (free_abs(abs_local), free(pathErr), free(Poisson),
                    free(global_k), free(pResultError), free(result),
                    free(diag), free(distinct_state_rewards_array),
                    free(index_), free_bitset(valid_states), NULL));
    }

    /* uniformize : lambda>0*/
    if ( 0.0 != lambda )
    {
        if ( err_state_iserror(mult_mtx_const(abs_local, 1 / lambda))
            || err_state_iserror(add_mtx_cons_diagonal( abs_local, 1.0)) )
        {
            err_msg_7(err_CALLBY, "uniformization_qureshi_sanders("
                    "%p[%d],%p[%d],%g,%g,%p)", (const void *) phi,
                    bitset_size(phi), (const void *) psi,
                    bitset_size(psi), supi, supj,
                    (void *) ppResultError,
                    (free_abs(abs_local), free(pathErr),
                    free(Poisson),free(global_k),free(pResultError),
                    free(result), free(diag),
                    free(distinct_state_rewards_array),free(index_),
                    free_bitset(valid_states), NULL));
        }
        /* print_mtx_sparse(get_sparse()); */

        /*print_mtx_sparse(abs_local);
        printf("lambda = %f, w = %f\n", lambda, w);*/
        for ( ki = 0 ; ki < rows ; ki++ )
        {
            int kIndex;
            kjnode *trav;

            total_prob = 0.0; tot_err = 0.0;
            kjstruct = get_new_kjstruc(number_of_distinct_state_rewards);
            /* Store the state reward index in the array of distinct state rewards */
            kIndex = index_[ki];

            /* Clean the recursion index array */
            memset(global_k, '\0', number_of_distinct_state_rewards * sizeof(int));

            /* print_mtx_sparse(abs_local); */
            global_k[kIndex]++;
            dfpg(global_k, abs_local, rewi_matrix, 0, ki, 1.0,
                            0.0, w, lambda, supi, supj, &psize,
                            &tot_err, kjstruct, psi);
            global_k[kIndex]--;     /* This does not do a thing !!!!! */

            /* testMaxPoisson(lambda, supi, &psize); */
            /* print_kjstorage(kjstruct); */
            trav = kjstruct->kjnodes;
            while(trav!=NULL)
            {
                double r_dash = supj/supi-trav->acir/supi;
                int pivot = findPivot(r_dash, distinct_state_rewards_array, number_of_distinct_state_rewards);
                int i, n=0;
                int * tmp_k=trav->k;
                double pr;
                double om;
                int cols;

                for(i=0; i<number_of_distinct_state_rewards; i++)
                    n+=tmp_k[i];
                --n;
                maxPoisson(n, lambda, supi, &psize);
                pr = trav->prob *  Poisson[n];
                om = omega(tmp_k, r_dash, pivot, number_of_distinct_state_rewards, distinct_state_rewards_array);
                cols = 0;
                for(i=0; i<pivot; i++) {
                    cols+=tmp_k[i];
                }
                while(pr>=w && r_dash>=0.0)
                {
                    total_prob += pr * (om);
                    tmp_k[number_of_distinct_state_rewards-1]++;
                    pr = pr * ((lambda*supi)/++n);
                    if(pr>=w && cols>0)
                    {
                        omIt(number_of_distinct_state_rewards, tmp_k, r_dash, number_of_distinct_state_rewards-1, pivot, cols, distinct_state_rewards_array);
                        om = omIteration[cols-1];
                    }
                }
                errorBounds(n, trav->prob, &tot_err, lambda, supi, &psize);
                trav=trav->next;
                free(omIteration); omIteration=NULL;
            }
            free_kjstorage(kjstruct);
            kjstruct=NULL;
            /* NOTE: indexes of states from the user point of view start from 1 */
            /* printf("State = %d, Probabilty = %e, Error Bound = %e \n", ki+1, total_prob, tot_err); */
            if( tot_err >= 0.1 ) printf("WARNING: The Probability Threshold w = %e is to big!\n", tot_err);
            result[ki] = total_prob;
            pResultError[ki] = tot_err;
        }

        /*Reset the matrix to its original state
        NOTE: operations on diagonals are not required */
        if ( err_state_iserror(mult_mtx_const(abs_local, lambda)) ) {
            err_msg_7(err_CALLBY, "uniformization_qureshi_sanders("
                    "%p[%d],%p[%d],%g,%g,%p)", (const void *) phi,
                    bitset_size(phi), (const void *) psi,
                    bitset_size(psi), supi, supj,
                    (void *) ppResultError,
                    (free_abs(abs_local), free(pathErr),
                    free(Poisson),free(global_k),free(pResultError),
                    free(result), free(diag),
                    free(distinct_state_rewards_array),free(index_),
                    free_bitset(valid_states), NULL));
        }
    }else{
        printf("ERROR: Uniformization rate is 0 in Qureshi Sanders\n");
    }

    /*free allocated memory*/
    free(distinct_state_rewards_array);
    free(diag);
    free_bitset(valid_states);
    free(Poisson); Poisson = NULL;
    free(pathErr); pathErr = NULL;
    free(global_k); global_k = NULL;
    free(index_); index_ = NULL;
    if ( err_state_iserror(free_abs(abs_local)) ) {
        err_msg_7(err_CALLBY, "uniformization_qureshi_sanders(%p[%d],"
                "%p[%d],%g,%g,%p)", (const void*) phi, bitset_size(phi),
                (const void *) psi, bitset_size(psi), supi, supj,
                (void *) ppResultError,
                (free(pResultError), free(result), NULL));
    }

    /* Assign the return argument value*/
    /* WARNING: Do not free pResultError, we copy a pointer here! */
    *ppResultError = pResultError;

    return result;
}

/***********************************************************************************/
/******************************UNIFORMIZATION SERICOLA******************************/
/***********************************************************************************/

/**
* Solve t-r-bounded until formula in MRM without impulse rewards
*  using uniformization.
* @param: int state: the starting state.
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: bitset *psi: satisfaction relation for psi formula.
* @param: double supi: upper time bound
* @param: double supj: upper reward bound
* @param: int: print.
* @return: double: result of the until formula for one state.
* NOTE: Impulse rewards are not supported by the algorithm.
*/
/* static double uniformization_sericola(int state, bitset *phi, bitset *psi, double supi, double supj, int print){
    sparse *r_matrix = get_state_space();
        const double *row_sums = get_row_sums();        // Note: DO NOT FREE the
                                                // array you have acquired!
                        // It is globally stored in runtime.c
    double *rew_matrix = getStateRewards();
    double eps=get_error_bound(), u=get_underflow(), o=get_overflow(), lambda=0.0, total_prob=0.0;
        int rows = bitset_size(phi), idx1, idx2;

    // TODO: IMPLEMENT

    return total_prob;
} */
