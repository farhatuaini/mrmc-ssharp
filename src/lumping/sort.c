/**
*	WARNING: Do Not Remove This Section
*
*       $LastChangedRevision: 256 $
*       $LastChangedDate: 2010-01-08 15:32:21 +0100 (Fri, 08 Jan 2010) $
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
*       Author: David N. Jansen
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
*       Source description: Sort partition structure for lumping.
*       Uses: DEF: partition.h
*/

/*@+looploopbreak@*/

#include "sort.h"

#include "runtime.h"

/**
* What is the maximum difference between two elements to put them into the same
* equivalence class?
*/
#define EPSILON 3E-12

/**
* The function finds the median of 3 entries.
* Parameters: P = partition
*       key = array of double values containing, for each state, a key
*       a, b, c = three positions in P->id
* Result: a, b, or c; depending on which state has the middle key.
*/
static pos_index median3_noerror(/*@observer@*/ /*@temp@*/ const partition * P,
                /*@observer@*/ const double * key, pos_index a, pos_index b,
                pos_index c)
                /*@requires maxRead(key) >= a /\ maxRead(key) >= b
                        /\ maxRead(key) >= c@*/
                /*@modifies nothing@*/
{
        if( /*@-realcompare@*/ key[P->ib[a].id]==key[P->ib[b].id]
                                                        /*@=realcompare@*/ )
        {
                return a;
        }
        if ( key[P->ib[a].id] > key[P->ib[b].id] ) {
                pos_index temp = a;
                a = b;
                b = temp;
        }
        /* a < b */
        if ( key[P->ib[b].id] <= key[P->ib[c].id] ) {
                /* a < b <= c */
                return b;
        }
        /* a < b and b > c */
        if ( key[P->ib[a].id] <= key[P->ib[c].id] ) {
                /* a <= c < b */
                return c;
        }
        /* c < a < b */
        return a;
}


/**
* The function finds a pivot candidate. For smaller arrays, it takes a median
* of 3 elements, and for larger ones, it takes a pseudomedian of 9 elements.
* It moves the pivot element to the last state of the block, but does not adapt
* P->pos[...].
* It is taken from
*       Bentley, Jon L.; McIlroy, M. Douglas: Engineering a sort function.
*       Software-practice and experience 23(11)1993, 1249-1265.
*/
static void find_pivot_noerror(partition * P, block * B,
                /*@observer@*/ const double * key)
                /*@requires notnull P->blocks@*/ /*@modifies P->ib[].id@*/;

static void find_pivot_noerror(partition * P, block * B,
                /*@observer@*/ const double * key)
{
        pos_index pivot_start, pivot_middle, pivot_end, temp;

        pivot_start = part_block_begin_nt(B);
        pivot_end = part_block_end_nt(B) - 1;
        pivot_middle = (pivot_start + pivot_end) / 2;
        if ( pivot_end - pivot_start >= 7 ) {
                if ( pivot_end - pivot_start >= 40 ) {
                        pos_index step = (pivot_end - pivot_start) / 8;
                        pivot_end = median3_noerror(P, key, pivot_end - 2*step,
                                        pivot_end - step, pivot_end);
                        pivot_middle = median3_noerror(P, key,pivot_middle-step,
                                        pivot_middle, pivot_middle + step);
                        pivot_start = median3_noerror(P, key, pivot_start,
                                        pivot_start + step,
                                        pivot_start + 2 * step);
                }
                pivot_middle = median3_noerror(P, key, pivot_start,pivot_middle,
                                pivot_end);
                pivot_end = part_block_end_nt(B) - 1;
        }
        /* Now, pivot_middle contains the position of the pivot, and pivot_end
           its desired position. */
        temp = P->ib[pivot_middle].id;
        P->ib[pivot_middle].id = P->ib[pivot_end].id;
        P->ib[pivot_end].id = temp;
}


/**
* The function makes a quicksort-like pass through the block. It is assumed that
* the last element contains the pivot. The function moves states with
* sum < key[pivot] towards the beginning, states with sum > key[pivot]
* to the middle, and states with sum == key[pivot] to the end of the block. In
* *middle, the position of the first state > pivot is returned, and the function
* result is the position of the first state == pivot.
* Please note that this function does not adapt P->pos[...].
*/
static state_index pass_file_noerror(partition * P, block * B,
                /*@observer@*/ const double *key, /*@out@*/ state_index *middle)
                /*@requires notnull P->blocks@*/
                /*@modifies *middle, P->ib[].id@*/
{
        pos_index end, small, large;
        state_index id_small, id_large;
        double pivot_key_minus_epsilon, pivot_key_plus_epsilon, cur_key;

        small = part_block_begin_nt(B);
        large = part_block_end_nt(B) - 1;
        pivot_key_minus_epsilon = key[P->ib[large].id] - EPSILON;
        pivot_key_plus_epsilon = key[P->ib[large].id] + EPSILON;
        /* begin(B) == small <= large == end(B) - 1
           key[id[large]] == pivot
           key[id[large+1]] ... key[id[end(B)-1]] is empty */
        /* First iteration special: as large is at the end of the file, we do
           not have to swap as long as we find only == pivot elements. */
        do {
                /* begin(B) == small <= large < end(B)
                   key[id[large]] ... key[id[end(B)-1]] == pivot */
                if ( small == large ) {
                        *middle = small;
                        return small;
                }
                /* begin(B) == small < large < end(B)
                   key[id[large]] ... key[id[end(B)-1]] == pivot */
                --large;
                id_large = P->ib[large].id;
                /* begin(B) == small <= large < end(B)
                   id_large == id[large]
                   key[id[large+1]] ... key[id[end(B)-1]] == pivot */
        } while ( (cur_key = key[id_large]) <= pivot_key_plus_epsilon
                        && cur_key >= pivot_key_minus_epsilon );
        /* begin(B) == small <= large < end(B)
           id_large == id[large]
           key[id_large] != pivot
           key[id[large+1]] ... key[id[end(B)-1]] == pivot */
        if ( ! (cur_key <= pivot_key_plus_epsilon) ) {
                /* begin(B) == small <= large < end(B)
                   key[id[large]] > pivot
                   key[id[large+1]] ... key[id[end(B)-1]] == pivot */
                end = large + 1;
                /* begin(B) == small <= large < end <= end(B)
                   key[id[large]] > pivot
                   key[id[large+1]] ... key[id[end-1]] is empty
                   key[id[end]] ... key[id[end(B)-1]] == pivot */
                /* Note: Exactly the same loop (with slightly different
                   assertions) is also placed towards the end of the function.
                   Let's hope that the compiler can optimize this code. */
                do {
                        do {
                                /* begin(B) == small <= large <= end <= end(B)
                                   key[id[large]] ... key[id[end-1]] > pivot
                                   key[id[end]]...key[id[end(B)-1]] == pivot */
                                if ( small == large ) {
                                        *middle = small;
                                        return end;
                                }
                                /* begin(B) == small < large <= end <= end(B)
                                   key[id[large]] ... key[id[end-1]] > pivot
                                   key[id[end]]...key[id[end(B)-1]] == pivot */
                                --large;
                                id_large = P->ib[large].id;
                                /* begin(B) == small <= large < end <= end(B)
                                   id_large == id[large]
                                   key[id[large+1]] ... key[id[end-1]] > pivot
                                   key[id[end]]...key[id[end(B)-1]] == pivot */
                        } while ( (cur_key = key[id_large])
                                                > pivot_key_plus_epsilon );
                        /* begin(B) == small <= large < end <= end(B)
                           id_large == id[large]
                           cur_key == key[id_large] <= pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        if ( cur_key < pivot_key_minus_epsilon ) {
                                break;
                        }
                        /* begin(B) == small <= large < end <= end(B)
                           id_large == id[large]
                           key[id_large] == pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        --end;
                        /* begin(B) == small <= large <= end < end(B)
                           id_large == id[large]
                           key[id_large] == pivot
                           key[id[large+1]] ... key[id[end]] > pivot
                           key[id[end+1]] ... key[id[end(B)-1]] == pivot */
                        /* Swap id[large] with id[end] */
                        P->ib[large].id = P->ib[end].id;
                        P->ib[end].id = id_large;
                        /* begin(B) == small <= large <= end < end(B)
                           key[id[large]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                } while ( TRUE );
                /* begin(B) == small <= large < end <= end(B)
                   id_large == id[large]
                   key[id_large] < pivot
                   key[id[large+1]] ... key[id[end-1]] > pivot
                   key[id[end]] ... key[id[end(B)-1]] == pivot */
        } else {
                end = large + 1;
                /* begin(B) == small <= large < end <= end(B)
                   id_large == id[large]
                   key[id_large] < pivot
                   key[id[large+1]] ... key[id[end-1]] is empty
                   key[id[end]] ... key[id[end(B)-1]] == pivot */
        }
        do {
                do {
                        /* begin(B) <= small <= large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small-1]] < pivot
                           id_large == id[large]
                           key[id_large] < pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        if ( small == large ) {
                                *middle = small + 1;
                                return end;
                        }
                        id_small = P->ib[small].id;
                        /* begin(B) <= small < large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small-1]] < pivot
                           id_small == id[small]
                           id_large == id[large]
                           key[id_large] < pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        if( (cur_key=key[id_small]) >= pivot_key_minus_epsilon )
                        {
                                /*@innerbreak@*/ break;
                        }
                        /* begin(B) <= small < large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small]] < pivot
                           id_large == id[large]
                           key[id_large] < pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        ++small;
                        /* begin(B) < small <= large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small-1]] < pivot
                           id_large == id[large]
                           key[id_large] < pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                } while ( TRUE );
                /* begin(B) <= small < large < end <= end(B)
                   key[id[begin(B)]] ... key[id[small-1]] < pivot
                   id_small == id[small]
                   cur_key == key[id_small] >= pivot
                   id_large == id[large]
                   key[id_large] < pivot
                   key[id[large+1]] ... key[id[end-1]] > pivot
                   key[id[end]] ... key[id[end(B)-1]] == pivot */
                P->ib[small].id = id_large;
                /* begin(B) <= small < large < end <= end(B)
                   key[id[begin(B)]] ... key[id[small]] < pivot
                   cur_key == key[id_small] >= pivot
                   id[large] has to be overwritten
                   key[id[large+1]] ... key[id[end-1]] > pivot
                   key[id[end]] ... key[id[end(B)-1]] == pivot */
                if ( cur_key <= pivot_key_plus_epsilon ) {
                        /* begin(B) <= small < large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small]] < pivot
                           key[id_small] == pivot
                           id[large] has to be overwritten
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        --end;
                        /* begin(B) <= small < large <= end < end(B)
                           key[id[begin(B)]] ... key[id[small]] < pivot
                           key[id_small] == pivot
                           id[large] has to be overwritten
                           key[id[large+1]] ... key[id[end]] > pivot
                           key[id[end+1]] ... key[id[end(B)-1]] == pivot */
                        P->ib[large].id = P->ib[end].id;
                        /* begin(B) <= small < large <= end < end(B)
                           key[id[begin(B)]] ... key[id[small]] < pivot
                           key[id_small] == pivot
                           key[id[large]] ... key[id[end-1]] > pivot
                           id[end] has to be overwritten
                           key[id[end+1]] ... key[id[end(B)-1]] == pivot */
                        P->ib[end].id = id_small;
                        /* begin(B) <= small < large <= end < end(B)
                           key[id[begin(B)]] ... key[id[small]] < pivot
                           key[id[large]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                } else {
                        /* begin(B) <= small < large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small]] < pivot
                           key[id_small] > pivot
                           id[large] has to be overwritten
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        P->ib[large].id = id_small;
                        /* begin(B) <= small < large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small]] < pivot
                           key[id[large]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                }
                /* begin(B) <= small < large <= end <= end(B)
                   key[id[begin(B)]] ... key[id[small]] < pivot
                   key[id[large]] ... key[id[end-1]] > pivot
                   key[id[end]] ... key[id[end(B)-1]] == pivot */
                ++small;
                /* begin(B) < small <= large <= end <= end(B)
                   key[id[begin(B)]] ... key[id[small-1]] < pivot
                   key[id[large]] ... key[id[end-1]] > pivot
                   key[id[end]] ... key[id[end(B)-1]] == pivot */

                /* Note: Exactly the same loop (with slightly different
                   assertions) is also placed above. Let's hope that the
                   compiler can optimize this code. */
                do {
                        do {
                                /* begin(B) < small <= large <= end <= end(B)
                                   key[id[begin(B)]]...key[id[small-1]] < pivot
                                   key[id[large]] ... key[id[end-1]] > pivot
                                   key[id[end]]...key[id[end(B)-1]] == pivot */
                                if ( small == large ) {
                                        *middle = small;
                                        return end;
                                }
                                /* begin(B) < small < large <= end <= end(B)
                                   key[id[begin(B)]]...key[id[small-1]] < pivot
                                   key[id[large]] ... key[id[end-1]] > pivot
                                   key[id[end]]...key[id[end(B)-1]] == pivot */
                                --large;
                                id_large = P->ib[large].id;
                                /* begin(B) < small <= large < end <= end(B)
                                   key[id[begin(B)]]...key[id[small-1]] < pivot
                                   id_large == id[large]
                                   key[id[large+1]] ... key[id[end-1]] > pivot
                                   key[id[end]]...key[id[end(B)-1]] == pivot */
                        } while ( (cur_key = key[id_large])
                                                > pivot_key_plus_epsilon );
                        /* begin(B) < small <= large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small-1]] < pivot
                           id_large == id[large]
                           cur_key == key[id_large] <= pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        if ( cur_key < pivot_key_minus_epsilon ) {
                                /*@innerbreak@*/ break;
                        }
                        /* begin(B) < small <= large < end <= end(B)
                           key[id[begin(B)]] ... key[id[small-1]] < pivot
                           id_large == id[large]
                           key[id_large] == pivot
                           key[id[large+1]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                        --end;
                        /* begin(B) < small <= large <= end < end(B)
                           key[id[begin(B)]] ... key[id[small-1]] < pivot
                           id_large == id[large]
                           key[id_large] == pivot
                           key[id[large+1]] ... key[id[end]] > pivot
                           key[id[end+1]] ... key[id[end(B)-1]] == pivot */
                        P->ib[large].id = P->ib[end].id;
                        P->ib[end].id = id_large;
                        /* begin(B) < small <= large <= end < end(B)
                           key[id[begin(B)]] ... key[id[small-1]] < pivot
                           key[id[large]] ... key[id[end-1]] > pivot
                           key[id[end]] ... key[id[end(B)-1]] == pivot */
                } while ( TRUE );
                /* begin(B) < small <= large < end <= end(B)
                   key[id[begin(B)]] ... key[id[small-1]] < pivot
                   id_large == id[large]
                   key[id_large] <= pivot
                   key[id[large+1]] ... key[id[end-1]] > pivot
                   key[id[end]] ... key[id[end(B)-1]] == pivot */
        } while ( TRUE );
        /*@notreached@*/
}


/**
* The function refines a block according to the floating point numbers in key[].
* It splits the block into subsets; each subset consists of states that have
* (almost) the same key value.
* The "_internal" function assumes that block B is a splitter and inserts all
* newly generated blocks into the list of splitters after B.
* The arrays P->ib[...].b and P->pos[...] are not changed by this
* function. Because it does swap elements and split B, they will have to be
* corrected afterwards.
* Parameters: P = partition
*       B = block to be refined;
*       key = array of doubles (key[s] = value for state s).
* Result: err_OK if everything went fine; err_ERROR otherwise.
*/
static err_state sort_and_split_block_internal(partition * P, block * B,
                /*@observer@*/ const double * key)
                /*@requires notnull P->blocks@*/
                /*@modifies P->num_blocks, P->ib[], B->u.next_Sp, B->next@*/
{
        /*@dependent@*/ /*@null@*/ block * B_left = NULL, * B_right = NULL;
        pos_index middle, equal;

        /* 1. qpartition the block. The "sorting" algorithm used here is similar
           to quicksort; with the difference that equal-to-pivot elements are
           collected at the end of the block B. */
        find_pivot_noerror(P, B, key);
        equal = pass_file_noerror(P, B, key, &middle);

        /* 2. split the block into up to three parts, and if desired, correct
              the pointers. Note that always equal < part_block_end_nt(B), as
              there will be at least one element equal to the pivot in B. */
        if ( middle > part_block_begin_nt(B) ) {
                if ( (B_left = part_split_block(P, B, middle)) == NULL ) {
                        err_msg_5(err_CALLBY, "sort_and_split_block_internal("
                                        "%p,%p[%d..%d),%p)", (void*)P, (void*)B,
                                        part_block_begin_nt(B),
                                        part_block_end_nt(B), (const void *)key,
                                        err_ERROR);
                }
                B_left->flags = SPLITTER;
                B_left->u.next_Sp = B->u.next_Sp;
                B->u.next_Sp = B_left;
                if ( 1 == middle - part_block_begin_nt(B_left) /* i. e.
                                                1 == part_block_size(B_left)*/ )
                {
                        B_left = NULL;
                }
        }
        if ( equal > middle ) {
                if ( (B_right = part_split_block(P, B, equal)) == NULL ) {
                        err_msg_5(err_CALLBY, "sort_and_split_block_internal("
                                        "%p,%p[%d..%d),%p)", (void*)P, (void*)B,
                                        part_block_begin_nt(B),
                                        part_block_end_nt(B), (const void *)key,
                                        err_ERROR);
                }
                B_right->flags = SPLITTER;
                B_right->u.next_Sp = B->u.next_Sp;
                B->u.next_Sp = B_right;
                if ( 1 == equal - middle /* i.e. 1==part_block_size(B_right)*/ )
                {
                        B_right = NULL;
                }
        }
        /* B_left now contains the elements smaller than the pivot; B_right
           contains the elements larger than the pivot; B contains the elements
           equal to pivot. */

        /* 3. qpartition the "smaller-than-pivot" and the "larger-than-pivot"
              parts. */
        /* Here, I suppress most of the callstack messages to allow the compiler
           to replace tail recursion by iteration. */
        if ( NULL == B_left ) {
                if ( NULL == B_right ) {
                        return err_OK;
                }
                return sort_and_split_block_internal(P, B_right, key);
        }
        if ( NULL == B_right ) {
                return sort_and_split_block_internal(P, B_left, key);
        }
        if ( equal - middle < middle - part_block_begin_nt(B_left) ) {
                /*@dependent@*/ block * temp = B_right;
                B_right = B_left;
                B_left = temp;
        }
        if ( err_state_iserror(sort_and_split_block_internal(P, B_left, key)) ){
                err_msg_5(err_CALLBY, "sort_and_split_block_internal(%p,%p[%d.."
                                "%d),%p)", (void *) P, (void *) B,
                                part_block_begin_nt(B), part_block_end_nt(B),
                                (const void *) key, err_ERROR);
        }
        return sort_and_split_block_internal(P, B_right, key);
}



/**
* The function is a simple wrapper around sort_and_split_block_internal() that
* makes the block B a splitter and checks parameters before calling the real
* subroutine to have the hard work done.
* Precondition: block B is not empty, and it is a splitter.
*/
err_state sort_and_split_block(/*@i1@*/ /*@null@*/ partition * P,
                /*@dependent@*/ /*@i1@*/ /*@null@*/ block * B,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * key,
                BOOL set_b)
{
        pos_index begin, i;

        if ( part_is_invalid(P) || NULL==B || SPLITTER!=B->flags || NULL==key ){
                err_msg_6(err_PARAM,"sort_and_split_block(%p,%p[%d..%d),%p,%d)",
                                (void *) P, (void *) B,
                                NULL == B ? 0 : part_block_begin_nt(B),
                                NULL == B ? 0 : part_block_end_nt(B),
                                (const void *) key, (int) set_b, err_ERROR);
        }

        begin = part_block_begin_nt(B);
        if ( part_block_size_nt(B) > 1
                        && err_state_iserror(sort_and_split_block_internal(P, B,
                                                                        key)) )
        {
                err_msg_6(err_CALLBY, "sort_and_split_block(%p,%p[%d..%d),%p,"
                                "%d)", (void *) P, (void *) B,
                                part_block_begin_nt(B), part_block_end_nt(B),
                                (const void*) key, (int) set_b, err_ERROR);
        }

        /* correct P->ib[...].b and P->pos[...] */
        i = part_block_end_nt(B);
        if ( ! set_b ) {
                pos_index new_begin = part_block_begin_nt(B);
                do {
                        --i;
                        P->pos[P->ib[i].id] = i;
                } while ( i > new_begin );
                if ( i <= begin ) {
                        return err_OK;
                }
                if ( (B = B->next) == NULL ) {
                        err_msg_3(err_INCONSISTENT,
                                        "sort_and_split_block(%p,...,%p,%d)",
                                        (void *) P, (const void *) key,
                                        (int) set_b, err_ERROR);
                }
        }
        do {
                pos_index new_begin = part_block_begin_nt(B);
                do {
                        --i;
                        P->ib[P->ib[i].id].b = B;
                        P->pos[P->ib[i].id] = i;
                } while ( i > new_begin );
                if ( i <= begin ) {
                        break;
                }
                if ( (B = B->next) == NULL ) {
                        err_msg_3(err_INCONSISTENT,
                                        "sort_and_split_block(%p,...,%p,%d)",
                                        (void *) P, (const void *) key,
                                        (int) set_b, err_ERROR);
                }
        } while ( TRUE );
        return err_OK;
}
