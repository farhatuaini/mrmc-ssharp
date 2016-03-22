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
*	Source description: Lumping of Markov chains.
*	Uses: DEF: lump.h, partition.h, splay.h, sparse.h, bitset.h,
*		label.h, runtime.h, macro.h,
*		LIB: partition.c, splay.c, sparse.c, bitset.c, label.c,
*/

#include "lump.h"
#include "sort.h"

#include "runtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
* Predeclarations:
*/
static /*@only@*/ /*@null@*/ sparse *calculate_lumped_probabilities(partition*P,
                /*@observer@*/ const sparse * Q)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@ensures notnull P->blocks@*/
                /*@modifies *P->blocks@*/;
static void quicksort(int * cols, double * vals, int left, int right)
                /*@modifies *cols, *vals@*/;

/**
* This internal subroutine finds all predecessors of splitter and treats them
* as follows:
* If some block C has to be regarded as absorbing or contains (at most) one
* state, it is not changed.
* Otherwise, if C contains predecessors, a new block C_pred is inserted near
* begin(C). C_pred is added to the list of predecessor classes, but not to the
* list of potential splitters.
* C_pred starts out empty. All predecessor states in C are then moved from C to
* C_pred. (At the end of the subroutine, C_pred will be nonempty, but C may have
* become empty.) However, the pointers block_of(P, s) are NOT changed,
* so that (temporarily) the data structure becomes inconsistent in this point.
* begin(C_pred) is the old value of begin(C); end(C) is not changed.
* For the states that are moved to some C_pred, the subroutine also calculates
* the probability to enter the splitter. This probability is saved in P->sum[s].
*/
static err_state find_predecessors(partition *P, /*@observer@*/ const sparse *Q,
                /*@observer@*/ const block * splitter)
                /*@requires notnull P->blocks, P->sum, P->pos@*/
                /*@requires isnull P->first_PredCl@*/
                /*@modifies P->first_PredCl, *P->sum, *P->pos, P->ib[],
                        P->num_blocks, splitter->next@*/
{
        pos_index s_pos, begin = part_block_begin_nt(splitter);

        for ( s_pos = part_block_end_nt(splitter) ; s_pos-- > begin ; ) {
                state_index s = P->ib[s_pos].id;

                mtx_walk_column(Q, t, (const int) s, val) {
                        /*@dependent@*/ block * C;
                        /*@dependent@*/ /*@null@*/ block * C_pred;
                        pos_index t_pos;

                        C = block_of(P, t);
                        /* if ( weak CTMC lumping && C == splitter )
                                continue; */
                        if ( 0 != (C->flags & ABSORBING) ) {
                                continue;
                        }

                        C_pred = C->next;
                        if ( NULL == C_pred || PARTITIONED != C_pred->flags ) {
                                /* C has not yet been split into C and C_pred */
                                if ( part_block_size_nt(C) <= 1 ) {
                                        continue;
                                }
                                if ( (C_pred = part_split_block(P, C,
                                                part_block_begin_nt(C)))==NULL )
                                {
                                        err_msg_5(err_CALLBY, "find_"
                                                "predecessors(%p,%p[%dx%d],%p)",
                                                (void *) P, (const void *) Q,
                                                mtx_rows(Q), mtx_cols(Q),
                                                (const void *) splitter,
                                                err_ERROR);
                                }
                                C_pred->flags = PARTITIONED;
                                C_pred->u.next_PredCl = P->first_PredCl;
                                P->first_PredCl = C_pred;
                        }
                        t_pos = P->pos[t];
                        if ( part_block_end_nt(C_pred) <= t_pos ) {
                                /* t is not yet in C_pred */
                                if ( /* splitter == C && */ part_block_end_nt(
                                                                C_pred) < s_pos
                                                        && s_pos <= t_pos )
                                {
                                        state_index ss;
                                        /* t is a state in the splitter; it has
                                           already been handled as predecessor
                                           of s; but if we would just swap it
                                           without thinking, then it would be
                                           handled again.
                                           Therefore we move t to position begin
                                           (guaranteed to be different states);
                                           the state at begin to end(C_pred)
                                           (may be the same); the state at
                                           end(C_pred) to s_pos (guaranteed to
                                           be different); the state at s_pos to
                                           t (may be the same).
                                           Note that after these swaps, s is no
                                           longer at s_pos. */
                                        if ( t_pos != s_pos ) {
                                                ss= P->ib[t_pos].id
                                                        = P->ib[s_pos].id;
                                                P->pos[ss] = t_pos;
                                        }
                                        ss = P->ib[s_pos].id
                                        = P->ib[part_block_end_nt(C_pred)].id;
                                        P->pos[ss] = s_pos++;
                                        if ( begin!=part_block_end_nt(C_pred) ){
                                                ss = P->ib[part_block_end_nt(
                                                C_pred)].id = P->ib[begin].id;
                                                P->pos[ss]
                                                = part_block_end_nt(C_pred);
                                        }
                                        P->ib[begin].id = t;
                                        P->pos[t] = begin++;
                                /* }else if ( part_block_end_nt(C_pred) >= s_pos
                                                        && s_pos > t_pos ) {
                                        That is impossible, as
                                        t_pos >= part_block_end_nt(C_pred). */
                                } else {
                                        /* swap state t with state at
                                           part_block_end_nt(C_pred) */
                                        state_index tt = P->ib[t_pos].id
                                        = P->ib[part_block_end_nt(C_pred)].id;
                                        P->pos[tt] = t_pos;
                                        P->ib[part_block_end_nt(C_pred)].id = t;
                                        P->pos[t] = part_block_end_nt(C_pred);
                                }
                                /* part_block_end_nt */ C_pred->end++;
                                P->sum[t] = 0.0;
                        }
                        P->sum[t] += val;
                } end_mtx_walk_column;
        }
        return err_OK;
}


/**
* This internal subroutine splits all blocks in P according to P->sum[...].
* It assumes that find_predecessors() has just finished.
* It also corrects the block pointers block_of(P, ...) that have been
* destroyed by find_predecessors().
*/
static err_state refine(partition * P)
                /*@requires notnull P->blocks, P->sum, P->pos@*/
                /*@ensures isnull P->first_PredCl@*/
                /*@modifies P->first_Sp, P->first_PredCl, *P->blocks,
                        P->num_blocks, P->ib, *P->pos@*/
{
        /*@dependent@*/ /*@null@*/ block * C_pred, * next_PredCl;

        C_pred = P->first_PredCl;
        P->first_PredCl = NULL;
        for ( ; NULL != C_pred ; C_pred = next_PredCl ) {
                /*@dependent @*/ block * C;
                BOOL C_is_splitter;
                pos_index C_begin = part_block_begin_nt(C_pred);

                next_PredCl = C_pred->u.next_PredCl;
                /* C_pred->u.next_PredCl = NULL; */
                /* C_pred->flags = 0; */
                /* find C */
                C = block_of(P, part_unlump_state_block_nt(P, C_pred));
                if ( C->next != C_pred || PARTITIONED != C_pred->flags
                                || 0 != (C->flags & ~ SPLITTER) )
                {
                        err_msg_1(err_INCONSISTENT, "refine(%p)", (void *) P,
                                        err_ERROR);
                } else if ( 0 == part_block_size_nt(C) ) {
                        /* C is empty: replace the contents of C by the contents
                           of C_pred and delete C_pred (C cannot be deleted
                           directly because there may be pointers to C.) */
                        /*@only@*/ /*@null@*/ block * temp = C_pred->next;
                        if ( NULL != C->next ) {
                                if ( err_state_iserror(free_block(C->next)) ) {
                                        /*@-mustfreeonly@*/
                                        err_msg_1(err_CALLBY, "refine("
                                                "%p)", (void *) P,
                                                err_ERROR);
                                        /*@=mustfreeonly@*/
                                }
                        }
                        C->next = temp;
                        /* C->end = C_pred->end; are already the same */
                        /* C->u.next_Sp has to be preserved */
                        /* C->flags has to be preserved */
                        P->num_blocks--;
                        C_pred = C;
                } else {
                        /* insert C_pred into list of splitters */
                        C_pred->flags = SPLITTER;
                        /* If C is not (yet) a splitter, the following
                           assignment is incomplete; that will be
                           corrected below */
                        C_pred->u.next_Sp = C->u.next_Sp;
                        C->u.next_Sp = C_pred;
                }
                /* insert C into list of splitters */
                if ( SPLITTER == C->flags ) {
                        C_is_splitter = TRUE;
                } else {
                        C->flags = SPLITTER;
                        /* Here, we have to use C_pred instead of C because of
                           the incomplete assignment just above */
                        if ( NULL == C_pred->u.next_Sp ) {
                                C_pred->u.next_Sp = P->first_Sp;
                                P->first_Sp = C;
                        }
                        C_is_splitter = FALSE;
                }

                /* split C_pred according to the probabilities in sum. Note that
                   the pointers block_of(P, ...) still point to the
                   block C. That will also be corrected by
                   sort_and_split_block. */
                if ( err_state_iserror(sort_and_split_block(P, C_pred, P->sum,
                                                        C != C_pred)) )
                {
                        err_msg_1(err_CALLBY,"refine(%p)",(void *) P,err_ERROR);
                }
                /* If C was not a splitter: remove the largest subblock of the
                   original C from the list of splitters. */
                if ( ! C_is_splitter ) {
                        /*@dependent@*/ /*@null@*/ block * temp;
                        /*@dependent@*/ block * largest = C;
                        state_count largest_size = part_block_size_nt(C);

                        for ( temp = C->next ; NULL != temp
                                        && part_block_end_nt(temp) > C_begin ;
                                        temp = temp->next )
                        {
                                if ( part_block_size_nt(temp) > largest_size ) {
                                        largest_size = part_block_size_nt(temp);
                                        largest = temp;
                                }
                        }
                        largest->flags = 0;
                }
        }
        return err_OK;
}


/*****************************************************************************
name		: lump
role		: This method computes the optimal lumped matrix of the input matrix.
@param		: partition *P: The initial partition.
@param		: sparse *Q: The matrix to lump.
@return     : The lumped matrix.
******************************************************************************/
/*@only@*/ /*@null@*/ sparse * lump(/*@i1@*/ /*@null@*/partition * P,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * Q)
{
        block * largest_b;
        int largest_size;
        sparse * res = NULL;

        if ( part_is_invalid(P) || NULL == P->pos || NULL != P->sum || NULL == Q
                        || part_unlumped_state_space_size_nt(P) != mtx_rows(Q))
        {
                if ( NULL != P ) {
                        free(P->pos);
                        P->pos = NULL;
                        free(P->sum);
                        P->sum = NULL;
                }
                /*@-nullstate@*/
                err_msg_4(err_PARAM, "lump(%p,%p[%dx%d])", (void *) P,
                                (const void *) Q, NULL != Q ? mtx_rows(Q) : 0,
                                NULL != Q ? mtx_cols(Q) : 0, NULL);
                /*@=nullstate@*/
        }

        P->first_PredCl = NULL;
        if ( (P->sum = calloc((size_t) part_unlumped_state_space_size_nt(P),
                                sizeof(double))) == NULL )
        {
                err_msg_4(err_MEMORY, "lump(%p,%p[%dx%d])", (void *) P,
                                (const void *) Q, mtx_rows(Q), mtx_cols(Q),
                                (free(P->pos), P->pos = NULL, free(P->sum),
                                P->sum = NULL, P->first_Sp = NULL, NULL));
        }

        /* first mark all blocks as splitters */
        largest_b = NULL;
        largest_size = 0;
        P->first_Sp = P->blocks;
        part_walk_blocks(P, b) {
                if ( part_block_size_nt(b) > largest_size ) {
                        largest_size = part_block_size_nt(b);
                        largest_b = b;
                }
                b->u.next_Sp = b->next;
                b->flags = SPLITTER | (b->flags & ABSORBING);
        } end_part_walk_blocks;
        if ( 0 != isRunMode(DTMC_MODE | DMRM_MODE | CTMDPI_MODE)
                                                && NULL != largest_b )
        {
                /* delete the largest block from the list of splitters. This is
                   allowed in DTMCs and in uniform CTMCs. */
                largest_b->flags &= ~ SPLITTER;
        }

        /* As long as there are (potential) splitters: */
        while ( NULL != P->first_Sp ) {
                block * splitter = P->first_Sp;

                if ( 0 != (splitter->flags & PARTITIONED) ) {
                        err_msg_4(err_INCONSISTENT, "lump(%p,%p[%dx%d])",
                                        (void*) P, (const void*) Q, mtx_rows(Q),
                                        mtx_cols(Q),
                                        (free(P->sum), P->sum = NULL,
                                        free(P->pos), P->pos = NULL,
                                        P->first_Sp = NULL, NULL));
                }
                P->first_Sp = splitter->u.next_Sp;
                splitter->u.next_Sp = NULL;
                if ( 0 != (splitter->flags & SPLITTER) ) {
                        splitter->flags &= ~ SPLITTER;
                        if(err_state_iserror(find_predecessors(P, Q, splitter))
                                        || err_state_iserror(refine(P)) )
                        {
                                err_msg_4(err_CALLBY, "lump(%p,%p[%dx%d])",
                                        (void*) P, (const void*) Q, mtx_rows(Q),
                                        mtx_cols(Q),
                                        (free(P->sum), P->sum = NULL,
                                        free(P->pos), P->pos = NULL,
                                        P->first_Sp = NULL, NULL));
                        }
                }
        }
        P->first_Sp = NULL; /* splint does not see that at this point,
                                P->first_Sp == NULL always. */
        free(P->sum);
        P->sum = NULL;
        free(P->pos);
        P->pos = NULL;

        /* assign numbers to blocks */
        if ( err_state_iserror(part_number_blocks(P))
/*                        || err_state_iserror(print_partition(P)) // DEBUG */
                        || (res = calculate_lumped_probabilities(P,Q)) == NULL )
        {
                err_msg_4(err_CALLBY, "lump(%p,%p[%dx%d])", (void *) P,
                                (const void*)Q, mtx_rows(Q), mtx_cols(Q), NULL);
        }
        return res;
}


/**
* The function allocates a new sparse matrix that contains the lumped transition
* matrix, based on the partition and the original transition matrix.
* The function assumes that the blocks in the partition have been duly numbered
* (by calling part_number_blocks()).
*/
static /*@only@*/ /*@null@*/ sparse *calculate_lumped_probabilities(partition*P,
                /*@observer@*/ const sparse * Q)
{
        state_count n;
        state_index xi, row;
        /*@only@*/ /*@null@*/ int * ncolse, * L1;
        /*@only@*/ /*@null@*/ sparse * Q1 = NULL;

        n = part_lumped_state_space_size_nt(P);
        /* The following code is taken mostly from the old lump() routine. */

	/* determine number of unique outgoing transitions for each block */
        L1 = (int *) calloc((size_t) n, sizeof(int));
        ncolse = (int *) calloc((size_t) n, sizeof(int));
        if ( NULL == L1 || NULL == ncolse ) {
                err_msg_4(err_MEMORY, "calculate_lumped_probabilities(%p,%p[%dx"
                                "%d])", (const void *) P, (const void *) Q,
                                mtx_rows(Q), mtx_cols(Q),
                                (free(ncolse), free(L1), NULL));
	}
	/* L1 is now used to store if a block is a successor */
	memset(L1, -1, n * sizeof(int));	/* init L1[i] with -1 for all i<n */
        part_walk_blocks(P, B) {
                if ( 0 != (B->flags & ~ ABSORBING) ) {
                        err_msg_4(err_INCONSISTENT, "calculate_lumped_"
                                        "probabilities(%p,%p[%dx%d])",
                                        (const void *) P, (const void *) Q,
                                        mtx_rows(Q), mtx_cols(Q),
                                        (free(ncolse), free(L1), NULL));
                }
                if ( 0 != test_flag((unsigned) B->flags, ABSORBING) ) {
			continue;
		}
                xi = part_unlump_state_block_nt(P, B);
                row = part_lump_state_block_nt(P, B);
                mtx_walk_row_nodiag(Q, (const int) xi, col, UNUSED(val)) {
                        if ( L1[part_lump_state_id_nt(P, col)] != row )
                        {
                                L1[part_lump_state_id_nt(P, col)] = row;
				++ncolse[row];
			}
		}
                end_mtx_walk_row_nodiag;
	}
        end_part_walk_blocks;
	free(L1);

	/* line 6-7 */
	Q1 = allocate_sparse_matrix_ncolse(n, n, ncolse);
	if(Q1 == NULL){
                err_msg_4(err_MEMORY, "calculate_lumped_probabilities(%p,%p[%dx"
                                "%d])", (const void *) P, (const void *) Q,
                                mtx_rows(Q), mtx_cols(Q), (free(ncolse), NULL));
	}
	free(ncolse);

	/* line 8-12 */
        part_walk_blocks(P, B) {
                row = part_lump_state_block_nt(P, B);
                if ( 0 != test_flag((unsigned) B->flags, ABSORBING) ) {
			/* make this block absorbing (only needed for formula-dependent lumping) */
                        if ( 0 != isRunMode(DTMC_MODE | DMRM_MODE) ) {
                                if ( err_state_iserror(mtx_set_diag_val(Q1, row,
                                                                        1.0)) )
                                {
                                        err_msg_4(err_CALLBY,"calculate_lumped_"
                                                "probabilities(%p,%p[%dx%d])",
                                                (const void*) P,(const void*) Q,
                                                mtx_rows(Q), mtx_cols(Q),
                                                (free_sparse_ncolse(Q1), NULL));
                                }
			}
			continue;
		}
		/* line 9 */
                xi = part_unlump_state_block_nt(P, B);
		/* line 10-12 */
                mtx_walk_row(Q, (const int) xi, col, val) {
                        if ( err_state_iserror(add_mtx_val_ncolse(Q1, row,
                                        part_lump_state_id_nt(P, col), val)) )
                        {
                                err_msg_4(err_CALLBY, "calculate_lumped_"
                                                "probabilities(%p,%p[%dx%d])",
                                                (const void*) P,(const void*) Q,
                                                mtx_rows(Q), mtx_cols(Q),
                                                (free_sparse_ncolse(Q1), NULL));
                        }
		}
                end_mtx_walk_row;
		/* order row by column index */
                if ( 1 < mtx_next_num(Q1, row) ) {
                        /*@-type@*/
                        /* sparse matrix internals used */
                        quicksort(Q1->valstruc[row].col, Q1->valstruc[row].val,
                                        0, mtx_next_num(Q1, row) - 1);
                        /*@=type@*/
		}
	}
        end_part_walk_blocks;
        printf("Lumping: The number of partition blocks is %d\n", mtx_rows(Q1));

	/* line 13 */
	return Q1;
}


/*****************************************************************************
name		: quicksort
role		: This method sorts the row elements by column index.
@param		: int *cols: column indices of non-zero elements.
@param		: int *vals: values of non-zero elements.
@param		: int left: first index.
@param		: int right: last index.
******************************************************************************/
static void quicksort(int *cols, double* vals, int left, int right)
{
	int i, j, x, y;
	double z;

	i = left; j = right;
	x = cols[(left+right)/2];

	do {
		while((cols[i] < x) && (i < right)) i++;
		while((x < cols[j]) && (j > left)) j--;

		if(i <= j) {
			y = cols[i];
			cols[i] = cols[j];
			cols[j] = y;

			z = vals[i];
			vals[i] = vals[j];
			vals[j] = z;

			i++; j--;
		}
	} while(i <= j);

	if(left < j) quicksort(cols, vals, left, j);
	if(i < right) quicksort(cols, vals, i, right);
}


/**
* This method changes the labelling structure such that it corresponds to the partition.
* @param	: labelling *labellin: the labelling structure.
* @param	: partition *P: the partition.
*/
err_state change_labelling(/*@i1@*/ /*@null@*/ labelling * labellin,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P)
{
        /*@only@*/ /*@null@*/ bitset_ptr * b = NULL;
        int label;
        BOOL all_clear;

        if ( NULL == labellin || part_is_invalid(P) )
                err_msg_2(err_PARAM, "change_labelling(%p,%p)", (void*)labellin,
                                (const void *) P, err_ERROR);

        b = (bitset_ptr *) calloc((size_t) labellin->n, sizeof(bitset_ptr));
        if ( NULL == b ) {
                err_msg_2(err_MEMORY,"change_labelling(%p,%p)", (void*)labellin,
                                (const void *) P, err_ERROR);
        }

	for(label = 0; label < labellin->n; label++){
                b[label] = lump_bitset(P, labellin->b[label]);
                if ( NULL == b[label] ) {
                        while ( label-- > 0 ) {
                                (void) /*@-nullpass@*/ free_bitset(b[label])
                                                        /*@=nullpass@*/;
                        }
                        err_msg_2(err_MEMORY, "change_labelling(%p,%p)",
                                        (void *) labellin, (const void *) P,
                                        (free(b), err_ERROR));
                }
	}

        all_clear = TRUE;
	for(label = 0; label < labellin->n; label++){
                if ( err_state_iserror(free_bitset(labellin->b[label])) ) {
                        all_clear = FALSE;
                }
	}
	free(labellin->b); /* Have to free the array itself also. */
        if ( ! all_clear ) {
                labellin->b = NULL;
                for ( label = 0 ; label < labellin->n ; label++ ) {
                        (void) free_bitset(b[label]);
                }
                /*@-nullstate@*/
                err_msg_2(err_CALLBY, "change_labelling(%p,%p)",
                                (void *) labellin, (const void *) P,
                                (free(b), err_ERROR));
                /*@=nullstate@*/
        }

	labellin->b = b;
        labellin->ns = part_lumped_state_space_size_nt(P);
        return err_OK;
}


/**
* This method changes the state rewards structure such that it
* corresponds to the partition. The actual change is done only
* in the case when the provided p_old_rew pointer is not NULL.
* NOTE:		We do not free the old_rewards structure because
*		it should be done on the outer level.
* @param	: double * p_old_rew: the the rewards structure.
* @param	: partition *P: the partition.
* @return	: The rewards structure after lumping or NULL if the original
*		  rewards structure is NULL
*/
/*@only@*/ /*@null@*/ double * change_state_rewards(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * p_old_rew,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P)
{
        /*@only@*/ /*@null@*/
	double * p_new_rew = NULL;

        if ( part_is_invalid(P) ) {
                err_msg_2(err_PARAM, "change_state_rewards(%p,%p)",
                                (const void*) p_old_rew, (const void*) P, NULL);
        }

	/* Test whether the state rewards are present. */
	if ( p_old_rew ){
                p_new_rew = (double *) calloc(
                                (size_t) part_lumped_state_space_size_nt(P),
                                sizeof(double));
                if ( NULL == p_new_rew ) {
                        err_msg_2(err_MEMORY, "change_state_rewards(%p,%p)",
                                        (const void *) p_old_rew,
                                        (const void *) P, NULL);
                }

		/* Take an arbitrary state of block, because all states */
		/* of the same block must have the same state reward */
                part_walk_blocks(P, pBlock) {
                        p_new_rew[part_lump_state_block_nt(P, pBlock)] =
                                p_old_rew[part_unlump_state_block_nt(P,pBlock)];
		}
                end_part_walk_blocks;
	}
	return p_new_rew;
}


/*****************************************************************************
name		: unlump_vector
role		: unlump a probability vector with respect to the partition.
@param		: partition *P: the partition of the original state space.
@param		: int size: the number of states in the original state space.
@param		: double *result_lumped: the probability vector of the lumped state space.
@return	: the probability vector of the original state space.
******************************************************************************/
/*@only@*/ /*@null@*/ double * unlump_vector(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P,
                state_count size,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * result_lumped)
{
	int i;
        /*@only@*/ /*@null@*/ double * result = NULL;

        if ( part_is_invalid(P) || part_unlumped_state_space_size_nt(P) != size
                        || NULL == result_lumped )
        {
                err_msg_3(err_PARAM, "unlump_vector(%p,%d,%p)", (const void*) P,
                                size, (const void *) result_lumped, NULL);
        }

        result = (double *) calloc((size_t) size, sizeof(double));
        if ( NULL == result ) {
                err_msg_3(err_MEMORY, "unlump_vector(%p,%d,%p)", (const void*)P,
                                size, (const void *) result_lumped, NULL);
        }

	for(i=0;i<size;i++){
                result[i] = result_lumped[part_lump_state_id_nt(P, i)];
	}

	return result;
}


/*****************************************************************************
name		: lump_bitset
role		: lump a bitset with respect to the partition.
@param		: partition *P: the partition of the original state space.
@param		: bitset *phi: the original bitset.
@return	: the lumped bitset.
******************************************************************************/
/*@only@*/ /*@null@*/ bitset * lump_bitset(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * phi)
{
        /*@only@*/ /*@null@*/ bitset * phi_lumped = NULL;

        if ( part_is_invalid(P) || NULL == phi
                                        || part_unlumped_state_space_size_nt(P)
                                                != bitset_size(phi) )
        {
                err_msg_2(err_PARAM, "lump_bitset(%p,%p)", (const void *) P,
                                (const void *) phi, NULL);
        }

        phi_lumped = get_new_bitset(part_lumped_state_space_size_nt(P));
        if ( NULL == phi_lumped ) {
                err_msg_2(err_MEMORY, "lump_bitset(%p,%p)", (const void *) P,
                                (const void *) phi, NULL);
        }
        part_walk_blocks(P, B) {
		/* if B is labeled with phi */
                if ( get_bit_val(phi, part_unlump_state_block_nt(P, B)) ) {
                        if ( err_state_iserror(set_bit_val(phi_lumped,
                                        part_lump_state_block_nt(P,B),BIT_ON)) )
                        {
                                err_msg_2(err_CALLBY, "lump_bitset(%p,%p)",
                                                (const void *) P,
                                                (const void *) phi,
                                                (free_bitset(phi_lumped),NULL));
                        }
		}
	}
        end_part_walk_blocks;

	return phi_lumped;
}
