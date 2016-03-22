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
* 	Copyright (C) The University of Twente, 2004-2008.
* 	Copyright (C) RWTH Aachen, 2008-2009.
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
*	Source description: Manage partition structure for lumping.
*	Uses: DEF: partition.h
*/

#include "partition.h"
#include "sort.h"

#include "runtime.h"

#include <errno.h>
#include <string.h>

/*****************************************************************************
* The function splits a block into two parts. It creates a new block and moves
* the states in ib[begin(B) ... pos-1].id to the new block. The new block may
* be empty, but it is an error if the old block becomes empty. It also is an
* error if the block was absorbing.
* The function does not change the block pointers P->ib[...].b.
* Parameters:   P = partition
*               B_old = block
*               pos = position where the block has to be split.
* Result: a pointer to the new block if everything went fine; NULL if some error
* has happened.
******************************************************************************/
/*@dependent@*/ /*@null@*/ block * part_split_block(
                /*@i1@*/ /*@null@*/ partition * P,
                /*@i1@*/ /*@null@*/ block * B_old,
                pos_index pos)
{
	block *b;

        if ( part_is_invalid(P) || NULL == B_old
                        || 0 != (B_old->flags & ABSORBING)
                        || pos < part_block_begin_nt(B_old)
                        || pos >= part_block_end_nt(B_old) )
        {
                err_msg_3(err_PARAM, "part_split_block(%p,%p,%d)", (void *) P,
                                (void *) B_old, pos, NULL);
        }

        b = malloc(sizeof(block));
	if(b == NULL){
                err_msg_3(err_MEMORY, "part_split_block(%p,%p,%d)", (void *) P,
                                (void *) B_old, pos, NULL);
	}
        b->next = B_old->next;
        b->u.next_Sp = NULL;
        b->end = pos;
	b->flags = 0;
        B_old->next = b;

        P->num_blocks /* lumped size */++;
	return b;
}


/**
* This method prints the position of the bits in the original
* bitset whose value is 1.
* @param P: The partition.
* @param b: The bitset corresponding to the lumped state space.
*/
err_state print_original_states(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * b)
{
        int ns;
	int xj;
	BOOL bFirst = TRUE;

        if ( part_is_invalid(P) || NULL == b )
                err_msg_2(err_PARAM, "print_original_states(%p,%p)",
                                (const void *) P, (const void *) b, err_ERROR);

        errno = 0;

        ns = part_unlumped_state_space_size_nt(P);

		printf("{ ");
		for(xj = 0; xj < ns; xj++){
                        if ( get_bit_val(b, part_lump_state_id_nt(P, xj)) ) {
				if( ! bFirst ){
					printf(", ");
				}else{
					bFirst = FALSE;
				}
				printf("%d", xj+1);
			}
		}
		printf(" }");

        if ( 0 != errno ) {
                err_msg_2(err_FILE, "print_original_states(%p,%p)",
                                (const void *) P, (const void *) b, err_ERROR);
        }
        return err_OK;
}

/**
* This method prints TRUE if index is in the unlumped bitset 'b',
* otherwise false.
* The error messages must be printed exactly as shown because the caller does
* not check the parameters b or index.
* @param P: The partition.
* @param b: The bitset corresponding to the lumped state space.
* @param const: The original-state-space state index
*/
err_state print_original_state(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * b,
                state_index index_local)
{
        state_count original_state_space_size;

        if ( part_is_invalid(P) || NULL == b ) {
                err_msg_4(err_PARAM, "print_original_state(%p[%p],%p,%d)",
                                (const void *) P,
                                NULL == P ? NULL : (void *) P->blocks,
                                (const void *) b, index_local, err_ERROR);
        }

        errno = 0;
        original_state_space_size = part_unlumped_state_space_size_nt(P);

                if ( 0<=index_local && index_local<original_state_space_size ) {
                        if ( get_bit_val(b, part_lump_state_id_nt(P,
                                                                index_local)) )
                        {
				printf("TRUE");
			}else{
				printf("FALSE");
			}
		}else{
			printf("??\n");
			printf("WARNING: Invalid index %d, required to be in the [1, %d] interval.",
                                index_local + 1, original_state_space_size);
		}

        if ( 0 != errno ) {
                err_msg_3(err_FILE, "print_original_state(%p,%p,%d)",
                                (const void *) P, (const void *) b, index_local,
                                err_ERROR);
        }
        return err_OK;
}

/**
* Print the vector of probability errors.
* @param P the partition.
* @param pErrorValues the array error values.
* NOTE: We assume the the size of pValues equal to the number of partition blocks
*/
err_state print_error_probs_partition(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pErrorValues)
{
	int j, i;
        int original_state_space_size;

        if ( part_is_invalid(P) || NULL == pErrorValues ) {
                err_msg_2(err_PARAM, "print_error_probs_partition(%p,%p)",
                                (const void *) P, (const void *) pErrorValues,
                                err_ERROR);
        }

        errno = 0;
        original_state_space_size = part_unlumped_state_space_size_nt(P);

		printf("( ");
		for( j = 0 ; j < original_state_space_size; j++ ){
                        if ( 0 != j ) {
				printf(", ");
			}

                        i = part_lump_state_id_nt(P, j);
			printf( "%1.2e", pErrorValues[i] );
		}
		printf(" )");

        if ( 0 != errno ) {
                err_msg_2(err_FILE, "print_error_probs_partition(%p,%p)",
                                (const void *) P, (const void *) pErrorValues,
                                err_ERROR);
        }
        return err_OK;
}

/**
* Print the unlumped vector of probabilities.
* @param P the partition.
* @param vec the array of probabilities to be printed.
* @param error_bound the error bound for all probabilities
* @param pErrorBounds the error bounds for all probabilities
* NOTE: We assume that pErrorBounds has the size equal to the number of blocks in the partition
*/
err_state print_state_probs_partition(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const partition * P,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * vec,
                double error_bound,
                /*@observer@*/ /*@null@*/ const double * pErrorBounds)
{
	int j, i;
        int original_state_space_size;
        int lumped_state_space_size;
        int precision;

        if ( part_is_invalid(P) || NULL == vec ) {
                err_msg_4(err_PARAM, "print_state_probs_partition(%p,%p,%g,%p)",
                                (const void*) P, (const void*) vec, error_bound,
                                (const void *) pErrorBounds, err_ERROR);
        }

        errno = 0;
        original_state_space_size = part_unlumped_state_space_size_nt(P);
        lumped_state_space_size = part_lumped_state_space_size_nt(P);
	/* Get the printing pattern that agrees with the error bound precision */
        precision = get_error_bound_precision(-1, error_bound,
                        lumped_state_space_size, pErrorBounds);

		printf("( ");
		for( j = 0 ; j < original_state_space_size; j++ ){
                        if ( 0 != j ) {
				printf(", ");
			}

                        i = part_lump_state_id_nt(P, j);
                        printf("%1.*f", precision, vec[i]);
		}
		printf(" )");

        if ( 0 != errno ) {
                err_msg_4(err_FILE, "print_state_probs_partition(%p,%p,%g,%p)",
                                (const void*) P, (const void*) vec, error_bound,
                                (const void *) pErrorBounds, err_ERROR);
        }
        return err_OK;
}

/**
* This method prints the probability of one state only, with
* the unlumping done first.
* The error messages have to be printed as shown because the caller does not
* check the parameters.
* @param P the partition.
* @param index the state index in the original state space.
* @param vec the array of probabilities to be printed.
* @param error_bound the error bound for the given state probability
* @param pErrorBounds the error bounds for all probabilities
* NOTE: We assume that pErrorBounds has the size equal to the number of blocks in the partition
*/
err_state print_state_prob_partition(/*@observer@*/ const partition * P,
                state_index index_local, /*@observer@*/ const double * vec,
                double error_bound, /*@observer@*/ const double * pErrorBounds)
{
	int lumped_index;
        int original_state_space_size;
        int lumped_state_space_size;

        if ( part_is_invalid(P) || NULL == vec )
                err_msg_6(err_PARAM, "print_state_prob_partition(%p[%p],%d,%p,"
                                "%g,%p)", (const void *) P,
                                NULL == P ? NULL : (void *) P->blocks,
                                index_local, (const void *) vec, error_bound,
                                (const void *) pErrorBounds, err_ERROR);

        errno = 0;
        original_state_space_size = part_unlumped_state_space_size_nt(P);
        lumped_state_space_size = part_lumped_state_space_size_nt(P);

                if ( 0<=index_local && index_local<original_state_space_size ) {
			/* Get the state index in the lumped state space */
                        lumped_index = part_lump_state_id_nt(P, index_local);

			/* Get the printing pattern that agrees with the error bound precision */
                        printf("%1.*f", get_error_bound_precision(lumped_index,
                                        error_bound,
										lumped_state_space_size,
                                        pErrorBounds),

                                vec[lumped_index]);
		}else{
			printf("??\n");
			printf("WARNING: Invalid index %d, required to be in the [1, %d] interval.",
                                index_local + 1, original_state_space_size);
		}

        if ( 0 != errno ) {
                err_msg_5(err_FILE, "print_state_prob_partition(%p,%d,%p,%g,"
                                "%p)", (const void *) P, index_local,
                                (const void *) vec, error_bound,
                                (const void *) pErrorBounds, err_ERROR);
        }
        return err_OK;
}


/**
* The function creates a new partition with a single block. The block is NOT
* inserted into the list of potential splitters.
* Parameter: size = number of states to be members of the partition
* Result: Pointer to a new partition; NULL if some error has happened.
*/
static /*@null@*/ /*@only@*/ partition * part_new(state_count size)
                /*@ensures notnull result->blocks@*/
                /*@ensures isnull result->first_Sp, result->first_PredCl,
                        result->sum, result->pos@*/
                /*@ensures maxRead(result->blocks) == 1 /\
                        maxRead(result->ib) == size@*/
                /*@defines result->blocks, result->first_Sp,
                        result->first_PredCl, result->num_blocks, result->sum,
                        result->pos, result->ib[0]@*/
                /*@modifies internalState@*/
{
        pos_index i;
        /*@null@*/ /*@only@*/ partition * P = malloc(sizeof(partition)
                        - sizeof(struct id_block)
                        + sizeof(struct id_block) * size);
        /*@null@*/ /*@only@*/ block * B = malloc(sizeof(block));
        if ( NULL == P || NULL == B ) {
                err_msg_1(err_MEMORY, "part_new(%d)", size,
                                (free(B), free(P), NULL));
        }
        B->next = NULL;
        B->u.next_Sp = NULL;
        B->end = size;
        B->flags = 0;

        P->first_Sp = NULL;
        P->first_PredCl = NULL;
        P->num_blocks /* lumped size */ = 1;
        P->sum = NULL;
        P->pos = NULL;
        P->blocks = B;
        P->ib[0].id = 0;
        P->ib[0].b = B;
        for ( i = 1 ; i < size ; i++ ) {
                P->ib[i].id = i;
                P->ib[i].b = B;
        }
        return P;
}


/**
* The function sorts the states in a block into two parts: the not-bs-states and
* the bs-states. The block is NOT split.
* Parameters: P = partition to be split
*       B = block to be split
*       bs = bitset that indicates how to split
* Result: the array P->ib[].id is sorted: P->ib[begin(B) ... return value-1].id
* is not in the bitset; P->ib[return value ... end(B)-1].id is in the
* bitset. If some error happens, the function returns state_index_ERROR.
*/
static pos_index block_sort_bitset(partition * P, block * B,
                /*@observer@*/ const bitset * bs)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->pos@*/ /*@modifies P->ib[].id@*/
{
        pos_index left, right;
        state_index temp;

        left = part_block_begin_nt(B);
        right = B->end - 1;
        /* we now go through the states in the block quicksort-like: look for a
           wrong element in the left half; look for a wrong element in the right
           half; swap those two elements. Just one pass is needed, because we
           only split into two parts. */
        do {
                /* left <= right;
                   ib[block_begin(B) ... left-1].id are NOT in bs;
                   ib[right+1 ... block_end(B)-1].id are in bs. */
                while ( ! get_bit_val(bs, P->ib[left ].id) ) {
                        /* ib[block_begin(B) ... left].id are NOT in bs; */
                        ++left;
                        if ( left > right )
                                return left;
                }
                /* left <= right;
                   ib[begin(B) ... left-1].id are NOT in bs;
                   ib[left].id is in bs;
                   ib[right+1 ... end(B)-1].id are in bs. */
                do {
                        if ( left >= right )
                                return left;
                        if ( ! get_bit_val(bs, P->ib[right].id) )
                                break;
                        /* ib[right ... block_end(B)-1].id are in bs. */
                        --right;
                } while ( TRUE );
                /* left < right;
                   ib[block_begin(B) ... left-1].id are NOT in bs;
                   ib[left].id is in bs;
                   ib[right].id is NOT in bs;
                   ib[right+1 ... block_end(B)-1].id are in bs. */
                /* Now swap the states at ib[left].id and ib[right].id. */
                temp = P->ib[left].id;
                P->ib[left ].id = P->ib[right].id;
                P->ib[right].id = temp;
                ++left, --right;
        } while ( left <= right );
        return left;
}


/**
* The function adapts the block pointers for the states in block B and its
* predecessors.
*/
static /*@dependent@*/ /*@null@*/ block * part_split_and_adapt_block(
                partition * P, block * B, pos_index pos)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->pos@*/
                /*@modifies B->next, P->num_blocks, P->ib[].b@*/
{
        pos_index i;
        /*@dependent@*/ block * B_new;

        i = part_block_begin_nt(B);
        if ( (B_new = part_split_block(P, B, pos)) == NULL ) {
                err_msg_3(err_CALLBY, "part_split_and_adapt_block(%p,%p,%d)",
                                (void *) P, (void *) B, pos, NULL);
        }

        for ( ; i < pos ; i++ ) {
                P->ib[P->ib[i].id].b = B_new;
        }
        return B_new;
}

/**
* The function splits all blocks of the partition into two parts: the bs-states
* and the not-bs-states.
* Parameters: P = partition to be split
*       bs = bitset that serves to split the single sets
* Result: err_OK if everything went fine; err_ERROR otherwise.
*/
static err_state part_sort_bitset(partition *P, /*@observer@*/ const bitset *bs)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->pos@*/
                /*@modifies P->ib, P->num_blocks, *P->blocks@*/
{
        part_walk_blocks(P, b) {
                if ( 0 == (b->flags & ABSORBING) ) {
                        pos_index split;

                        if ( (split = block_sort_bitset(P, b, bs))
                                                        > part_block_end_nt(b)
                                || (split < part_block_end_nt(b)
                                    && (split < part_block_begin_nt(b)
                                        || (split > part_block_begin_nt(b)
                                            && part_split_and_adapt_block(P, b,
                                                        split) == NULL))) )
                        {
                                err_msg_2(err_CALLBY,"part_sort_bitset(%p,%p)",
                                                (void *) P, (const void *) bs,
                                                err_ERROR);
                        }
                }
        } end_part_walk_blocks;
        return err_OK;
}


/**
* The function refines each block in the partition according to the rewards.
* Parameters: P = partition to be refined
*       reward = pointer to array of doubles; reward[s] = reward rate of state s
* Result: err_OK if everything went fine; err_ERROR otherwise.
*/
err_state part_refine_rewards(/*@i1@*/ /*@null@*/ partition * P,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * reward)
{
        if ( part_is_invalid(P) || NULL == reward )
                err_msg_2(err_PARAM, "part_refine_rewards(%p,%p)", (void *) P,
                                (const void *) reward, err_ERROR);

        part_walk_blocks(P, b) {
                if ( 0 != (b->flags & PARTITIONED) ) {
                        err_msg_2(err_INCONSISTENT,"part_refine_rewards(%p,%p)",
                                        (void *) P, (const void *) reward,
                                        err_ERROR);
                }

                if ( 0 == (b->flags & ABSORBING) ) {
                        /*sort_and_split_block() requires that b be a splitter*/
                        if ( b->flags != SPLITTER ) {
                                b->flags = SPLITTER;
                                if ( NULL == b->u.next_Sp ) {
                                        b->u.next_Sp = P->first_Sp;
                                        P->first_Sp = b;
                                }
                        }
                        if ( err_state_iserror(sort_and_split_block(P, b,
                                                        reward, FALSE)) )
                        {
                                err_msg_2(err_CALLBY, "part_refine_rewards(%p,"
                                                "%p)", (void *) P,
                                                (const void*)reward, err_ERROR);
                        }
                }
        } end_part_walk_blocks;
        return err_OK;
}


/**
* The function creates a new partition that can serve as initial partition for
* formula-independent lumping: blocks are created depending on the labelling.
* Parameter: labellin = labelling that indicates how blocks have to be split.
* Result: Pointer to a new partition; NULL if some error has happened.
* All blocks in the new partition are splitters.
*/
/*@only@*/ /*@null@*/ partition * init_partition(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const labelling * labellin)
{
        partition * P;
        state_count size;
        int i;
        /*@observer@*/ /*@dependent@*/ const double * pStateRewards;

        if ( NULL == labellin ) {
                err_msg_1(err_PARAM, "init_partition(%p)",
                                (const void *) labellin, NULL);
        }

        size = labellin->ns;
        if ( (P = part_new(size)) == NULL ) {
                err_msg_1(err_CALLBY, "init_partition(%p)",
                                (const void *) labellin, NULL);
        }
        /* now split all blocks according to the bitsets */
        for ( i = labellin->n ; i-- > 0 ; ) {
                if ( err_state_iserror(part_sort_bitset(P, labellin->b[i])) ) {
                        err_msg_1(err_CALLBY, "init_partition(%p)",
                                        (const void *) labellin,
                                        (free_partition(P), NULL));
                }
        }

        /* set P->pos to the correct values */
        P->pos = (pos_index *) calloc((size_t)
                part_unlumped_state_space_size_nt(P), sizeof(pos_index));
        if ( NULL == P->pos ) {
                err_msg_1(err_MEMORY, "init_partition(%p)",
                                (const void *) labellin,
                                (free_partition(P), NULL));
        }
        for ( i = 1 ; i < size ; i++ ) {
                P->pos[P->ib[i].id] = i;
                /* calloc() already does the 0-th iteration through this loop */
        }

        pStateRewards = getStateRewards();
        if ( NULL != pStateRewards ) {
                /* split according to state rewards */
                if ( err_state_iserror(part_refine_rewards(P, pStateRewards)) ){
                        err_msg_1(err_CALLBY, "init_partition(%p)",
                                        (const void *) labellin,
                                        (free_partition(P), NULL));
                }
        }
        return P;
}


/**
* The function creates a new partition that will serve as initial partition for
* formula-dependent lumping: three or four blocks are created: a ``phi and not
* psi'' block, a ``nor phi nor psi'' block (absorbing), and finally:
* if interval == FALSE: a ``phi'' block (absorbing)
* if interval == TRUE:  a ``phi and psi'' block and a ``psi and not phi'' block
* Parameters: phi, psi, interval = bitsets and boolean as indicated
* Result: Pointer to a new partition; NULL if some error has happened.
* All blocks in the new partition are splitters.
*/
/*@only@*/ /*@null@*/ partition * init_partition_formula(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * phi,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * psi,
                BOOL interval)
{
        /*@only@*/ /*@null@*/ partition * P = NULL;
        pos_index split = 0;
        state_count size;
        /*@observer@*/ const double * pStateRewards;
        block * not_psi_block = NULL;
        state_index i;

        if ( NULL == phi || NULL==psi || bitset_size(phi) != bitset_size(psi) ){
                err_msg_3(err_PARAM, "init_partition_formula(%p,%p,%d)",
                                (const void *) phi, (const void *) psi,
                                (int) interval, NULL);
        }

        if ( (size = bitset_size(phi)) == state_count_ERROR
                                || (P = part_new(size)) == NULL
                        || (split = block_sort_bitset(P, P->blocks, psi)) < 0
                        || split > size )
        {
                err_msg_3(err_CALLBY, "init_partition_formula(%p,%p,%d)",
                                (const void *) phi, (const void *) psi,
                                (int) interval,
                                (NULL == P || (free_partition(P),FALSE), NULL));
        }
        if ( split < size ) {
                /* there are psi states */
                if ( 0 < split ) {
                        /* there are not-psi states */
                        if ( (not_psi_block = part_split_and_adapt_block(P,
                                                P->blocks, split)) == NULL )
                        {
                                err_msg_3(err_CALLBY, "init_partition_formula("
                                                "%p,%p,%d)", (const void *) phi,
                                                (const void*)psi, (int)interval,
                                                (free_partition(P), NULL));
                        }
                }
                if ( ! interval ) {
                        P->blocks->flags = ABSORBING;
                } else {
                        /* split the psi-states */
                        pos_index split_psi;

                        if ( (split_psi = block_sort_bitset(P, P->blocks, phi))
                                                < split || split_psi > size
                                        || (split_psi < size && split_psi>split
                                            && part_split_and_adapt_block(P,
                                                P->blocks, split_psi) == NULL) )
                        {
                                err_msg_3(err_CALLBY, "init_partition_formula("
                                                "%p,%p,%d)", (const void *) phi,
                                                (const void*)psi, (int)interval,
                                                (free_partition(P), NULL));
                        }
                }
        } else {
                not_psi_block = P->blocks;
        }
        if ( NULL != not_psi_block ) {
                /* there are not-psi states */
                pos_index split_not_psi;

                if ( (split_not_psi=block_sort_bitset(P,not_psi_block, phi)) < 0
                                                || split_not_psi > split )
                {
                        err_msg_3(err_CALLBY,"init_partition_formula(%p,%p,%d)",
                                        (const void *) phi, (const void *) psi,
                                        (int)interval,(free_partition(P),NULL));
                }
                if ( 0 < split_not_psi ) {
                        block * not_phi_nor_psi_block;

                        if ( split_not_psi < split ) {
                                if ( (not_phi_nor_psi_block =
                                                part_split_and_adapt_block(P,
                                                        not_psi_block,
                                                        split_not_psi)
                                           ) == NULL )
                                {
                                        err_msg_3(err_CALLBY, "init_partition_"
                                                "formula(%p,%p,%d)",
                                                (const void *) phi,
                                                (const void*)psi, (int)interval,
                                                (free_partition(P), NULL));
                                }
                        } else {
                                not_phi_nor_psi_block = not_psi_block;
                        }
                        not_phi_nor_psi_block->flags = ABSORBING;
                }
        }

        /* set P->pos to the correct values */
        P->pos = (pos_index *) calloc((size_t)
                part_unlumped_state_space_size_nt(P), sizeof(pos_index));
        if ( NULL == P->pos ) {
                err_msg_3(err_MEMORY, "init_partition_formula(%p,%p,%d)",
                                (const void *) phi, (const void *) psi,
                                (int) interval, (free_partition(P), NULL));
        }
        for ( i = 1 ; i < size ; i++ ) {
                P->pos[P->ib[i].id] = i;
                /* calloc() already does the 0-th iteration through this loop */
        }

        pStateRewards = getStateRewards();
        if ( NULL != pStateRewards ) {
                /* split according to state rewards */
                if ( err_state_iserror(part_refine_rewards(P, pStateRewards)) ){
                        err_msg_3(err_CALLBY, "init_partition_formula(%p,%p,"
                                        "%d)", (const void *) phi,
                                        (const void *) psi, (int) interval,
                                        (free_partition(P), NULL));
                }
        }

        return P;
}


/*****************************************************************************
name		: free_partition
role		: This method frees the partition.
@param		: partition *P: The partition.
******************************************************************************/
err_state free_partition(/*@only@*/ /*@i1@*/ /*@null@*/ partition * P) {
        /*@only@*/ block * B = NULL;

        if ( part_is_invalid(P) ) {
                err_msg_1(err_PARAM, "free_partition(%p)", (void *) P,
                                err_ERROR);
        }

	B = P->blocks;
        do {
                /*@only@*/ /*@null@*/ block *
		next = B->next;
                if ( err_state_iserror(free_block(B)) ) {
                        /*@-mustfreeonly@*/
                        err_msg_1(err_CALLBY, "free_partition(%p)", (void *) P,
                                (free(P->pos),free(P->sum),free(P), err_ERROR));
                        /*@=mustfreeonly@*/
                }
		B = next;
	}
        while ( NULL != B );

        free(P->pos);
        free(P->sum);
	free(P);
        return err_OK;
}


/**
* The function numbers the blocks in partition P. As many states as possible get
* a number n such that state n is also a member of the block.
* Parameter: P = partition
* Result: err_OK if everything went fine; err_ERROR otherwise.
*/
err_state part_number_blocks(/*@i1@*/ /*@null@*/ partition * P) {
        state_count i;

        if ( part_is_invalid(P) || NULL != P->first_Sp
                                        || NULL != P->first_PredCl )
        {
                err_msg_1(err_PARAM, "part_number_blocks(%p)", (void *) P,
                                err_ERROR);
        }
        part_walk_blocks(P, b) {
                b->u.row = state_index_NONE;
                b->flags &= ~ (SPLITTER | PARTITIONED);
        } end_part_walk_blocks;

#       if 1
        {
                /* The following implementation numbers the blocks in the order
                   of their member state with the smallest index. */
                state_index num = 0;
                for ( i = 0 ; i < part_unlumped_state_space_size_nt(P) ; i++ ) {
                        if ( state_index_NONE == part_lump_state_id_nt(P, i) ) {
                                block_of(P, i)->u.row = num;
                                if ( ++num>=part_lumped_state_space_size_nt(P) )
                                {
                                        break;
                                }
                        }
                }
        }
#       else
        {
                /* The following implementation tries to give each block the
                   number of a state in the block, if possible. */
                state_count free_number;
                state_index * free_id = calloc(
                                (size_t) part_lumped_state_space_size_nt(P),
                                sizeof(state_index));
                if ( NULL == free_id ) {
                        err_msg_1(err_MEMORY, "part_number_blocks(%p)",
                                        (void *) P, err_ERROR);
                }
                free_number = 0;
                for ( i = 0 ; i < part_lumped_state_space_size_nt(P) ; i++ ) {
                        if ( state_index_NONE == part_lump_state_id_nt(P, i) ) {
                                block_of(P, i)->u.row = i;
                        } else {
                                free_id[free_number++] = i;
                        }
                }
                if ( free_number > 0 ) {
                        const state_index * used_id = free_id;

                        for ( ; i < part_unlumped_state_space_size_nt(P) ; i++ )
                        {
                                if ( state_index_NONE ==
                                                part_lump_state_id_nt(P, i) )
                                {
                                        block_of(P, i)->u.row = *used_id++;
                                        if ( --free_number <= 0 )
                                                break;
                                }
                        }
                }
                free(free_id);
        }
#       endif

        return err_OK;
}


/**
* The function is a compare function for qsort(). It interprets its two
* arguments as pointers to integers and compares them directly.
*/
static int comp_state_index(const void * i1, const void * i2) {
        return *(const state_index *) i1 - *(const state_index *) i2;
}

/*****************************************************************************
name		: print_block_states
role		: This method prints state indices.
@param		: state *states: The list of states.
******************************************************************************/
static err_state print_block_states(/*@observer@*/ /*@temp@*/ const partition*P,
                /*@observer@*/ const block * B)
                /*@requires notnull P->blocks@*/
                /*@modifies fileSystem@*/
{
        /*@only@*/ /*@null@*/ state_index * ssort;
        /*@dependent@*/ state_index * ssort_temp;
        const struct id_block * id_temp;
        pos_index j;

        /* create a sorted list of the states in block B */
        ssort = (state_index *) calloc((size_t) part_block_size_nt(B),
                        sizeof(state_index));
        if ( NULL == ssort ) {
                err_msg_2(err_MEMORY, "print_block_states(%p,%p)",
                                (const void *) P, (const void *) B, err_ERROR);
        }
        ssort_temp = ssort;
        id_temp = &P->ib[part_block_begin_nt(B)];
        j = part_block_size_nt(B);
        printf("%d state%c: { ", j, 1 == j ? ' ' : 's');
        for ( ; j > 0 ; j-- ) {
                *ssort_temp++ = id_temp++->id;
        }
        qsort(ssort, (size_t) part_block_size_nt(B), sizeof(state_index),
                        comp_state_index);
        for ( j = 0 ; j < part_block_size_nt(B) ; j++ ) {
                if ( 0 != j ) {
			printf(", ");
		}
                printf("%d", ssort[j] + 1);
                if ( block_of(P, ssort[j]) != B ) {
                        /* block_of is inconsistent */
                        (void) putchar('B');
                }
                if ( NULL != P->pos && P->ib[P->pos[ssort[j]]].id != ssort[j] ){
                        /* pos is inconsistent */
                        (void) putchar('P');
                }
	}
        free(ssort);
	printf(" }\n");
        return err_OK;
}


/*****************************************************************************
name		: print_partition
role		: This method prints the blocks and states of the partition.
@param		: partition *P: The partition.
******************************************************************************/
err_state print_partition(/*@observer@*/ /*@i1@*/ /*@null@*/ const partition *P)
{
        /*@only@*/ /*@null@*/ const block ** bsort = NULL;
                                        /* pointer to array of block pointers */
        pos_index i;

        if ( part_is_invalid(P) ) {
                err_msg_1(err_PARAM, "print_partition(%p)", (const void *) P,
                                err_ERROR);
        }

        errno = 0;
	printf("Partition:\n");

        /* create a sorted list of the blocks, so that they will be printed in
           order */
        bsort = (const block **) calloc((size_t)
                part_lumped_state_space_size_nt(P), sizeof(const block *));
        if ( NULL == bsort ) {
                err_msg_1(err_MEMORY, "print_partition(%p)", (const void *) P,
                                err_ERROR);
        }
        part_walk_blocks(P, B) {
                if ( 0 > part_lump_state_block_nt(P, B)
                                || part_lumped_state_space_size_nt(P)
                                                <= part_lump_state_block_nt(P,B)
                                || NULL!=bsort[part_lump_state_block_nt(P, B)] )
                {
                        err_msg_1(err_INCONSISTENT, "print_partition(%p)",
                                        (const void *) P,
                                        (free(bsort), err_ERROR));
                }
                bsort[part_lump_state_block_nt(P, B)] = B;
        } end_part_walk_blocks;

        /* now print the blocks in order */
        for ( i = 0 ; i < part_lumped_state_space_size_nt(P) ; i++ ) {
                /*@dependent@*/ const block * B = bsort[i];

                printf("state %d: ", part_lump_state_block_nt(P, B) + 1);
                if ( err_state_iserror(print_block_states(P, B)) ) {
                        err_msg_1(err_CALLBY, "print_partition(%p)",
                                        (const void *) P,
                                        (free(bsort), err_ERROR));
                }
	}
        free(bsort);
        if ( 0 != errno ) {
                err_msg_1(err_FILE, "print_partition(%p)", (const void *) P,
                                err_ERROR);
        }
        return err_OK;
}


/**
* The function prints an extensive view of the partition (with more information
* about the structure of the partition than the above function).
*/
/* The function is, however, never used, and therefore, I commented it out.
#define LINELEN 90

static char space20[] = "                    ";

err_state part_print_extended(const partition * P) {
        char line1[LINELEN+10], line2[LINELEN+10];
        int len1, format_len;
        state_index j, i;
        const block ** bsort = NULL;

        if ( part_is_invalid(P) ) {
                err_msg_1(err_PARAM, "print_partiton_extended(%p)",
                                (const void *) P, err_ERROR);
        }

        bsort = (const block **) calloc((size_t)
                part_lumped_state_space_size_nt(P), sizeof(const block *));
        if ( NULL == bsort ) {
                err_msg_1(err_MEMORY, "print_partiton_extended(%p)",
                                (const void *) P, err_ERROR);
        }

        j = part_lumped_state_space_size_nt(P);
        part_walk_blocks(P, b) {
                bsort[--j] = b;
        } end_part_walk_blocks;

        strcpy(line1, "{ ");
        strcpy(line2, "  ");
        len1 = 2;
        for ( j = 0 ; j < part_lumped_state_space_size_nt(P) ; j++ ) {
                const block * b = bsort[j];
                if ( NULL == b ) {
                        err_msg_1(err_INCONSISTENT, "print_partiton_extended("
                                        "%p)", (const void *) P,
                                        (free(bsort), err_ERROR));
                }

                if ( 0 != j ) {
                        strcpy(&line1[len1], ",  ");
                        strncat(&line2[len1], space20, 3-strlen(&line2[len1]));
                        if ( (len1 += 3) >= LINELEN ) {
                                printf("%s\n%s\n\n", line1, line2);
                                line1[0] = line2[0] = '\0';
                                len1 = 0;
                        }
                }
                strcpy(&line1[len1], "{");
                if ( NULL == P->first_Sp && 0 <= part_lump_state_block_nt(P, b)
                                && part_lump_state_block_nt(P, b)
                                        < part_lumped_state_space_size_nt(P) )
                {
                        sprintf(&line2[len1], "%d", b->u.row + 1);
                } else {
                        strcpy(&line2[len1], &space20[19]);
                }
                ++len1;
                if ( 0 != (b->flags & SPLITTER) ) {
                        strcat(&line2[len1], "S");
                }
                if ( 0 != (b->flags & ABSORBING) ) {
                        strcat(&line2[len1], "A");
                }
                for ( i = part_block_begin_nt(b); i < b->end ; ++i ) {
                        if ( i > part_block_begin_nt(b) ) {
                                strcpy(&line1[len1], ",");
                                if ( '\0' == line2[len1] ) {
                                        strcpy(&line2[len1], &space20[19]);
                                }
                                if ( ++len1 >= LINELEN ) {
                                        printf("%s\n%s\n\n", line1, line2);
                                        line1[0] = line2[0] = '\0';
                                        len1 = 0;
                                }
                        }
                        format_len = sprintf(&line1[len1],"%d",P->ib[i].id+1);
                        if ( format_len > (int) sizeof(line1) - 1 - len1 ) {
                                format_len = (int) sizeof(line1) - 1 - len1;
                        }
                        strncat(&line2[len1], space20,
                                        format_len - strlen(&line2[len1]));
                        len1 += format_len;
                        if ( NULL != P->pos && P->pos[P->ib[i].id] != i ) {
                                line2[len1 - 1] = '*';
                        }
                }
                strcpy(&line1[len1], "}");
                if ( '\0' == line2[len1] ) {
                        strcpy(&line2[len1], &space20[19]);
                }
                ++len1;
        }
        printf("%s }\n%s\n", line1, line2);
        free(bsort);
        return err_OK;
}


Another format to print the partition in detail:
{
        const block * B;
        state_index i;

        if ( part_is_invalid(P) ) {
                err_msg_1(err_PARAM, "print_partiton_extended(%p)",
                                (const void *) P, err_ERROR);
        }

        errno = 0;
        printf("Partition: blocks=%p, first_Sp=%p, first_PredCl=%p, "
                        "num_blocks=%d, sum=%p, pos=%p\n", (void*) P->blocks,
                        (void *) P->first_Sp, (void *) P->first_PredCl,
                        part_lumped_state_space_size_nt(P), (void *) P->sum,
                        (void *) P->pos);
        printf("States in the partition:\n");
        for ( i = 0 ; i < part_unlumped_state_space_size_nt(P) ; i++ ) {
                pos_index pos_i = NULL == P->pos ? -1 : P->pos[i];
                printf("State %d: pos=%d, block=%p, sum=%g, ib[pos[%d]].id==%d\n",
                                i, pos_i, (const void *) block_of(P, i),
                                P->sum[i], i, P->ib[pos_i].id);
        }
        printf("Blocks in the partition:\n");
        i = part_lumped_state_space_size_nt(P);
        part_walk_blocks(P, B) {
                printf("Block %d, address=%p: next=%p, next_Sp=%p/row=%d, "
                                "end=%d, splitter=%d, end_Pred=%d\n", --i,
                                (const void *) B, (const void *) B->next,
                                (void *) B->u.next_Sp,
                                part_lump_state_block_nt(P, B), B->end,
                                0!=(B->flags&SPLITTER), B->end_Pred_xor_abs);
        } end_part_walk_blocks;
        if ( 0 != errno ) {
                err_msg_1(err_FILE, "print_partition_extended(%p)",
                                (const void *) P, err_ERROR);
        }
        return err_OK;
} */
