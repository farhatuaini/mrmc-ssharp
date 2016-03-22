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
*	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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

#ifndef PARTITION_H
#define PARTITION_H

#include "label.h"

/* flags for block */
#define SPLITTER	0x01
#define PARTITIONED	0x02
#define ABSORBING	0x04

typedef state_index pos_index;

/**
*	STRUCTURE
* A block in a partition.
* next_Sp               : pointer to next block in list of potential splitters;
*                         if this pointer is == NULL, the block is the last in
*                         the list of splitters, or it is not a splitter.
* next_PredCl           : pointer to next block in list PredCl (i.e. predecessor
*                         classes of the current splitter); if this pointer is
*                         == NULL, the block is the last in the list of
*                         predecessor classes, or it is not a predecessor.
* @member row		: row index in new (lumped) matrix.
* end                   : indicates the end of the block. The end is the first
*                         position in P->id that does not belong to this block;
*                         the beginning of the block is the end of this->next.
*                         So, the block contains the states
*                         P->ib[next->end ... end-1].id.
* @member next		: pointer to next block.
*/
typedef struct block_elem block;
struct block_elem {
        union next_Sp_PredCl_and_row {
                /* The following three fields are never used at the same time,
                   so we save space by declaring them as a union. */
                /*@null@*/ /*@dependent@*/ block * next_Sp;
                /*@null@*/ /*@dependent@*/ block * next_PredCl;
                state_index row;
        } u;
	int flags;
        pos_index end;
        /*@null@*/ /*@owned@*/
	block *next;
};


/**
*	STRUCTURE
* @member blocks		: linked-list of blocks in the partition.
* first_Sp              : pointer to first class in the list of (potential)
*                         splitters; see the field next_Sp in struct block_elem.
* first_PredCl          : pointer to first class in the list of predecessors of
*                         the current splitter; see the field next_PredCl in
*                         struct block_elem.
* @member num_blocks	: number of blocks in the partition.
* sum                   : pointer to array of doubles that contains, for each
*                         state, the probability to enter the splitter.
* ib[pos].id            : array that contains state indexes in a special order,
*                         as explained with struct block_elem.
* ib[state].b           : array that contains, for each state, the block to
*                         which it belongs.
* Note that the order of the two subarrays of ib is different; they are put
* together anyway to make optimal use of the struct hack.
*/
typedef struct{
        /*@owned@*/
	block *blocks;
        /*@dependent@*/ /*@null@*/ block * first_Sp;
        /*@dependent@*/ /*@null@*/ block * first_PredCl;
        state_count num_blocks;
        /*@only@*/ /*@null@*/ double * sum;
        /*@only@*/ /*@null@*/ pos_index * pos;

        /* struct hack: saves one pointer dereference. In practice, the
           following array contains an element for each state, not only one
           element. */
        struct id_block {
                state_index id;
                /*@dependent@*/ block * b;
        } ib[1];
} partition;

/**
* This macro returns the position index of the last state in the block PLUS ONE.
* Parameter: B = block
* Result: position index.
*/
extern pos_index part_block_end_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const block * B)
                /*@modifies nothing@*/;
#define part_block_end_nt(B) ((void) (B)->u.next_Sp, (void) (B)->u.next_PredCl,\
                (void) ((B)->u.row + (B)->flags), (void) (B)->next, \
                (const pos_index) (B)->end)
/* The superfluous references to fields of B act as an additional, hidden type
   check. They should be optimized away by the compiler. */

/**
* This macro returns the position index of the first state in the block.
* The suffix "_nt" indicates that the macro does not test its parameters for
* validity.
* Parameter: B = block
* Result: position index.
*/
extern pos_index part_block_begin_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const block * B)
                /*@modifies nothing@*/;
#define part_block_begin_nt(B) ((const pos_index) \
                (NULL == (B)->next ? 0 : part_block_end_nt((B)->next)))

/**
* This macro returns the size of the block.
* Parameter: B = block
* Result: size; state_count_ERROR if some error has happened.
*/
extern state_count part_block_size_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const block * B)
                /*@modifies nothing@*/;
#define part_block_size_nt(B) ((const state_count) \
                ((B)->end - part_block_begin_nt((B))))

/**
* This macro makes a few simple tests to check whether the partition is (in)-
* valid; it can be used to check a parameter.
*/
extern /*@nullwhentrue@*/ BOOL part_is_invalid(
                /*@observer@*//*@temp@*//*@sef@*//*@null@*/ const partition *P);
#define part_is_invalid(P) (sizeof(partition) != sizeof(*(P)) || NULL == (P) \
                || NULL == (P)->blocks || ((void) (P)->blocks->u.next_Sp, \
                (void) (P)->first_Sp, (void) (P)->first_PredCl, \
                (void) ((P)->num_blocks + (P)->ib[0].id), (void) (P)->sum, \
                (void) (P)->pos, (void) (P)->ib[0].b->flags, \
                FALSE))

/**
* This macro finds some state in the block B.
* The suffix "_nt" indicates that the macro does not test its parameters for
* validity.
*/
extern state_index part_unlump_state_block_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition * P,
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const block * B)
                /*@requires notnull P->blocks@*/
                /*@modifies nothing@*/;
#define part_unlump_state_block_nt(P,B) ((void) (P)->blocks->u.next_Sp, \
                (void) (P)->first_Sp, (void) (P)->first_PredCl, \
                (void) ((P)->num_blocks + (P)->ib[0].id), (void) (P)->sum, \
                (void) (P)->pos, (void) (P)->ib[0].b->flags, \
                (void) (B)->u.next_Sp, (void) (B)->u.next_PredCl, \
                (void) ((B)->u.row + (B)->flags), (void) (B)->next, \
                (const state_index) (P)->ib[(B)->end - 1].id)

extern state_index part_unlump_state_block(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition * P,
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const block * B)
                /*@requires notnull P->blocks@*/
                /*@modifies nothing@*/;
#define part_unlump_state_block(P,B) (part_is_invalid((P)) \
                        || sizeof(block) != sizeof(*(B)) || NULL == (B) \
                ? err_macro_2(err_PARAM, "part_unlump_state_block(%p,%p)", \
                  (const void *) (P), (const void *) (B), state_index_ERROR) \
                : part_unlump_state_block_nt((P), (B)))

/**
* This macro returns the block of a state.
* @param	: partition *P: The partition.
* @param	: int i: The state index.
*/
extern /*@dependent@*/ block * block_of(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition * P,
                state_index state) /*@requires notnull P->blocks@*/
                /*@modifies nothing@*/;
#define block_of(P,state) ((void) (P)->blocks->u.next_Sp, \
                (void) (P)->first_Sp, (void) (P)->first_PredCl, \
                (void) ((P)->num_blocks + (P)->ib[0].id), (void) (P)->sum, \
                (void) (P)->pos, (void) (P)->ib[0].b->flags, \
                (block * const) (P)->ib[(state)].b)

/**
* This macro finds the number of block B. (It is not used internally because it
* tightly checks its parameters.)
*/
extern state_index part_lump_state_block_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition * P,
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const block * B)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies nothing@*/;
#define part_lump_state_block_nt(P,B) ((void) (P)->blocks->u.next_Sp, \
                (void) (P)->first_Sp, (void) (P)->first_PredCl, \
                (void) ((P)->num_blocks + (P)->ib[0].id), (void) (P)->sum, \
                (void) (P)->pos, (void) (P)->ib[0].b->flags, \
                (void) (B)->u.next_Sp, (void) (B)->u.next_PredCl, \
                (void) ((B)->u.row + (B)->flags), (void) (B)->next, \
                (const state_index) (B)->u.row)

extern state_index part_lump_state_block(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition * P,
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const block * B)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies nothing@*/;
#define part_lump_state_block(P,B) (part_is_invalid((P)) \
                        || sizeof(block) != sizeof(*(B)) || NULL == (B) \
                ? err_macro_2(err_PARAM, "part_lump_state_block(%p,%p)", \
                  (const void *) (P), (const void *) (B), state_index_ERROR) \
                : part_lump_state_block_nt((P), (B)))

/**
* This method prints the position of the bits in the original
* bitset whose value is 1.
* @param P: The partition.
* @param b: The bitset corresponding to the lumped state space.
*/
extern err_state print_original_states(/*@observer@*/ const partition * P,
                /*@observer@*/ const bitset * b)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies fileSystem@*/;

/**
* This method prints TRUE if index is in the unlumped bitset 'b',
* otherwise false.
* @param P: The partition.
* @param b: The bitset corresponding to the lumped state space.
* @param index_local: The original-state-space state index
*/
extern err_state print_original_state(/*@observer@*/ const partition * P,
                /*@observer@*/ const bitset * b, state_index index_local)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies fileSystem@*/;

/**
* Print the vector of probability errors.
* @param P the partition.
* @param pErrorValues the array error values.
* NOTE: We assume that the size of pValues is equal to the number of partition
* blocks.
*/
extern err_state print_error_probs_partition(/*@observer@*/ const partition * P,
                /*@observer@*/ const double * pErrorValues)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies fileSystem@*/;

/**
* Print the unlumped vector of probabilities.
* @param P the partition.
* @param vec the array of probabilities to be printed.
* @param error_bound the error bound for all probabilities
* @param pErrorBounds the error bounds for all probabilities
* NOTE: We assume that pErrorBounds has the size equal to the number of blocks in the partition
*/
extern err_state print_state_probs_partition(/*@observer@*/ const partition * P,
                /*@observer@*/ const double * vec, double error_bound,
                /*@observer@*/ /*@null@*/ const double * pErrorBounds)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies fileSystem@*/;

/**
* This method prints the probability of one state only, with
* the unlumping done first.
* @param P the partition.
* @param index_local the state index in the original state space.
* @param vec the array of probabilities to be printed.
* @param error_bound the error bound for the given state probability
* @param pErrorBounds the error bounds for all probabilities
* NOTE: We assume that pErrorBounds has the size equal to the number of blocks in the partition
*/
extern err_state print_state_prob_partition(/*@observer@*/ const partition * P,
                state_index index_local, /*@observer@*/ const double *vec,
                double error_bound,
                /*@observer@*/ /*@null@*/ const double * pErrorBounds)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies fileSystem@*/;

/**
* The function splits a block into two parts. It creates a new block and moves
* the states in ib[block_begin(B) ... pos-1].id to the new block.
* Parameters: P = partition
*               B = block
*               pos = position where the block has to be split.
*               set_b = if TRUE, block_split also sets P->ib[...].b to the
*                       new block.
* Result: err_OK if everything went fine; err_ERROR if some error has happened.
*/
extern /*@dependent@*/ /*@null@*/ block * part_split_block(partition * P,
                block * B, pos_index pos)
                /*@requires notnull P->blocks@*/
                /*@modifies P->num_blocks, B->next@*/;

/**
* The function refines each block in the partition according to the rewards.
* Parameters: P = partition to be refined
*       reward = pointer to array of doubles; reward[s] = reward rate of state s
* Result: err_OK if everything went fine; err_ERROR otherwise.
*/
extern err_state part_refine_rewards(partition * P,
                /*@observer@*/ const double * reward)
                /*@requires notnull P->blocks, P->pos@*/
                /*@modifies P->first_Sp, *P->blocks, P->ib, P->num_blocks,
                        *P->pos @*/;

/**
* The function creates a new partition that can serve as initial partition for
* formula-independent lumping: blocks are created depending on the labelling.
* Parameter: labellin = labelling that indicates how blocks have to be split.
* Result: Pointer to a new partition; NULL if some error has happened.
*/
extern /*@only@*/ /*@null@*/ partition * init_partition(
                /*@observer@*/ const labelling * labellin)
                /*@ensures notnull result->blocks, result->pos@*/
                /*@ensures isnull result->first_PredCl, result->sum@*/
                /*@modifies nothing@*/;

extern /*@only@*/ /*@null@*/ partition * init_partition_formula(
                /*@observer@*/ const bitset * phi,
                /*@observer@*/ const bitset * psi, BOOL interval)
                /*@ensures notnull result->blocks, result->pos@*/
                /*@ensures isnull result->first_PredCl, result->sum@*/
                /*@modifies nothing@*/;

/**
* This method frees the partition.
* @param	: partition *P: The partition.
*/
extern err_state free_partition(/*@only@*/ partition * P)
                /*@requires notnull P->blocks@*/ /*@modifies *P@*/;

/**
* This method frees the block.
* @param	: block *B: The block.
*/
extern err_state free_block(/*@only@*/ /*@sef@*/ block *B) /*@modifies B@*/;
#define free_block(B) (sizeof(block) != sizeof(*(B)) || NULL == (B) \
                        || ((void) (B)->u.next_Sp, (void) (B)->u.next_PredCl, \
                            (void) ((B)->u.row + (B)->flags), (void) (B)->next,\
                            FALSE) \
                ? err_macro_1(err_PARAM,"free_block(%p)",(void*)(B),err_ERROR) \
                : (free((B)), err_OK))

/**
* The function numbers the blocks in partition P. As many states as possible get
* a number n such that state n is also a member of the block.
* Parameter: P = partition
* Result: err_OK if everything went fine; err_ERROR otherwise.
*/
extern err_state part_number_blocks(partition * P)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@ensures notnull P->blocks@*/
                /*@modifies P->ib[].b->u, *P->blocks->next@*/;

/**
* This method prints the blocks and states of the partition.
* @param	: partition *P: The partition.
*/
extern err_state print_partition(/*@observer@*/ /*@temp@*/ const partition * P)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@ensures notnull P->blocks@*/
                /*@modifies fileSystem@*/;

/**
* The macro gets the size of the lumped state space.
*/
extern state_count part_lumped_state_space_size_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition *P)
                /*@requires notnull P->blocks@*/
                /*@modifies nothing@*/;
#define part_lumped_state_space_size_nt(P) ((void) (P)->blocks->u.next_Sp, \
                (void) (P)->first_Sp, (void) (P)->first_PredCl, \
                (void) ((P)->num_blocks + (P)->ib[0].id), (void) (P)->sum, \
                (void) (P)->pos, (void) (P)->ib[0].b->flags, \
                (const state_count) (P)->num_blocks)

/**
* Here we get the original state-space size,
* that has changed because of lumping.
* The suffix "_nt" indicates that the macro does not test its parameters for
* validity.
* @return the original state-space size
*/
extern state_count part_unlumped_state_space_size_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition * P)
                /*@requires notnull P->blocks@*/
                /*@modifies nothing@*/;
#define part_unlumped_state_space_size_nt(P) ((void) (P)->blocks->u.next_Sp, \
                (void) (P)->first_Sp, (void) (P)->first_PredCl, \
                (void) ((P)->num_blocks + (P)->ib[0].id), (void) (P)->sum, \
                (void) (P)->pos, (void) (P)->ib[0].b->flags, \
                (const state_count) (P)->blocks->end)

extern state_count get_unlumped_state_space_size(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition *pPartition)
                /*@requires notnull pPartition->blocks@*/
                /*@modifies nothing@*/;
#define get_unlumped_state_space_size(pPartition) (part_is_invalid(pPartition) \
                ? err_macro_1(err_PARAM, "get_unlumped_state_space_size(%p)", \
                                (const void *) (pPartition), state_count_ERROR)\
                : part_unlumped_state_space_size_nt((pPartition)))

/**
* This method is used to retrieve the state index in the lumped state space.
* The suffix "_nt" indicates that the macro does not test its parameters for
* validity.
* @param pPartition the state-space partition
* @param state_idx the state index in the unlumped state space
* @return "-1" if "state_idx" is not within the state-index range of
*	the original state space, and the state index in the lumped space otherwise.
*/
extern state_index part_lump_state_id_nt(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition * P,
                state_index state_idx)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl@*/
                /*@modifies nothing@*/;
#define part_lump_state_id_nt(P,state_idx) ((void) (P)->blocks->u.next_Sp, \
                (void) (P)->first_Sp, (void) (P)->first_PredCl, \
                (void) ((P)->num_blocks + (P)->ib[0].id), (void) (P)->sum, \
                (void) (P)->pos, (void) (P)->ib[0].b->flags, \
                (const state_index) (P)->ib[(state_idx)].b->u.row)

extern state_index getLumpedStateIndex(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const partition *pPartition,
                /*@sef@*/ state_index state_idx)
                /*@requires notnull pPartition->blocks@*/
                /*@requires isnull pPartition->first_Sp,
                        pPartition->first_PredCl@*/
                /*@modifies nothing@*/;
#define getLumpedStateIndex(pPartition,state_idx) ( \
                part_is_invalid((pPartition)) || 0 > (state_idx) \
                || (pPartition)->blocks->end <= (state_idx) \
                ? err_macro_2(err_PARAM, "getLumpedStateIndex(%p,%d)", \
                                (const void *) (pPartition), (state_idx), \
                                (state_index) -1) \
                : part_lump_state_id_nt((pPartition), (state_idx)))

/**
* The following iterator macros walk through all blocks of a given partition.
*/

/*@iter part_walk_blocks(sef partition * P, yield dependent block * m_B)@*/
#define part_walk_blocks(P,m_B) \
                { \
                        block * m_B, * m_next__blocks; \
                        if ( part_is_invalid((P)) ) \
                                exit(err_macro_1(err_PARAM, \
                                        "part_walk_block(%p,B)", \
                                        (const void *) (P), EXIT_FAILURE)); \
                        m_next__blocks = (P)->blocks; \
                        for ( m_B = m_next__blocks ; NULL != (m_B) ; \
                                                m_B = m_next__blocks ) { \
                                m_next__blocks = (m_B)->next;

#define end_part_walk_blocks \
                        } \
                        /*@-noeffect@*/ (void) m_next__blocks; /*@=noeffect@*/ \
                }

#endif
