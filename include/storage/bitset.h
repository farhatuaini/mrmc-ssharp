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
*	Authors: Maneesh Khattri, Ivan Zapreev, David N. Jansen
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
*	Source description: Manage a Bitset - For the satisfaction relation
*	in MCC.
*	Uses: Definition of bitset - bitset.h
*/

#ifndef BITSET_H
#define BITSET_H

#include "error.h"

#include <stdlib.h>

/* This is the type of one block in the bitset */
typedef unsigned int BITSET_BLOCK_TYPE;
/* This is the size of the bitset block */
#define BITSET_BLOCK_SIZE (sizeof(BITSET_BLOCK_TYPE) * 8)

/**
*This are the constant bits
*/
#define BIT_ON	((BITSET_BLOCK_TYPE) 1)
#define BIT_OFF	((BITSET_BLOCK_TYPE) 0)

/**
* The bitset structure
* @member n			: The number of bits in this bitset.
* @member bytesp		: An array of BITSET_BLOCK_TYPE.
*/
typedef /*@abstract@*/ /*@immutable@*/ struct bitset
{
	state_count n;
	BITSET_BLOCK_TYPE bytesp[1];
} bitset;

typedef /*@only@*/ bitset * bitset_ptr;

/**
* Get a new bitset.
* A bitset is composed of an array of BITSET_BLOCK_SIZE-bit unsigned integers.
* @param n indicates the size of the required bitset.
* @return returns a pointer to a new bitset. NULL if some error has occurred.
*/
extern /*@null@*/ /*@only@*/ bitset * get_new_bitset(state_count n)
                /*@modifies nothing@*/;

/**
* This function extends a bitset, similar to realloc().
* We assume that the bitset size may only be increased.
* @param pBitset the bitset to increase
* @param new_size the new bitset size, should be: "new_size >= pBitset->n"
* @param isFillNewBlocksWithOne if TRUE then the new bits will be set.
* @result the new address of the bitset. If some error occurred, the function
*		returns NULL and does not change anything. In particular, the
*		old bitset is still allocated.
* (former function name: extend_bitset)
*/
extern /*@null@*/ /*@only@*/ bitset * bs_extend(
		/*@null@*/ /*@only@*/ bitset * pBitset,
		state_count new_size, BOOL isFillNewBlocksWithOne)
                /*@modifies *pBitset@*/;

/**
* This function frees a bitset.
* @param a: the bitset that has to be freed.
*/

/* prototype for macro: helps splint to check for correct usage */
extern err_state free_bitset(/*@only@*/ /*@out@*/ /*@sef@*/ bitset * a)
                /*@modifies a@*/;
#define free_bitset(a) (sizeof(bitset) != sizeof(*(a)) || NULL == (a) \
                ? err_macro_2(err_PARAM, "free_bitset(%p[%d])", (void *) (a), \
                                (NULL) != (a) ? /*@-usedef@*/ \
                                                ((void) (a)->bytesp[0], (a)->n)\
                                                /*@=usedef@*/ \
                                              : 0, err_ERROR)\
                : (free((a)), err_OK))
/* The superfluous reference to (a)->bytesp[0] serves as an additional, hidden
   type check. It will be optimized away by the compiler. David N. Jansen. */

/**
* This method just copies elements from one bitset into another.
* @param from_bitset the bitset to be copied .
* @param to_bitset the bitset to copy to.
* NOTE: The two bitsets must have the same size.
*/
extern err_state copy_bitset(/*@observer@*/ const bitset * from_bitset,
                bitset * to_bitset)
		/*@modifies *to_bitset@*/;

/**
* Allocate a bitset and set it to the disjunction of two bitsets.
* @param a one of the operands to the disjunction operation.
* @param b the second operand to the disjuction operation.
* @return returns a pointer to a new bitset that contains the result of
*		disjunction. NULL if some error occurred, e. g. a and b have
*		different sizes.
*/
extern /*@null@*/ /*@only@*/ bitset * or(/*@observer@*/ const bitset * a,
                /*@observer@*/ const bitset * b) /*@modifies nothing@*/;

/**
* Obtain the disjunction of two bitsets and put result in arg2.
* @param a one of the operands to the disjunction operation.
* @param b the second operand to the disjunction operation. Upon return, it will
*		contain the disjunction.
*/
extern err_state or_result(/*@observer@*/ const bitset * a, bitset * reslt)
		/*@modifies *reslt@*/;

/**
* Allocate a bitset and set it to the exclusive-or of two bitsets.
* @param a one of the operands to the exclusive-or operation.
* @param b the second operand to the exclusive-or operation.
* @return returns a pointer to a new bitset that contains the result of
*		exclusive-or. NULL if some error occurred, e. g. a and b have
*		different sizes.
*/
extern /*@null@*/ /*@only@*/ bitset * xor(/*@observer@*/ const bitset * a,
                /*@observer@*/ const bitset * b) /*@modifies nothing@*/;

/**
* Obtain the exclusive-or of two bitsets and put result in arg2.
* @param a one of the operands to the exclusive-or operation.
* @param b the second operand to the exclusive-or operation. Upon return, it
*		will contain the exclusive-or.
*/
extern err_state xor_result(/*@observer@*/ const bitset * a, bitset * reslt)
		/*@modifies *reslt@*/;

/**
* Allocate a bitset and set it to the conjunction of two bitsets.
* @param a one of the operands to the conjunction operation.
* @param b the second operand to the conjuction operation.
* @return returns a pointer to a new bitset that contains the result of
*		conjunction. NULL if some error occurred, e. g. a and b have
*		different sizes.
*/
extern /*@null@*/ /*@only@*/ bitset * and(/*@observer@*/ const bitset * a,
                /*@observer@*/ const bitset * b) /*@modifies nothing@*/;

/**
* Obtain the conjunction of two bitsets and put result in arg2.
* @param a one of the operands to the conjunction operation.
* @param b the second operand to the conjunction operation. Upon return, it will
*		contain the conjunction.
*/
extern err_state and_result(/*@observer@*/ const bitset * a, bitset * reslt)
		/*@modifies *reslt@*/;

/**
* Allocate a bitset and set it to the complement of a bitset.
* @param a the operand to the complement operation.
* @return returns a pointer to a new bitset that contains the result of
*		complementation. NULL if some error occurred.
*/
extern /*@null@*/ /*@only@*/ bitset * not(/*@observer@*/ const bitset * a)
                /*@modifies nothing@*/;

/**
* Obtain the complement of a bitsets and put result in the argument.
* @param reslt the operand to the complement operation. Upon return, it will
*		contain the complement.
*/
extern err_state not_result(bitset * reslt)
		/*@modifies *reslt@*/;

/**
* Gets the value of a specified bit of a bitset.
* @param a One of the bits of this bitset is of interest.
* @param i The position of the bit of interest in the bitset.
* @return TRUE if the bit is set to 1; FALSE otherwise; -1 if some error
*		occurred.
*/

/* prototype for macro: helps splint to check for correct usage */
extern BOOL get_bit_val(/*@sef@*/ /*@observer@*/ const bitset * a,
                /*@sef@*/ state_index i)
		/*@modifies nothing@*/;
#define get_bit_val(a,i) \
                (sizeof(bitset) != sizeof(*(a)) || NULL == (a) || \
                                        (unsigned) (i) >= (unsigned) (a)->n \
		 ? err_macro_3(err_PARAM, "get_bit_val(%p[%d],%d)", \
				(const void *) (a), NULL != (a) ? (a)->n : 0, \
				(i), -1) \
		 : BIT_OFF != ((a)->bytesp[(i) / BITSET_BLOCK_SIZE] \
			       & ((BITSET_BLOCK_TYPE) 1 \
				  << ((unsigned) (i) % BITSET_BLOCK_SIZE))))

/**
* Sets the value of a specified bit of a bitset.
* @param a One of the bits of this bitset is to be changed.
* @param i The position of the bit of interest in the bitset.
* @param val The new value for the bit of interest: 0 = clear bit; other value =
*		set bit.
*/

/* prototype for macro: helps splint to check for correct usage */
extern err_state set_bit_val(/*@sef@*/ bitset * a, /*@sef@*/ state_index i,
		/*@sef@*/ BITSET_BLOCK_TYPE val) /*@modifies *a@*/;
#define set_bit_val(a,i,val) \
                (sizeof(bitset) != sizeof(*(a)) || NULL == (a) || \
                                        (unsigned) (i) >= (unsigned) (a)->n \
		 ? err_macro_4(err_PARAM, "set_bit_val(%p[%d],%d,%u)", \
				(void *) (a), NULL != (a) ? (a)->n : 0, (i), \
				(val), err_ERROR) \
		 : (BIT_OFF != (val) \
		    ? (void) ((a)->bytesp[(i) / BITSET_BLOCK_SIZE] |= \
			(BITSET_BLOCK_TYPE) 1 \
				<< ((unsigned) (i) % BITSET_BLOCK_SIZE)) \
		    : (void) ((a)->bytesp[(i) / BITSET_BLOCK_SIZE] &= \
			~((BITSET_BLOCK_TYPE) 1 \
				<< ((unsigned) (i) % BITSET_BLOCK_SIZE))), \
		    err_OK))

/**
* Print the position of the bits in the given bitset whose value
* is 1.
* @param a the bitset to be printed.
*/
extern err_state print_bitset_states(/*@observer@*/ const bitset *)
		/*@modifies fileSystem@*/;

/**
* Fills the given bitset with one.
* @param a the bitset to be filled.
* @return returns a pointer to the result of the filling operation.
*/
extern err_state fill_bitset_one(bitset * a) /*@modifies *a@*/;

/**
* Fills the given bitset with zero.
* @param a the bitset to be filled.
* @return returns a pointer to the result of the filling operation.
*/
extern err_state fill_bitset_zero(bitset * a) /*@modifies *a@*/;

/**
* Checks if the bitset contains only zeros.
* @param a the bitset to be checked.
* @return 0: bitset contains a non-zero, 1: bitset contains only zeros.
*/
extern BOOL is_bitset_zero(/*@observer@*/ const bitset*) /*@modifies nothing@*/;

/**
* Get the Index of the next non-zero element.
* @param a the bitset to be checked.
* @param idx the Present index. If set to -1 then it is the first time call
*		for this method on the given bitset i.e. we have
*		to start searching from the very beginning
* @return the next index, or -1 if there are no more elements left
*/
extern state_index get_idx_next_non_zero(/*@observer@*/ const bitset *,
                state_index)
		/*@modifies nothing@*/;

/**
* Count the number of non-zero elements in the given bitset.
* @param a the bitset to be checked.
* @return the count.
*/
extern state_count count_non_zero(/*@observer@*/ const bitset *)
                /*@modifies nothing@*/;

/**
* Allocate an array that contains the numbers of all bits that are set in the
* bitset.
* the first element of the array is the count.
* the rest are the bit ids of set bits.
* @param	: bitset *toCount: the bitset in question
* @return	: the count of bits that are set and their ids.
*/
extern /*@only@*/ /*@null@*/ state_index * count_set(
                /*@observer@*/ const bitset * toCount)
                /*@modifies nothing@*/;

/**
* Returns the number of bits in the bitset.
*/
extern state_count bitset_size(
                /*@observer@*/ /*@temp@*/ /*@sef@*/ const bitset * bs)
                /*@modifies nothing@*/;
#define bitset_size(bs) (sizeof(bitset) != sizeof(*(bs)) || NULL == (bs) \
                ? err_macro_2(err_PARAM, "bitset_size(%p[%d])", \
                                (const void *) (bs), NULL!=(bs) \
                                ? ((void) (bs)->bytesp[0], (bs)->n) \
                                : 0, state_count_ERROR) \
                : (bs)->n)

#endif
