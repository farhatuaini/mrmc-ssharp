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

#include "bitset.h"

#include <string.h>

/* This macro comutes the number of bytes structure */
/* needed in order to hold the LENGTH bits */
state_count NUMBER_OF_BLOCKS(state_count LENGTH) /*@*/;
#define NUMBER_OF_BLOCKS(LENGTH) ((state_count) \
		((BITSET_BLOCK_SIZE - 1 + (LENGTH)) / BITSET_BLOCK_SIZE))

#define BYTE_OF_ONES ((char) ~'\0')
#define BLOCK_OF_ONES (~(BITSET_BLOCK_TYPE) 0)

/**
* Get a new bitset.
* A bitset is composed of an array of BITSET_BLOCK_SIZE-bit unsigned integers.
* @param n indicates the size of the required bitset.
* @return returns a pointer to a new bitset. NULL if some error has occurred.
*/
/*@null@*/ bitset * get_new_bitset(state_count n) {
	state_count d;
	/*@only@*/ /*@null@*/ bitset * newbitset = NULL;

        if ( n < 0 ) {
                err_msg_1(err_PARAM, "get_new_bitset(%d)", n, (bitset *) NULL);
        }

	d = NUMBER_OF_BLOCKS(n);
	/* Allocate the bitset */
        newbitset = (bitset *) calloc((size_t) 1, sizeof(bitset)
		- sizeof(BITSET_BLOCK_TYPE) + d * sizeof(BITSET_BLOCK_TYPE));

        if ( (bitset *) NULL == newbitset ) {
                err_msg_1(err_MEMORY, "get_new_bitset(%d)", n, (bitset *) NULL);
        }

	/* Init the other bitset fields */
	newbitset->n = n;

	return newbitset;
}


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
/*@null@*/ bitset * bs_extend(/*@only@*/ /*@null@*/ bitset * pBitset,
			state_count new_size, BOOL isFillNewBlocksWithOne)
{
	state_count old_size;
	state_count new_blocks_num, old_blocks_num;

        if ( (bitset *) NULL == pBitset ) {
		old_size = 0;
        } else {
		old_size = pBitset->n;
        }
	/* We check for the number of valid bits here but not the block_num */
	/* because we can have the same number of blocks but still less valid
	   elements */
        if ( old_size > new_size ) {
                /*@-mustfreeonly@*/
		err_msg_4(err_PARAM, "bs_extend(%p[%d],%d,%d)",
				(void *)pBitset, NULL!=pBitset ? pBitset->n : 0,
				new_size, (int) isFillNewBlocksWithOne,
                                (bitset *) NULL);
                /*@=mustfreeonly@*/
        }
	new_blocks_num = NUMBER_OF_BLOCKS(new_size);
	old_blocks_num = NUMBER_OF_BLOCKS(old_size);
	if ( new_blocks_num > old_blocks_num || (bitset *) NULL == pBitset ) {
		/* If pBitset == NULL, then we have to try and allocate memory
		   even if new_size == 0; returning NULL would lead the caller
		   to think the function has failed. */
		/* Add new blocks */
                /*@only@*/ /*@null@*/ bitset * pBitset_new = realloc(pBitset,
			sizeof(bitset) - sizeof(BITSET_BLOCK_TYPE) +
			new_blocks_num * sizeof(BITSET_BLOCK_TYPE));
                if ( NULL == pBitset_new ) {
                        /*@-usereleased@*/
			err_msg_4(err_MEMORY, "bs_extend(%p[%d],%d,%d)",
				(void *) pBitset, pBitset->n,
                                new_size, (int) isFillNewBlocksWithOne,
                                (bitset *) NULL);
                        /*@=usereleased@*/
                }
                pBitset = pBitset_new;
		/* Fill the new blocks with 0 or 1 as desired */
		memset(&pBitset->bytesp[old_blocks_num],
			isFillNewBlocksWithOne ? (int) BYTE_OF_ONES : 0,
			sizeof(BITSET_BLOCK_TYPE)
			* (new_blocks_num - old_blocks_num));
        /*@i1@*/ }
	/* Otherwise there is no need to reallocate the memory */

	/* Fill the new bits in the old block with 1 or set them to 0 */
	if ( old_size % BITSET_BLOCK_SIZE != 0 ) {
		BITSET_BLOCK_TYPE old_invalid_bit_mask = BLOCK_OF_ONES
				<< ((unsigned) old_size % BITSET_BLOCK_SIZE);
                if ( isFillNewBlocksWithOne ) {
			pBitset->bytesp[old_blocks_num - 1]
				|= old_invalid_bit_mask;
                } else {
			pBitset->bytesp[old_blocks_num - 1]
				&= ~old_invalid_bit_mask;
                }
	}

	/* set unused bits to 0 */
        if ( new_size % BITSET_BLOCK_SIZE != 0 ) {
		pBitset->bytesp[new_blocks_num - 1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) new_size % BITSET_BLOCK_SIZE));
        }
	pBitset->n = new_size;
	return pBitset;
}


/**
* This method just copies elements from one bitset into another.
* @param from_bitset the bitset to be copied .
* @param to_bitset the bitset to copy to.
* NOTE: The two bitsets must have the same size.
*/
err_state copy_bitset(/*@observer@*//*@i1@*//*@null@*/const bitset *from_bitset,
		/*@i1@*/ /*@null@*/ bitset * to_bitset)
{
	if ( (bitset *) NULL == from_bitset
				|| (bitset *) NULL == to_bitset
				|| from_bitset->n != to_bitset->n )
        {
		err_msg_4(err_PARAM, "copy_bitset(%p[%d],%p[%d])",
			(const void *) from_bitset,
				NULL != from_bitset ? from_bitset->n : 0,
			(void *) to_bitset, NULL!=to_bitset ? to_bitset->n : 0,
			err_ERROR);
        }

	/* copy the bit values */
	memcpy(to_bitset->bytesp, from_bitset->bytesp,
		NUMBER_OF_BLOCKS(to_bitset->n) * sizeof(BITSET_BLOCK_TYPE));

	/* set unused bits to 0 */
        if ( to_bitset->n % BITSET_BLOCK_SIZE != 0 ) {
		to_bitset->bytesp[(to_bitset->n /*- 1*/) / BITSET_BLOCK_SIZE] &=
			~(BLOCK_OF_ONES
			<< ((unsigned) to_bitset->n % BITSET_BLOCK_SIZE));
        }
	return err_OK;
}


/**
* Allocate a bitset and set it to the disjunction of two bitsets.
* @param a one of the operands to the disjunction operation.
* @param b the second operand to the disjunction operation.
* @return returns a pointer to a new bitset that contains the result of
*		disjunction. NULL if some error occurred, e. g. a and b have
*		different sizes.
*/
/*@null@*/ bitset * or(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a,
                /*@observer@*/
		/*@i1@*/ /*@null@*/ const bitset * b)
{
	state_count i;
	/*@only@*/ /*@null@*/ bitset * reslt = (bitset *) NULL;
	BITSET_BLOCK_TYPE * pres;
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa, * pb;

        if ( (bitset *) NULL == a || (bitset *) NULL == b || a->n != b->n ) {
		err_msg_4(err_PARAM, "or(%p[%d],%p[%d])", (const void *) a,
				NULL != a ? a->n : 0, (const void *) b,
                                NULL != b ? b->n : 0, (bitset *) NULL);
        }

	/* allocate space for the result */
	reslt = get_new_bitset(a->n);
        if ( (bitset *) NULL == reslt ) {
		err_msg_4(err_CALLBY, "or(%p[%d],%p[%d])", (const void *) a,
                                a->n, (const void *) b, b->n, (bitset *) NULL);
        }

	/* calculate the bit values of the result */
	pres = reslt->bytesp;
	pa = a->bytesp, pb = b->bytesp;
        for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres++ = *pa++ | *pb++;
        }

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return reslt;
}

/**
* Obtain the disjunction of two bitsets and put result in arg2.
* @param a one of the operands to the disjunction operation.
* @param b the second operand to the disjunction operation. Upon return, it will
*		contain the disjunction.
*/
err_state or_result(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a,
		/*@i1@*/ /*@null@*/ bitset * reslt)
{
	state_count i;
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa;
	BITSET_BLOCK_TYPE * pres;

	if ( (bitset *) NULL == a || (bitset *) NULL == reslt
			|| a->n != reslt->n )
        {
		err_msg_4(err_PARAM, "or_result(%p[%d],%p[%d])",
				(const void *) a, NULL != a ? a->n : 0,
				(void *) reslt, NULL != reslt ? reslt->n : 0,
				err_ERROR);
        }

	/* calculate the bit values of the result */
	pa = a->bytesp;
	pres = reslt->bytesp;
        for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres++ |= *pa++;
        }

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return err_OK;
}


/**
* Allocate a bitset and set it to the exclusive-or of two bitsets.
* @param a one of the operands to the exclusive-or operation.
* @param b the second operand to the exclusive-or operation.
* @return returns a pointer to a new bitset that contains the result of
*		exclusive-or. NULL if some error occurred, e. g. a and b have
*		different sizes.
*/
/*@null@*/ bitset * xor(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a,
                /*@observer@*/
		/*@i1@*/ /*@null@*/ const bitset * b)
{
	state_count i;
        /*@null@*/ /*@only@*/ bitset * reslt = (bitset *) NULL;
	BITSET_BLOCK_TYPE *pres;
        /*@observer@*/
	const BITSET_BLOCK_TYPE *pa, *pb;

        if ( (bitset *) NULL == a || (bitset *) NULL == b || a->n != b->n ) {
		err_msg_4(err_PARAM, "xor(%p[%d],%p[%d])", (const void *) a,
				NULL != a ? a->n : 0, (const void *) b,
                                NULL != b ? b->n : 0, (bitset *) NULL);
        }

	/* allocate space for the result */
	reslt = get_new_bitset(a->n);
        if ( (bitset *) NULL == reslt ) {
		err_msg_4(err_CALLBY, "xor(%p[%d],%p[%d])", (const void *) a,
                                a->n, (const void *) b, b->n, (bitset *) NULL);
        }

	/* calculate the bit values of the result */
	pres = reslt->bytesp;
	pa = a->bytesp, pb = b->bytesp;
        for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres++ = *pa++ ^ *pb++;
        }

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return reslt;
}

/**
* Obtain the exclusive-or of two bitsets and put result in arg2.
* @param a one of the operands to the exclusive-or operation.
* @param b the second operand to the exclusive-or operation. Upon return, it
*		will contain the exclusive-or.
*/
err_state xor_result(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a,
		/*@i1@*/ /*@null@*/ bitset * reslt)
{
	state_count i;
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa;
	BITSET_BLOCK_TYPE * pres;

	if ( (bitset *) NULL == a || (bitset *) NULL == reslt
			|| a->n != reslt->n )
        {
		err_msg_4(err_PARAM, "xor_result(%p[%d],%p[%d])",
				(const void *) a, NULL != a ? a->n : 0,
				(void *) reslt, NULL != reslt ? reslt->n : 0,
				err_ERROR);
        }

	/* calculate the bit values of the result */
	pa = a->bytesp;
	pres = reslt->bytesp;
        for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres++ ^= *pa++;
        }

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return err_OK;
}

/**
* Allocate a bitset and set it to the conjunction of two bitsets.
* @param a one of the operands to the conjunction operation.
* @param b the second operand to the conjunction operation.
* @return returns a pointer to a new bitset that contains the result of
*		conjunction. NULL if some error occurred, e. g. a and b have
*		different sizes.
*/
/*@null@*/ bitset * and(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a,
                /*@observer@*/
		/*@i1@*/ /*@null@*/ const bitset * b)
{
	state_count i;
	bitset * reslt;
	BITSET_BLOCK_TYPE * pres;
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa, * pb;

        if ( (bitset *) NULL == a || (bitset *) NULL == b || a->n != b->n ) {
		err_msg_4(err_PARAM, "and(%p[%d],%p[%d])", (const void *) a,
				NULL != a ? a->n : 0, (const void *) b,
                                NULL != b ? b->n : 0, (bitset *) NULL);
        }

	/* allocate space for the result */
	reslt = get_new_bitset(a->n);
        if ( (bitset *) NULL == reslt ) {
		err_msg_4(err_CALLBY, "and(%p[%d],%p[%d])", (const void *) a,
                                a->n, (const void *) b, b->n, (bitset *) NULL);
        }

	/* calculate the bit values of the result */
	pres = reslt->bytesp;
	pa = a->bytesp, pb = b->bytesp;
        for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres++ = *pa++ & *pb++;
        }

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return reslt;
}

/**
* Obtain the conjunction of two bitsets and put result in arg2.
* @param a one of the operands to the conjunction operation.
* @param b the second operand to the conjunction operation. Upon return, it will
*		contain the conjunction.
*/
err_state and_result(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a,
		/*@i1@*/ /*@null@*/ bitset * reslt)
{
	state_count i;
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa;
	BITSET_BLOCK_TYPE * pres;

	if ( (bitset *) NULL == a || (bitset *) NULL == reslt
				|| a->n != reslt->n )
        {
		err_msg_4(err_PARAM, "and_result(%p[%d],%p[%d])",
				(const void *) a, NULL != a ? a->n : 0,
				(void *) reslt, NULL != reslt ? reslt->n : 0,
				err_ERROR);
        }

	/* calculate the bit values of the result */
	pa = a->bytesp;
	pres = reslt->bytesp;
        for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres++ &= *pa++;
        }

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return err_OK;
}

/**
* Allocate a bitset and set it to the complement of a bitset.
* @param a the operand to the complement operation.
* @return returns a pointer to a new bitset that contains the result of
*		complementation. NULL if some error occurred.
*/
/*@null@*/ bitset * not(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a)
{
	state_count i;
	/*@only@*/ /*@null@*/ bitset * reslt = (bitset *) NULL;
	BITSET_BLOCK_TYPE * pres;
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa;

        if ( (bitset *) NULL == a ) {
		err_msg_2(err_PARAM, "not(%p[%d])", (const void *) a,
                                NULL != a ? a->n : 0, (bitset *) NULL);
        }

	/* allocate space for the result */
	reslt = get_new_bitset(a->n);
        if ( (bitset *) NULL == reslt ) {
		err_msg_2(err_CALLBY, "not(%p[%d])", (const void *) a,
                                a->n, (bitset *) NULL);
        }

	/* calculate the bit values of the result */
	pres = reslt->bytesp;
	pa = a->bytesp;
        for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres++ = ~*pa++;
        }

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return reslt;
}

/**
* Obtain the complement of a bitsets and put result in the argument.
* @param reslt the operand to the complement operation. Upon return, it will
*		contain the complement.
*/
err_state not_result(/*@i1@*/ /*@null@*/ bitset * reslt)
{
	state_count i;
	BITSET_BLOCK_TYPE * pres;

        if ( (bitset *) NULL == reslt ) {
		err_msg_2(err_PARAM, "not_result(%p[%d])", (void *) reslt,
				NULL != reslt ? reslt->n : 0, err_ERROR);
        }

	/* calculate the bit values of the result */
	pres = reslt->bytesp;
	for ( i = NUMBER_OF_BLOCKS(reslt->n) ; i > 0 ; --i ) {
		*pres = ~*pres;
		++pres;
	}

	/* set unused bits to 0 */
        if ( reslt->n % BITSET_BLOCK_SIZE != 0 ) {
		pres[-1] &= ~(BLOCK_OF_ONES
				<< ((unsigned) reslt->n % BITSET_BLOCK_SIZE));
        }
	return err_OK;
}


/**
* Print the position of the bits in the given bitset whose value
* is 1.
* @param a the bitset to be printed.
*/
err_state print_bitset_states(/*@observer@*//*@i1@*//*@null@*/ const bitset * a)
{
	BOOL bFirst;
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa;
	state_index i;
	BITSET_BLOCK_TYPE val;
	state_index idx;

        if ( (bitset *) NULL == a ) {
                err_msg_2(err_PARAM, "print_bitset_states(%p[%d])",
                                (const void *)a, NULL!=a ? a->n : 0, err_ERROR);
        }
	printf("{ ");
	bFirst = TRUE;
	pa = a->bytesp;
	for ( i = 0 ; i < NUMBER_OF_BLOCKS(a->n) ; ++i ) {
		val = *pa++;

		/* if this is the last block, then ignore unused bits */
                if ( a->n / (state_count) BITSET_BLOCK_SIZE == i ) {
			/* The above test is equivalent to:
			 * (a->n - 1) / BITSET_BLOCK_SIZE == i
			 * && a->n % BITSET_BLOCK_SIZE != 0 */
			val &= ~(BLOCK_OF_ONES
				<< ((unsigned) a->n % BITSET_BLOCK_SIZE));
                }

		/* check all bits in block i */
                for ( idx = 0 ; 0 != val ; val >>= 1, ++idx ) {
			if ( 0 != (val & 1) ) {
                                if ( ! bFirst ) {
					printf(", ");
                                } else {
					bFirst = FALSE;
                                }
				printf("%u", (unsigned)
					(i * BITSET_BLOCK_SIZE + idx + 1));
			}
                }
	}
	printf(" }");
	return err_OK;
}

/**
* Fills the given bitset with one.
* @param a the bitset to be filled.
* @return returns a pointer to the result of the filling operation.
*/
err_state fill_bitset_one(/*@i1@*/ /*@null@*/ bitset * a)
{
	/* Sets all bytes in the bitset to 255 (converted into a char).
	   It is equal to setting all bits to 1 as 255 is the greatest
	   value of unsigned char. */
        if ( (bitset *) NULL == a ) {
		err_msg_2(err_PARAM, "fill_bitset_one(%p[%d])", (const void *) a,
				NULL != a ? a->n : 0, err_ERROR);
        }
	memset(a->bytesp, (int) BYTE_OF_ONES,
			sizeof(BITSET_BLOCK_TYPE) * NUMBER_OF_BLOCKS(a->n));

	/* set unused bits to 0 */
        if ( a->n % BITSET_BLOCK_SIZE != 0 ) {
		a->bytesp[(a->n /*- 1*/) / BITSET_BLOCK_SIZE] &= ~(BLOCK_OF_ONES
				<< ((unsigned) a->n % BITSET_BLOCK_SIZE));
        }
	return err_OK;
}

/**
* Fills the given bitset with zero.
* @param a the bitset to be filled.
* @return returns a pointer to the result of the filling operation.
*/
err_state fill_bitset_zero(/*@i1@*/ /*@null@*/ bitset * a)
{
	/* Sets all bytes in the bitset to 0 (converted into an unsigned char)*/
        if ( (bitset *) NULL == a ) {
		err_msg_2(err_PARAM, "fill_bitset_zero(%p[%d])", (const void *) a,
				NULL != a ? a->n : 0, err_ERROR);
        }
	memset(a->bytesp, 0,
			sizeof(BITSET_BLOCK_TYPE) * NUMBER_OF_BLOCKS(a->n));
	/* set unused bits to 0 -- no separate action is needed, as this has
	   already been done by the above call to memset() */
	return err_OK;
}

/**
* Checks if the bitset contains only zeros.
* @param a the bitset to be checked.
* @return 0: bitset contains a non-zero, 1: bitset contains only zeros.
*/
BOOL is_bitset_zero(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a)
{
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa;
	state_count i;

        if ( (bitset *) NULL == a ) {
		err_msg_2(err_PARAM, "is_bitset_zero(%p[%d])", (const void *) a,
				NULL != a ? a->n : 0, (BOOL) -1);
        }
	/* check all blocks except the last one */
	pa = a->bytesp;
        for ( i = NUMBER_OF_BLOCKS(a->n) ; i > 1 ; --i ) {
                if ( 0 != *pa++ ) {
			return FALSE;
                }
        }

	/* check the last block, ignoring the unused bits */
	if ( 0 != (a->n % BITSET_BLOCK_SIZE != 0
			? *pa & ~(BLOCK_OF_ONES
				<< ((unsigned) a->n % BITSET_BLOCK_SIZE))
			: *pa) )
        {
		return FALSE;
        }
	return TRUE;
}

/**
* Get the Index of the next non-zero element.
* @param a the bitset to be checked.
* @param idx the Present index. If set to -1 then it is the first time call
*		for this method on the given bitset i.e. we have
*		to start searching from the very beginning
* @return the next index, or -1 if there are no more elements left
*/
state_index get_idx_next_non_zero(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a,
                state_index idx)
{
        /*@observer@*/
	const BITSET_BLOCK_TYPE *pa;
	BITSET_BLOCK_TYPE val;
	state_index i;

	if ( (bitset *) NULL == a || ((idx < 0 || idx >= a->n)
						&& idx != state_index_NONE) )
        {
		err_msg_3(err_PARAM, "get_idx_next_non_zero(%p[%d],%d)",
				(const void *) a, NULL != a ? a->n : 0, idx,
				state_index_ERROR);
        }
        if ( idx < 0 ) {
		idx = 0;
        } else {
		++idx;
        }
        if ( idx >= a->n ) {
		return state_index_NONE;
        }
	/* Now, idx == bit that has to be checked next */
	i = idx / (state_index) BITSET_BLOCK_SIZE;
	idx = idx % (state_index) BITSET_BLOCK_SIZE;
	/* Now, i*BITSET_BLOCK_SIZE + idx == bit that has to be checked next */
	pa = &a->bytesp[i];
	do {
		val = *pa++ >> (unsigned) idx;
		/* Now, LSB of val = bit that has to be checked next */
		if ( 0 != val ) {
			state_index result;
                        for ( ; 0 == (val & 1) ; val >>= 1, ++idx ) {
				;
                        }
			result = i * (state_index) BITSET_BLOCK_SIZE + idx;
                        if ( result >= a->n ) {
				/* if the result bit is too large, it is unused
				   and should be ignored */
				return state_index_NONE;
                        }
			return result;
		}
		++i;
		idx = 0;
	} while ( i < NUMBER_OF_BLOCKS(a->n) );
	return state_index_NONE;
}

/**
* Count the number of non-zero elements in the given bitset.
* @param a the bitset to be checked.
* @return the count.
*/
state_count count_non_zero(/*@observer@*/ /*@i1@*/ /*@null@*/ const bitset * a)
{
        /*@observer@*/
	const BITSET_BLOCK_TYPE * pa;
	state_index i;
	state_count count;

        if ( (bitset *) NULL == a ) {
		err_msg_2(err_PARAM, "count_non_zero(%p[%d])",
				(const void *) a, NULL != a ? a->n : 0,
				state_count_ERROR);
        }
	count = 0;

	pa = a->bytesp;
	for ( i = NUMBER_OF_BLOCKS(a->n) ; i > 0 ; --i ) {
		BITSET_BLOCK_TYPE val = *pa++;

		/* ignore unused bits, if this is the last block */
                if ( 1 == i && a->n % BITSET_BLOCK_SIZE != 0 ) {
			val &= ~(BLOCK_OF_ONES <<
					((unsigned) a->n % BITSET_BLOCK_SIZE));
                }
                for ( ; 0 != val ; val >>= 1 ) {
                        if ( 0 != (val & 1) ) {
				++count;
                        }
                }
	}
	return count;
}

/**
* Allocate an array that contains the numbers of all bits that are set in the
* bitset.
* the first element of the array is the count.
* the rest are the bit ids of set bits.
* @param	: bitset *toCount: the bitset in question
* @return	: the count of bits that are set and their ids.
*/
/*@only@*/ /*@null@*/ state_index * count_set(
                /*@observer@*/
		/*@i1@*/ /*@null@*/ const bitset *toCount)
{
        /*@observer@*/
	const BITSET_BLOCK_TYPE *pa;
	state_index i;
	/*@null@*/ /*@only@*/ state_index *pValidSet = (state_index *) NULL;

        if ( (bitset *) NULL == toCount ) {
		err_msg_2(err_PARAM, "count_set(%p[%d])",
				(const void *) toCount,
                                NULL != toCount ? toCount->n : 0,
                                (state_index *) NULL);
        }

	/* allocate an empty result set */
        pValidSet = (state_index *) calloc((size_t) 1, sizeof(state_index));
        if ( (state_index *) NULL == pValidSet ) {
		err_msg_2(err_MEMORY, "count_set(%p[%d])",
				(const void *) toCount,
                                toCount->n, (state_index *) NULL);
        }
	pa = toCount->bytesp;
	for ( i = 0 ; i < NUMBER_OF_BLOCKS(toCount->n) ; ++i ) {
		BITSET_BLOCK_TYPE val = *pa++;
		state_index idx;

		/* ignore unused bits, if this is the last block */
                if ( toCount->n / (state_count) BITSET_BLOCK_SIZE == i ) {
			/* The above condition is equivalent to:
			   (toCount->n - 1) / BITSET_BLOCK_SIZE == i
			   && toCount->n % BITSET_BLOCK_SIZE != 0 */
			val &= ~(BLOCK_OF_ONES <<
				((unsigned) toCount->n % BITSET_BLOCK_SIZE));
                }
                for ( idx = 0 ; 0 != val ; val >>= 1, ++idx ) {
			if ( 0 != (val & 1) ) {
				/* some set bit has been found -- append to the
				   result */
                                /*@only@*/
				state_index * temp;
				++pValidSet[0];
				temp = (state_index *) realloc(pValidSet,
						(pValidSet[0] + 1)
						* sizeof(state_index));
                                if ( (state_index *) NULL == temp ) {
                                        /*@-usereleased@*/
					err_msg_2(err_MEMORY,
						"count_set(%p[%d])",
						(const void *) toCount,
						NULL!=toCount ? toCount->n : 0,
                                                (/*@i1@*/free(pValidSet),
                                                 (state_index *) NULL));
                                }
				pValidSet = temp;
                                /*@=usereleased@*/
				pValidSet[pValidSet[0]] =
					i*(state_index)BITSET_BLOCK_SIZE + idx;
			}
                }
	}
	return pValidSet;
}
