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
*       Authors: Maneesh Khattri, Ivan Zapreev, David N. Jansen
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
*	Source description: Manage Sparse Matrix.
*	Uses: DEF: sparse.h
*	Definition of sparse matrix - sparse.h
*/

#ifndef SPARSE_H
#define SPARSE_H

#include "macro.h"
#include "error.h"

/*****************************************************************************
			STRUCTURE
name		: values
purpose         : information about a single state in a transition probability
                  matrix.
@member col	: column indices of non-zero elements.
@member val	: values of non-zero elements.
@member back_set: states from which this state is reachable.
@member ncols   : number of next-states, i.e. entries in the row, not counting
                  the diagonal
@member pcols   : number of previous-states, i.e. entries in the column, not
                  counting the diagonal.
@member diag	: The value of the diagonal element in this row.
remark		: There is a one-one relation between col and val. For instance
		  col[0] indicates the column index of val[0]. This structure
		  is similar to the compressed-row storage for sparse matrices.
******************************************************************************/
/* The following splint comments only hold for normal and contiguous (ncolse)
sparse matrices. In wrapper matrices, ``owned'' would have to be replaced by
``dependent''. However, this is not achievable without very deep changes of the
sparse matrix implementation. */
typedef struct values
{
        /*@owned@*/ /*@relnull@*/
	int *col;    /* the column ids */
                                        /* is required to be the first field in
                                           the structure. See
                                           allocate_sparse_matrix_ncolse(). */
        /*@owned@*/ /*@relnull@*/
	double *val; /* the corresponding values */
        /*@only@*/ /*@null@*/
	int *back_set;
                                        /* list of pointers to previous-states*/
        state_count ncols;              /* number of next-states, i.e. size of
                                           *col and *val */
        state_count pcols;              /* number of previous-states, i.e. size
                                           of *back_set */
	double diag;
                                        /* probability of self-loop */
}values;

/*****************************************************************************
			STRUCTURE
name			: sparse
@member rows            : number of rows in the matrix.
@member cols            : number of columns in the matrix.
@member valstruc        : this is a vector of rows, each element is a structure,
			containing values of non-zero elements.
remark			: There is a one-one relation between col and val. For instance
			col[0] indicates the column index of val[0]. This structure
			is similar to the compressed-row storage for sparse matrices.
******************************************************************************/
typedef /*@abstract@*/ struct sparse
{
	int rows;
                                        /* number of rows */
	int cols;
                                        /* numer of columns */
        struct values valstruc[1];      /* struct hack: actually a rows-sized
                                           array of information about the matrix
                                           rows */
}sparse;

/************************************************************************/
/**************      SIMPLE MATRIX INFORMATION ACCESS      **************/
/************************************************************************/
        /**
        * The four following functions provide some simple information on a
        * matrix.
        * mtx_rows(p_mtx) -- the number of rows of *p_mtx
        * mtx_cols(p_mtx) -- the number of columns of *p_mtx
        * mtx_next_num(p_mtx,row) -- the number of off-diagonal elements in the
        *                       given row of *p_mtx
        * mtx_prev_num(p_mtx,col) -- the number of off-diagonal elements in the
        *                       given column of *p_mtx
        */
        extern state_count mtx_rows(
                        /*@sef@*/
                        /*@observer@*/ /*@temp@*/ const sparse * p_mtx)
                        /*@modifies nothing@*/;
#       define mtx_rows(p_mtx) ((void) (p_mtx)->cols, \
                        (void) (p_mtx)->valstruc[0].back_set, \
                        (void) ((p_mtx)->valstruc[0].ncols + \
                                (p_mtx)->valstruc[0].pcols), \
                        (void) (p_mtx)->valstruc[0].diag, \
                        (const state_count) (p_mtx)->rows)
        extern state_count mtx_cols(
                        /*@sef@*/
                        /*@observer@*/ /*@temp@*/ const sparse * p_mtx)
                        /*@modifies nothing@*/;
#       define mtx_cols(p_mtx) ((void) (p_mtx)->rows, \
                        (void) (p_mtx)->valstruc[0].back_set, \
                        (void) ((p_mtx)->valstruc[0].ncols + \
                                (p_mtx)->valstruc[0].pcols), \
                        (void) (p_mtx)->valstruc[0].diag, \
                        (const state_count) (p_mtx)->cols)
        extern state_count mtx_next_num(
                        /*@sef@*/
                        /*@observer@*/ /*@temp@*/ const sparse * p_mtx,
                        state_index row) /*@modifies nothing@*/;
#       define mtx_next_num(p_mtx,row) ((void) (p_mtx)->rows, \
                        (void) (p_mtx)->cols, \
                        (void) (p_mtx)->valstruc[0].back_set, \
                        (void) (p_mtx)->valstruc[0].pcols, \
                        (void) (p_mtx)->valstruc[0].diag, \
                        (const state_count) (p_mtx)->valstruc[(row)].ncols)
        extern state_count mtx_prev_num(
                        /*@observer@*/ /*@temp@*/ const sparse * p_mtx,
                        state_index col);
#       define mtx_prev_num(p_mtx,col) ((void) (p_mtx)->rows, \
                        (void) (p_mtx)->cols, \
                        (void) (p_mtx)->valstruc[0].back_set, \
                        (void) (p_mtx)->valstruc[0].ncols, \
                        (void) (p_mtx)->valstruc[0].diag, \
                        (const state_count) (p_mtx)->valstruc[(col)].pcols)

/*======================================================================*/
/************************************************************************/
/************************A trivial MATRIX ALLOCATION*********************/
/*======================================================================*/

	/**
	* Allocates a sparse matrix
	* @param rows the number of rows in the matrix.
	* @param cols the number of cols in the matrix.
	* @return the pointer to the newly created sparse matrix
        *               or NULL if some error happened.
	*/
        extern /*@only@*/ /*@null@*/ sparse * allocate_sparse_matrix(int rows,
                        int cols) /*@modifies nothing@*/;


/*======================================================================*/
/************************************************************************/
/************************An alternative MATRIX ALLOCATION****************/
/*======================================================================*/
/**
* WARNING: The way this sparse matrix is allocated in
*       "allocate_sparse_matrix_ncolse" imposes some
*	strict conditions on some operations on this matrix, such as adding a new element
*/

	/**
	* Allocate a sparse matrix in a special way. Here we assume that we know the number of non-zero
	* off-diagonal elements in the matrix. And we also assume that the matrix non-diagonal values are
	* unlikely to change.
	* @param rows the number of rows in the matrix.
	* @param cols the number of cols in the matrix.
	* @param ncolse the array with the number of non-zero cols in each row.
	* @return a newly created sparse matrix
        *       or NULL if some error happened.
	* NOTE: Matrices which were allocated with this method are faster to
	*	work with, at the same time they should be treated with care
	*	because you can not just add - reallocate rows in them
	*	So there is a specific set of methods to be used for this.
	*/
        extern /*@only@*/ /*@null@*/ sparse * allocate_sparse_matrix_ncolse(
                        int rows, int cols, /*@observer@*/ const int * ncolse)
                        /*@modifies nothing@*/;

	/**
	* Frees the sparse matrix.
	* WARNING: This method should be used with the matrices allocated with
	*	the allocate_sparse_matrix_ncolse(...) method.
	* @param pM the matrix to remove
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state free_sparse_ncolse(/*@only@*/ sparse * pM)
                        /*@modifies pM@*/;

	/**
	* Sets a value to a pre-allocated element in the matrix.
	* WARNING: This method should be used with the matrices allocated with
	*	the allocate_sparse_matrix_ncolse(...) method.
	* @param pM the matrix to add element to
	* @param row the row number in which the value is to be set.
	* @param col the column number in which the value is to be set.
	* @param val the value to be set.
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state set_mtx_val_ncolse(sparse * pM, int row, int col,
                        double val) /*@modifies *pM@*/;

	/**
	* Searches for a pre-allocated matrix element and adds a value to it.
	* WARNING: This method should be used with the matrices allocated with
	*	the allocate_sparse_matrix_ncolse(...) method.
	* @param pM the matrix
	* @param row the row number in which the value is to be added.
	* @param col the column number in which the value is to be added.
	* @param val the value to be added.
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state add_mtx_val_ncolse(sparse * pM, int row, int col,
                        double val) /*@modifies *pM@*/;

/*======================================================================*/
/************************************************************************/
/************************General Sparse matrix methods*******************/
/**************For matrices allocated with "init_matrix"*****************/
/*  NOTE: NOT ALL of these methods are applicable to matrices allocated
	with the "allocate_sparse_matrix_ncolse" method! For example you can not
	add a new row element with reallocating the row elements array!  */
/*======================================================================*/

	/**
	* Creates new R matrix and copies all the rows that belong to the pValidStates into it
	* @param pStateSpace the initial matrix
	* @param pQ the new matrix to where to copy rows
	* @param pValidStates the valid states i.e. ids of the rows to be copied
	*                     this array contains the number of nodes as the first element
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state initMatrix(/*@observer@*/ const sparse * pStateSpace,
                        sparse * pQ, /*@observer@*/ const int * pValidStates)
                        /*@modifies *pQ@*/;

	/**
	* This method removes all the non zero rows from the Q matrix of the BSCC
	* @param pQ the BSCC nodes of the diagoriginal matrix
	* @param pValidStates the valid states i.e. ids of the rows which have to be cleaned
	*                     this array contains the number of nodes as the first element
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state cleanMatrix(sparse * pQ,
                        /*@observer@*/ const int * pValidStates)
                        /*@modifies *pQ@*/;

	/*****************************************************************************
	name		: free_mtx_sparse
	role		: frees the sparse matrix.
	@param		: sparse * pM the matrix to remove
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state free_mtx_sparse(/*@only@*/ sparse * pM)
                        /*@modifies pM@*/;

	/**
	* ADD a NEW value to the matrix.
	* WARNING: Should NOT be used with the matrices allocated with
	* allocate_sparse_matrix_ncolse(,,) method!
	* @param	: sparse * pM the matrix to add element to
	* @param	: int row: The row number in which the value is to be set.
	* @param	: int col: The column number in which the value is to be set.
	* @param	: double val: The value to be set.
        * @return       : err_ERROR: fail, err_OK: success
        * NOTE: return err_ERROR indicates that row/col does not meet the bounds
        *               of matrix
	*/
        extern err_state set_mtx_val(sparse * pM, int row, int col,
                        double val) /*@modifies *pM@*/;

	/*****************************************************************************
	name		: free_mtx_wrap
	role		: frees the sparse matrix. Does not clean rows.
	@param		: sparse * pM the matrix to remove
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state free_mtx_wrap(/*@only@*/sparse * pM) /*@modifies pM@*/;

	/**
	* Add a value to an existing one in the matrix or inserts a new value if the is none.
	* WARNING: Should NOT be used with the matrices allocated with
	* allocate_sparse_matrix_ncolse(,,) method! Unless you are 100% sure that
	* the pM[row,col] value already exists in the matrix!
	* @param	: sparse * pM the matrix
	* @param     : int row: The row number in which the value is to be added.
	* @param	: int col: The column number in which the value is to be added.
	* @param	: double val: to be added.
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state add_mtx_val(sparse * pM, int row, int col,
                        double val) /*@modifies *pM@*/;

	/*****************************************************************************
	name		: get_mtx_diag_val
	role		: get the diag value in the row of a matrix.
	@param		: sparse * pM the matrix to get element from
	@param		: int row: The row number in which the value is to be set.
	@param		: double * val (by reference): The value.
        @return         : err_ERROR: fail, err_OK: success
        remark          : return err_ERROR: row/col does not meet bounds.
	******************************************************************************/
        extern double mtx_get_diag_val_nt(
                        /*@temp@*/ /*@observer@*/ const sparse * pM,
                        int row) /*@modifies nothing@*/;
#       define mtx_get_diag_val_nt(pM,row) ((void)(pM)->rows, (void)(pM)->cols,\
                        (void) (pM)->valstruc[0].back_set, \
                        (void) ((pM)->valstruc[0].ncols + \
                                (pM)->valstruc[0].pcols), \
                        (const double) (pM)->valstruc[(row)].diag)

        extern err_state get_mtx_diag_val(
                        /*@sef@*/ /*@temp@*/ /*@observer@*/ const sparse * pM,
                        /*@sef@*/ int row, /*@sef@*/ /*@out@*/ double * value)
                        /*@modifies *value@*/;
#       define get_mtx_diag_val(pM,row,value) \
                ( NULL == (pM) || (unsigned) (row) >= (unsigned) mtx_rows((pM))\
                 || (unsigned)(row)>=(unsigned)mtx_cols((pM)) || NULL==(value) \
                 ? err_macro_3(err_PARAM, "get_mtx_diag_val(%p,%d,%p)", \
                        (const void *) (pM), (row), (void *) (value), \
                        (*(value)=0.0, err_ERROR)) \
                 : (*(value) = mtx_get_diag_val_nt((pM), (row)), err_OK) )

        /**
        * This function sets a value on the diagonal of pM.
        * @param pM     the matrix whose value has to be set
        * @param row    row number of the element which has to be set (= column
        *               number)
        * @param value  new value
        * @result       err_OK if everything went fine; err_ERROR otherwise.
        */
        extern void mtx_set_diag_val_nt(/*@temp@*/ sparse * pM, int row,
                        double value) /*@modifies *pM@*/;
#       define mtx_set_diag_val_nt(pM,row,value) ((void) (pM)->rows, \
                        (void) (pM)->cols, \
                        (void) (pM)->valstruc[0].back_set, \
                        (void) ((pM)->valstruc[0].ncols + \
                                (pM)->valstruc[0].pcols), \
                        (void) ((pM)->valstruc[(row)].diag = (value)))

        extern err_state mtx_set_diag_val(/*@temp@*/ /*@sef@*/ sparse * pM,
                        /*@sef@*/ int row, double value) /*@modifies *pM@*/;
#       define mtx_set_diag_val(pM,row,value) \
                ( NULL == (pM) || (unsigned) (row) >= (unsigned) mtx_rows((pM))\
                                || (unsigned)(row) >= (unsigned)mtx_cols((pM)) \
                 ? err_macro_3(err_PARAM, "set_mtx_diag_val(%p,%d,%g)", \
                                (void *) (pM), (row), (value), err_ERROR) \
                 : (mtx_set_diag_val_nt((pM), (row), (value)), err_OK) )

	/*****************************************************************************
	name		: get_mtx_val
	role		: get a value in the matrix.
	@param		: sparse * pM the matrix to get element from
	@param		: int row: The row number in which the value is to be set.
	@param		: int col: The column number in which the value is to be set.
	@param		: double * val (by reference): The value.
        @return         : err_ERROR: fail, err_OK: success
        remark          : If the value was not found, err_OK is returned.
	******************************************************************************/
        extern err_state get_mtx_val(/*@observer@*/ const sparse *pM,
                        int row, int col, /*@out@*/ double * value)
                        /*@modifies *value@*/;

	/*****************************************************************************
	name		: get_mtx_row_sums
	role		: add the elements in each row.
	@param		: sparse * pM the matrix to work with
	@return		: double *: row sums for each row.
                          or NULL if some error happened.
	******************************************************************************/
        extern /*@only@*/ /*@null@*/ double * get_mtx_row_sums(
                        /*@observer@*/ const sparse * pM)
                        /*@modifies nothing@*/;

	/*****************************************************************************
	name		: get_mtx_cer_row_sums
	role		: add elements in certain rows and return the vector of values.
	@param		: sparse * pM the matrix to work with
	@param		: int * pIds the ids of the rows for which the sums should be computed
				the pIds[0] contains the number of row ids
	@return		: the rows' sums
                          or NULL if some error happened.
	******************************************************************************/
        extern /*@only@*/ /*@null@*/ double * get_mtx_cer_row_sums(
                        /*@observer@*/ const sparse * pM,
                        /*@observer@*/ const int * pIds)
                        /*@modifies nothing@*/;

	/*****************************************************************************
	name		: mult_mtx_const
	role		: multiply all elements in the matrix with a constant.
	@param		: sparse * pM the matrix to multiply
	@param		: double constant: The constant value.
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state mult_mtx_const(sparse * pM, double constant)
                        /*@modifies *pM@*/;

	/*****************************************************************************
	name		: mult_mtx_cer_const
	role		: multiply all elements in the matrix with a constant.
	@param		: sparse * pM the matrix to multiply
	@param		: double constant: The constant value.
	@param		: int * pValidStates the valid rows to multiply pValidStates[0]
			  is the number of states (rows) in the array
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state mult_mtx_cer_const(sparse * pM, double constant,
                        /*@observer@*/ const int * pValidStates)
                        /*@modifies *pM@*/;

	/**
	* Multiply all elements in all i'th rows of the matrix by a constant[i] (or 1/constant[i] ).
	* @param pM the matrix to multiply
	* @param constant The constant value.
	* @param pValidStates the valid rows to multiply pValidStates[0]
	*			is the number of states (rows) in the array
	* @param invert if == 0 then nothing, otherwise we invert the elements of constant
	*			before the multiplication.
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state mult_mtx_cer_const_array(sparse * pM,
                        /*@observer@*/ const double * constant,
                        /*@observer@*/ const int * pValidStates, BOOL invert)
                        /*@modifies *pM@*/;

	/**
	* This method makes an "embedded" DTMC out of the rate matrix.
	* WARNING: The diagonal elements remain untouched, and therefore should not be
	*		considered in any further operations on the resulting matrix.
	* WARNING: We do not just make an embedded DTMC, but we also eliminate the self loops
	*		As if we would not have the rate matrix but the generator matrix as the original
	* WARNING: The self loops on states are not modified. Thus the proper way
	*		to identify an absorbing state is to make sure that there are
	*		no non-diagonal elements in the corresponding row.
	* @param pM the matrix to operate on
	* @param pRowSums the array of row sums of the rate matrix.
	* @param pValidStates the valid rows of the matrix to operate on.
	*		NOTE: The element pValidStates[0] should store
	*		the number of states (rows) in the array.
	* @return the array of exit rates for the pValidStates states
        *               or NULL if some error happened.
	*	NOTE: the array dimension corresponds to the state-space dimensions
	*/
        extern /*@only@*/ /*@null@*/ double *
                        make_embedded_dtmc_out_of_rate_mtx_vs(sparse * pM,
                        /*@observer@*/ const double * pRowSums,
                        /*@observer@*/ const int * pValidStates)
                        /*@modifies *pM@*/;

	/**
        * This method restores the rate matrix from the "embedded" DTMC obtained
	* with the "make_embedded_dtmc_out_of_rate_mtx" method.
	* WARNING: The diagonal elements remain untouched, and therefor
	* we expect them to be preserved from before the matrix was turned in to an embedded DTMC.
	* @param pM the matrix to operate on
	* @param pValidStates the valid rows of the matrix to operate on.
	*		NOTE: The element pValidStates[0] should store
	*		the number of states (rows) in the array.
	* @param pExitRates the array of exit rates returned by make_embedded_dtmc_out_of_rate_mtx
        * @return       : err_ERROR: fail, err_OK: success
	* WARNING: The pExitRates has to be obtained with the same pValidStates on the same pM matrix
	*		using the make_embedded_dtmc_out_of_rate_mtx method
	*/
        extern err_state restore_rate_mtx_out_of_embedded_dtmc_vs(sparse * pM,
                        /*@observer@*/ const int * pValidStates,
                        /*@observer@*/ const double * pExitRates)
                        /*@modifies *pM@*/;

	/**
	* This method makes an "embedded" DTMC out of the rate matrix.
	* WARNING: The diagonal elements remain untouched, and therefore should not be
	*		considered in any further operations on the resulting matrix.
	* WARNING: We do not just make an embedded DTMC, but we also eliminate the self loops
	*		As if we would not have the rate matrix but the generator matrix as the original
	* WARNING: The self loops on states are not modified. Thus the proper way
	*		to identify an absorbing state is to make sure that there are
	*		no non-diagonal elements in the corresponding row.
	* @param pM the matrix to operate on
	* @param pRowSums the array of row sums of the rate matrix.
	* @return the array of exit rates for all states
        *               or NULL if some error happened.
	*	NOTE: the array dimension corresponds to the state-space dimensions
	*/
        extern /*@only@*/ /*@null@*/ double *
                        make_embedded_dtmc_out_of_rate_mtx_all(sparse * pM,
                        /*@observer@*/ const double * pRowSums)
                        /*@modifies *pM@*/;

	/**
        * This method restores the rate matrix from the "embedded" DTMC obtained
	* with the "make_embedded_dtmc_out_of_rate_mtx_all" method.
	* WARNING: The diagonal elements remain untouched, and therefor
	* we expect them to be preserved from before the matrix was turned in to an embedded DTMC.
	* @param pM the matrix to operate on
	* @param pExitRates the array of exit rates returned by make_embedded_dtmc_out_of_rate_mtx
        * @return       : err_ERROR: fail, err_OK: success
	*/
        extern err_state restore_rate_mtx_out_of_embedded_dtmc_all(sparse * pM,
                        /*@observer@*/ const double * pExitRates)
                        /*@modifies *pM@*/;

	/*****************************************************************************
        name            : add_mtx_diagonal
	role		: add the elements of the input array to the diagonal one-one.
	@param		: sparse * pM the matrix to work with
	@param		: double *diagonal: The addendum to the diagonal.
        @return         : err_ERROR: fail, err_OK: success
        remark          : return err_ERROR: sparse matrix is not valid.
			  size should be correct.
	******************************************************************************/
        extern err_state add_mtx_diagonal(sparse * pM,
                        /*@observer@*/ const double * diagonal)
                        /*@modifies *pM@*/;

	/*****************************************************************************
        name            : sub_mtx_diagonal
        role            : sub the elements of the input array from the diagonal
                          one-one.
	@param		: sparse * pM the matrix to work with
        @param          : double *diagonal: The subtrahend from the diagonal.
        @return         : err_ERROR: fail, err_OK: success
        remark          : return err_ERROR: sparse matrix is not valid.
			  size should be correct.
	******************************************************************************/
        extern err_state sub_mtx_diagonal(sparse * pM,
                        /*@observer@*/ const double * diagonal)
                        /*@modifies *pM@*/;

	/****************************************************************************
        name            : add_cer_cons_diagonal_mtx
        role            : add the elements of the input array to the diagonal
                          one-one.
	@param		: sparse * pM the matrix to work with
	@param		: int * pIds the ids of the rows in which diagonal elements
                          should be increased by the corresponding pE values
			  the pIds[0] contains the number of row ids
	@param		: int constant: an operand (cons)
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state add_cer_cons_diagonal_mtx(sparse * pM,
                        /*@observer@*/ const int * pIds, int constant)
                        /*@modifies *pM@*/;

	/****************************************************************************
        name            : sub_cer_cons_diagonal_mtx
	role		: diag elements = cons - diag elements
	@param		: sparse * pM the matrix to work with
	@param		: int * pIds the ids of the rows in which diagonal elements
                          should be subtracted from the corresponding pE values
			  the pIds[0] contains the number of row ids
	@param		: int constant: an operand (cons)
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        /* extern err_state sub_cer_cons_diagonal_mtx(sparse * pM,
                        / *@observer@* / const int * pIds, int constant)
                        / *@modifies *pM@* /; */

	/*****************************************************************************
	name		: sub_mtx_cer_diagonal
	role		: Subtract elements of the pE array from the predefined diagonal
			  elements of the matrix
	@param		: sparse * pM the matrix to work with
	@param		: int * pIds the ids of the rows in which diagonal elements
                          should be reduced by the corresponding pE values
			  the pIds[0] contains the number of row ids
	@param		: double * pE: The elements to be subtracted
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state sub_mtx_cer_diagonal(sparse * pM,
                        /*@observer@*/ const int * pIds,
                        /*@observer@*/ const double * pE) /*@modifies *pM@*/;

	/*****************************************************************************
	name		: add_mtx_cons_diagonal
        role            : add a constant to all elements of the diagonal.
	@param		: sparse * pM the matrix to work with
	@param		: double constant: The constant value.
        @return         : err_ERROR: fail, err_OK: success
        remark          : return err_ERROR: sparse matrix is not valid.
	******************************************************************************/
        extern err_state add_mtx_cons_diagonal(sparse * pM, double cons)
                        /*@modifies *pM@*/;

	/*****************************************************************************
        name            : mtx_copy_row
        role            : copy a row from one matrix to another.
                          NOTE: This function copies the values structure. That
                          is why changes of fields in the original row later on
                          will not affect the one which was set by this
                          function. In particular, changing the diagonal element
			  in the original matrix will not affect the row added to the other matrix
			  by this method.
        @param          : sparse * to: the goal matrix to copy to
	@param     	: int row: The row number in which the value is to be set.
        @param          : sparse * from: the source matrix to copy from
        @return         : err_ERROR: fail, err_OK: success
        remark          : return err_ERROR: indicates that row/col does not meet
                          the bounds
			  of matrix
	******************************************************************************/
        extern err_state mtx_copy_row(sparse * to, int row,
                        /*@observer@*/ const sparse * from)
                        /*@modifies *to@*/;

	/*****************************************************************************
	name		: multiply_mtx_MV
	role		: multiply a matrix with a vector.
	@param		: sparse * pM: operand matrix.
	@param		: double *vec: The operand vector.
	@param		: double *res: The resulting vector.
        @return         : err_ERROR: fail, err_OK: success
	remark		: size should be correct.
	******************************************************************************/
        extern err_state multiply_mtx_MV(/*@observer@*/ const sparse * pM,
                        /*@observer@*/ const double * vec,
                        /*@out@*/ double * res) /*@modifies *res@*/;

	/*****************************************************************************
	name		: multiply_mtx_cer_MV
	role		: multiply certain rows of a matrix with a vector.
	@param		: sparse * pM: operand matrix.
	@param		: double *vec: The operand vector.
	@param		: double *res: The resulting vector.
	@param		: int num: number of valid rows.
	@param		: int *valid_rows: the valid rows
        @return         : err_ERROR: fail, err_OK: success
	remark		: size should be correct.
	******************************************************************************/
        extern err_state multiply_mtx_cer_MV(/*@observer@*/ const sparse * pM,
                        /*@observer@*/ const double * vec,
                        /*@out@*/ double * res, int num,
                        /*@observer@*/ const int * valid_rows)
                        /*@modifies *res@*/;

	/*****************************************************************************
	name		: multiply_mtx_TMV
	role		: multiply a vector to a matrix.
	@param		: sparse * pM the matrix to be used
	@param		: double *vec: The operand vector.
	@param		: double *res: The resulting vector.
        @return         : err_ERROR: fail, err_OK: success
	remark		: size should be correct.
	******************************************************************************/
        extern err_state multiply_mtx_TMV(/*@observer@*/ const sparse * pM,
                        /*@observer@*/ const double * vec,
                        /*@out@*/ double * res) /*@modifies *res@*/;

	/*****************************************************************************
	name		: multiply_mtx_cer_TMV
	role		: multiply a vector to a matrix.
	@param		: sparse * pM the matrix to be used
        @param          : double * vec: The operand vector.
        @param          : double * res: The resulting vector.
	@param		: int * pIds the valid ids, pIds[0] contains the number of ids
        @return         : err_ERROR: fail, err_OK: success
	remark		: size should be correct.
	******************************************************************************/
        extern err_state multiply_mtx_cer_TMV(/*@observer@*/ const sparse * pM,
                        /*@observer@*/ const double * vec,
                        /*@out@*/ double * res, /*@observer@*/ const int * pIds)
                        /*@modifies *res@*/;

/************************************************************************/
/**************SUPPLEMENTARTY METHODS USED FOR PRINTING ETC**************/
/************************************************************************/

        extern err_state print_pattern_vec_double(
                        /*@observer@*/ const char * pattern, int length,
                        /*@observer@*/ /*@null@*/ const double * vec)
                        /*@modifies fileSystem@*/;

        extern err_state print_vec_double(int length,
                        /*@observer@*/ /*@null@*/ const double * vec)
                        /*@modifies fileSystem@*/;

        extern err_state print_vec_int(int length,
                        /*@observer@*/ const int * vec)
                        /*@modifies fileSystem@*/;

        /* extern err_state print_vec_short(int length,
                        / *@observer@* / / *@null@* / const short * vec)
                        / *@modifies fileSystem@* /; */

	/**
	* Write the matrix to a file in tra-format
	* @param pM the matrix to work with
	* @param fname: the name of the output file
        * @return       : err_ERROR: fail, err_OK: success
        * NOTE: This function returns err_ERROR if a sparse matrix is not valid,
        *       or the file
	*	could not be opened for writing. This is also indicated via a printf statement.
	*/
        /* extern err_state write_tra_file(/ *@observer@* / const sparse * pM,
                        / *@observer@* / const char * fname)
                        / *@modifies fileSystem@* /; */

	/*****************************************************************************
        name            : print_mtx_sparse
	role		: print the matrix.
	@param		: sparse * pM the matrix to work with
        @return         : err_ERROR: fail, err_OK: success
        remark          : return err_ERROR: sparse matrix is not valid.
	******************************************************************************/
        extern err_state print_mtx_sparse(/*@observer@*/ const sparse * pM)
                        /*@modifies fileSystem@*/;

	/*****************************************************************************
	name		: split_A_into_DI_LU
	role		: this method splits the given A matrix into two matrixes DI and L+U
			: where DI is an inverted diagonal matrix and L+U are all the rest elements of
			  the initial matrix A.
	@param		: sparse *pA the initial matrix A
	@param		: sparse *pDI the empty matrix which will be filled with the D elements
	@param		: sparse *pLU the empty matrix which will be filled with the LU elements
        @return         : err_ERROR: fail, err_OK: success
	*****************************************************************************/
        extern err_state split_A_into_DI_LU(/*@observer@*/ const sparse * pA,
                        sparse * pDI, sparse * pLU) /*@modifies *pDI, *pLU@*/;

	/*****************************************************************************
	name		: split_A_into_D_LU
	role		: this method splits the given A matrix into two matrixes DI and L+U
			: where D is a diagonal matrix and L+U are all the rest elements of
			  the initial matrix A.
	@param		: sparse *pA the initial matrix A
	@param		: double *pD the empty vector which will be filled with the diagonal elements
	@param		: sparse *pLU the empty matrix which will be filled with the LU elements
        @return         : err_ERROR: fail, err_OK: success
	*****************************************************************************/
        extern err_state split_A_into_D_LU(/*@observer@*/ const sparse * pA,
                        /*@out@*/ double * pD, sparse * pLU)
                        /*@modifies *pD, *pLU@*/;

	/*****************************************************************************
	name		: split_A_into_D_L_U
	role		: This method splits the given A matrix into two matrixes DI and L+U
			: where D is a diagonal matrix and L+U are all the rest elements of
			  the initial matrix A.
			  NOTE: This method presumes that row elements are ordered by column id
			  and that is why it does not copy the pA elements to pL, pU but
			  simply copy pointers. Take this in mind while deleting pL and pU
	@param		: sparse *pA the initial matrix A
	@param		: double *pD the empty vector which will be filled with the diagonal elements
	@param		: sparse *pL the empty matrix which will be filled with the L elements
	@param		: sparse *pU the empty matrix which will be filled with the U elements
        @return         : err_ERROR: fail, err_OK: success
	*****************************************************************************/
        extern err_state split_A_into_D_L_U(/*@observer@*/ const sparse * pA,
                        /*@out@*/ double * pD, sparse * pL, sparse * pU)
                        /*@modifies *pD, *pL, *pU@*/;

	/*****************************************************************************
	name		: sum_vec_vec
	role		: This method is used to compute the sum of two vectors A+B=Res.
	@param		: int length the length of vectors
	@param		: double *pA vector A
	@param		: double *pB vector B
	@param		: double *pRes vector Res
        @return         : err_ERROR: fail, err_OK: success
	*****************************************************************************/
        extern err_state sum_vec_vec(int length,
                        /*@observer@*/ const double * pA,
                        /*@observer@*/ const double * pB,
                        /*@out@*/ double * pRes) /*@modifies *pRes@*/;

	/*****************************************************************************
	name		: sum_cer_vec_vec
	role		: This method is used to compute the sum of two vectors A+B=Res.
	@param		: double *pA vector A
	@param		: double *pB vector B
	@param		: double *pRes vector Res
	@param		: int * pIds the valid ids, pIds[0] contains the number of ids
        @return         : err_ERROR: fail, err_OK: success
	*****************************************************************************/
        extern err_state sum_cer_vec_vec(/*@observer@*/ const double * pA,
                        /*@observer@*/ const double * pB,
                        /*@out@*/ double * pRes,
                        /*@observer@*/ const int * pIds) /*@modifies *pRes@*/;

	/*****************************************************************************
	name		: multiply_inv_diag_D_V
	role		: multiply v*diag(d)I. I.e. take vector of diagonal elements
			  create matrix diag(d) invert it and multiply by the vector v
			  from the left side.
	@param		: double * pD: The diagonal elements.
	@param		: double * pV: The vector to multiply with.
	@param		: double * pR: The resulting vector.
	@param		: int length : The length of vectors
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state multiply_inv_diag_D_V(/*@observer@*/ const double * pD,
                        /*@observer@*/ const double * pV, /*@out@*/ double * pR,
                        int length) /*@modifies *pR@*/;

	/*****************************************************************************
	name		: multiply_cer_inv_diag_D_V
	role		: multiply v*diag(d)I. I.e. take vector of diagonal elements
			  create matrix diag(d) invert it and multiply by the vector v
			  from the left side.
	@param		: double * pD: The diagonal elements.
	@param		: double * pV: The vector to multiply with.
	@param		: double * pR: The resulting vector.
	@param		: int * pValidStates : The set of valid states
        @return         : err_ERROR: fail, err_OK: success
	******************************************************************************/
        extern err_state multiply_cer_inv_diag_D_V(
                        /*@observer@*/ const double * pD,
                        /*@observer@*/ const double * pV, /*@out@*/ double * pR,
                        /*@observer@*/ const int * pValidStates)
                        /*@modifies *pR@*/;

/************************************************************************/
/**************              MATRIX ITERATORS              **************/
/************************************************************************/

        /* Usage: To visit all elements in a given row, write:
                const sparse * matrix = ...;
                state_index row = ...;
                mtx_walk_row(matrix,row,col,val) {
                        statement to be executed on element (row,col);
                } end_mtx_walk_row;
        To visit all elements in column 123, write:
                const sparse * matrix = ...;
                mtx_walk_column(matrix,row,123,val) {
                        statement1;
                        statement2;
                } end_mtx_walk_column;
        If only the position and not the value is relevant, there is a more
        efficient way to visit all elements in a column using a different macro:
                const sparse * matrix = ...;
                state_index row = ...;
                mtx_walk_column_noval(matrix,row,col) {
                        statement;
                } end_mtx_walk_column_noval;
        The iterator macro defines the variables col or row, respectively, and
        val. They should not be defined outside the iterator scope.

        The iterators work like a loop, so it is possible to use continue and
        break. */

        /*@iter mtx_walk_row(sef observer const sparse * p_mtx,
                        sef state_index row, yield state_index m_col,
                        yield double m_val)@*/
#       define mtx_walk_row(p_mtx,row,m_col,m_val) \
                { \
                        state_index m_col = (row); \
                        double m_val = 0.0; \
                        state_count m_i__row; \
                        const int * m_colrow; \
                        const double * m_valrow; \
                        if ( NULL == (p_mtx) || (unsigned) m_col >= \
                                                (unsigned) mtx_rows((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_walk_row(%p," \
                                        "%d,col,val)", (const void *) (p_mtx), \
                                        m_col, EXIT_FAILURE)); \
                        } \
                        m_colrow = (p_mtx)->valstruc[m_col].col; \
                        m_valrow = (p_mtx)->valstruc[m_col].val; \
                        m_i__row = mtx_next_num((p_mtx), m_col) + 1; \
                        if ( m_col >= mtx_cols((p_mtx)) || \
                                        (m_val = mtx_get_diag_val_nt((p_mtx), \
                                         m_col)) == 0.0 ) \
                        { \
                                /* The following expression should be as */ \
                                /* similar as possible to the iteration */ \
                                /* expression in the for statement; this */ \
                                /* helps the compiler to optimize the code. */ \
                                (void) (0 < --m_i__row && \
                                        (m_col=*m_colrow++, m_val=*m_valrow++, \
                                         TRUE)); \
                        } \
                        for ( ; 0 < m_i__row ; \
                                (void) (0 < --m_i__row && \
                                        (m_col=*m_colrow++, m_val=*m_valrow++, \
                                         TRUE)) ) \
                        {

#       define end_mtx_walk_row \
                        } \
                        /* The following superfluous statement makes sure */ \
                        /* that the iterators are not mixed (all m_i__... */ \
                        /* variables have distinct names). The compiler */ \
                        /* should optimize the statement away. */ \
                        /*@-noeffect@*/ (void) m_i__row; /*@=noeffect@*/ \
                 }
        
        /* The slightly awkward handling of m_i (in particular, using m_i-1 in
           the step expression of for) is done to make sure that the macro also
           works, should one decide to make the type state_count unsigned to
           handle even larger state spaces. */

        /*@iter mtx_walk_row_nodiag(sef observer const sparse * p_mtx,
                        sef state_index row, yield state_index m_col,
                        yield double m_val)@*/
#       define mtx_walk_row_nodiag(p_mtx,row,m_col,m_val) \
                { \
                        state_index m_col = (row); \
                        state_count m_i__row_nodiag; \
                        const int * m_colrow; \
                        const double * m_valrow; \
                        if ( NULL == (p_mtx) || (unsigned) m_col >= \
                                                (unsigned) mtx_rows((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_walk_row_" \
                                        "nodiag(%p,%d,col,val)", (const void *)\
                                        (p_mtx), m_col, EXIT_FAILURE)); \
                        } \
                        m_colrow = (p_mtx)->valstruc[m_col].col; \
                        m_valrow = (p_mtx)->valstruc[m_col].val; \
                        for ( m_i__row_nodiag = mtx_next_num((p_mtx), m_col) ; \
                                        0 < m_i__row_nodiag-- ; \
                                                m_colrow++, m_valrow++ ) \
                        { \
                                const double m_val = *m_valrow; \
                                m_col = *m_colrow;

#       define end_mtx_walk_row_nodiag \
                        } \
                        /*@-noeffect@*/ (void) m_i__row_nodiag; /*@=noeffect@*/\
                }

        /*@iter mtx_walk_row_sorted(sef observer const sparse * p_mtx,
                        sef state_index row, yield state_index m_col,
                        yield double m_val)@*/
#       define mtx_walk_row_sorted(p_mtx,row,m_col,m_val) \
                { \
                        state_index m_col = (row); \
                        state_count m_i__row_sorted; \
                        const int * m_colrow; \
                        const double * m_valrow; \
                        BOOL m_diag_passed = FALSE; \
                        if ( NULL == (p_mtx) || (unsigned) m_col >= \
                                                (unsigned) mtx_rows((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_walk_row_" \
                                        "sorted(%p,%d,col,val)", (const void *)\
                                        (p_mtx), m_col, EXIT_FAILURE)); \
                        } \
                        m_colrow = (p_mtx)->valstruc[m_col].col; \
                        m_valrow = (p_mtx)->valstruc[m_col].val; \
                        for ( m_i__row_sorted = mtx_next_num((p_mtx), m_col) ; \
                                                TRUE ; ) \
                        { \
                                double m_val; \
                                /* check if the diagonal */ \
                                /* element is the next one to be read */ \
                                if ( /* surely not if it's already been read */\
                                        m_diag_passed || \
                                     /* surely not if there is another */ \
                                     /* element left of the diagonal */ \
                                        (0 < m_i__row_sorted \
                                         && *m_colrow <= (row)) || \
                                     /* yes, at least in principle */ \
                                        (m_diag_passed = TRUE, \
                                     /* except if the diagonal element */ \
                                     /* doesn't exist */ \
                                         (m_col=(row)) >= mtx_cols((p_mtx))) ||\
                                     /* and except if the diagonal element */ \
                                     /* is zero */ \
                                        (m_val = mtx_get_diag_val_nt((p_mtx), \
                                                        m_col)) == 0.0 ) \
                                { \
                                        /* check whether there is an */ \
                                        /* off-diagonal element left */ \
                                        if ( 0 >= m_i__row_sorted ) \
                                                break; \
                                        /* yes, there is: */ \
                                        m_col = *m_colrow++; \
                                        m_val = *m_valrow++; \
                                        m_i__row_sorted--; \
                                }

#       define end_mtx_walk_row_sorted \
                        } \
                        /*@-noeffect@*/ (void) m_i__row_sorted; /*@=noeffect@*/\
                }

        /*@iter mtx_walk_column(sef observer const sparse * p_mtx,
                        yield state_index m_row, sef unused state_index clm,
                        yield double m_val)@*/
#       define mtx_walk_column(p_mtx,m_row,clm,m_val) \
                { \
                        state_index m_row = (clm); \
                        state_count m_i__column; \
                        const int * m_rows; \
                        double m_val = 0.0; \
                        if ( NULL == (p_mtx) || (unsigned) m_row >= \
                                                (unsigned) mtx_cols((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_walk_column(" \
                                        "%p,row,%d,val)", (const void*)(p_mtx),\
                                        m_row, EXIT_FAILURE)); \
                        } \
                        m_i__column = mtx_prev_num((p_mtx), m_row) + 1; \
                        m_rows = (p_mtx)->valstruc[m_row].back_set; \
                        if ( m_row >= mtx_rows((p_mtx)) || \
                                (m_val = mtx_get_diag_val_nt((p_mtx), m_row)) \
                                        == 0.0 ) \
                        { \
                                (void) (0 < --m_i__column \
                                && (m_row = *m_rows++, err_state_iserror( \
                                    get_mtx_val((p_mtx), m_row,(clm), &m_val)))\
                                && (exit(err_macro_2(err_CALLBY, "mtx_walk_" \
                                        "column(%p,row,%d,val)", (const void *)\
                                        (p_mtx), (clm), EXIT_FAILURE)), TRUE));\
                        } \
                        for ( ; 0 < m_i__column ; \
                                (void) (0 < --m_i__column \
                                && (m_row = *m_rows++, err_state_iserror( \
                                    get_mtx_val((p_mtx), m_row,(clm), &m_val)))\
                                && (exit(err_macro_2(err_CALLBY, "mtx_walk_" \
                                        "column(%p,row,%d,val)", (const void *)\
                                        (p_mtx),(clm), EXIT_FAILURE)), TRUE)) )\
                        {

#       define end_mtx_walk_column \
                        } \
                        /*@-noeffect@*/ (void) m_i__column; /*@=noeffect@*/ \
                }

        /*@iter mtx_walk_column_nodiag(sef observer const sparse * p_mtx,
                        yield state_index m_row, sef unused state_index clm,
                        yield double m_val)@*/
#       define mtx_walk_column_nodiag(p_mtx,m_row,clm,m_val) \
                { \
                        state_index m_row = (clm); \
                        state_count m_i__column_nodiag; \
                        const int * m_rows; \
                        if ( NULL == (p_mtx) || (unsigned) m_row >= \
                                                (unsigned) mtx_cols((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_walk_column_" \
                                        "nodiag(%p,row,%d,val)", (const void *)\
                                        (p_mtx), m_row, EXIT_FAILURE)); \
                        } \
                        m_rows = (p_mtx)->valstruc[m_row].back_set; \
                        for ( m_i__column_nodiag=mtx_prev_num((p_mtx), m_row) ;\
                                0 < m_i__column_nodiag ; m_rows++, \
                                                m_i__column_nodiag-- ) \
                        { \
                                double m_val; \
                                m_row = *m_rows; \
                                if ( err_state_iserror(get_mtx_val((p_mtx), \
                                                m_row, (clm), &m_val)) ) \
                                { \
                                        exit(err_macro_2(err_CALLBY, "mtx_walk"\
                                                "_column_nodiag(%p,row,%d,val)"\
                                                , (const void *)(p_mtx), (clm),\
                                                EXIT_FAILURE)); \
                                }

#       define end_mtx_walk_column_nodiag \
                        } \
                        /*@-noeffect@*/(void)m_i__column_nodiag;/*@=noeffect@*/\
                }

        /* The _noval variant of mtx_walk_column is more efficient because the
           call to get_mtx_val takes some time; if the value is not used, it is
           not necessary to spend that. */
        /*@iter mtx_walk_column_noval(sef observer const sparse * p_mtx,
                        yield state_index m_row, sef unused state_index clm)@*/
#       define mtx_walk_column_noval(p_mtx,m_row,clm) \
                { \
                        state_index m_row = (clm); \
                        state_count m_i__column_noval; \
                        const int * m_rows; \
                        if ( NULL == (p_mtx) || (unsigned) m_row >= \
                                                (unsigned) mtx_cols((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_walk_column_" \
                                        "noval(%p,row,%d)", (const void *) \
                                        (p_mtx), m_row, EXIT_FAILURE)); \
                        } \
                        m_i__column_noval = mtx_prev_num((p_mtx), m_row) + 1; \
                        m_rows = (p_mtx)->valstruc[m_row].back_set; \
                        if ( m_row >= mtx_rows((p_mtx)) || \
                                0.0 == mtx_get_diag_val_nt((p_mtx), m_row) ) \
                        { \
                                (void) (0 < --m_i__column_noval \
                                        && (m_row = *m_rows++, TRUE)); \
                        } \
                        for ( ; 0 < m_i__column_noval ; \
                                (void) (0 < --m_i__column_noval \
                                        && (m_row = *m_rows++, TRUE)) ) \
                        {

#       define end_mtx_walk_column_noval \
                        } \
                        /*@-noeffect@*/(void) m_i__column_noval;/*@=noeffect@*/\
                }

        /*@iter mtx_walk_column_nodiag_noval(sef observer const sparse * p_mtx,
                        yield state_index m_row, sef unused state_index clm)@*/
#       define mtx_walk_column_nodiag_noval(p_mtx,m_row,clm) \
                { \
                        state_index m_row = (clm); \
                        state_count m_i__column_nodiag_noval; \
                        const int * m_rows; \
                        if ( NULL == (p_mtx) || (unsigned) m_row >= \
                                                (unsigned) mtx_cols((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_walk_column_" \
                                        "nodiag_noval(%p,row,%d)",(const void*)\
                                        (p_mtx), m_row, EXIT_FAILURE)); \
                        } \
                        m_rows = (p_mtx)->valstruc[m_row].back_set; \
                        for ( m_i__column_nodiag_noval = mtx_prev_num((p_mtx), \
                                        m_row) ; \
                                0 < m_i__column_nodiag_noval ; m_rows++, \
                                                m_i__column_nodiag_noval-- ) \
                        { \
                                m_row = *m_rows;

#       define end_mtx_walk_column_nodiag_noval \
                        } \
                        /*@-noeffect@*/ (void) m_i__column_nodiag_noval; \
                        /*@=noeffect@*/ \
                }

        /*@iter mtx_walk_all(sef observer const sparse * p_mtx,
                        yield state_index m_row, yield state_index m_col,
                        yield double m_val)@*/
#       define mtx_walk_all(p_mtx,m_row,m_col,m_val) \
                { \
                        state_index m_row, m_col; \
                        double m_val; \
                        state_count m_i__all; \
                        const struct values * m_values; \
                        const int * m_colrow; \
                        const double * m_valrow; \
                        if ( NULL == (p_mtx) ) { \
                                exit(err_macro_1(err_PARAM, "mtx_walk_all(%p," \
                                        "row,col,val)", (const void *) (p_mtx),\
                                        EXIT_FAILURE)); \
                        } \
                        m_row = mtx_rows((p_mtx)); \
                        m_values = &(p_mtx)->valstruc[m_row]; \
                        /* The following lines should be as similar as */ \
                        /* possible to a part of the iteration expression */ \
                        /* in the for statement; this helps the compiler to */ \
                        /* optimize the code. Note, in particular, the use */ \
                        /* of commas instead of semicola. */ \
                        (void) (m_col = --m_row, --m_values, \
                        m_i__all = m_values->ncols + 1, \
                        m_colrow = m_values->col, m_valrow = m_values->val, \
                        m_val = m_col<mtx_cols((p_mtx)) ? m_values->diag : 0.0,\
                        TRUE); \
                        for( ; 0 < m_i__all ; \
                                0 >= --m_i__all \
                                     ? (void) (0 < m_row && \
                                       (m_col = --m_row, --m_values, \
                                        m_i__all = m_values->ncols + 1, \
                                        m_colrow = m_values->col, \
                                        m_valrow = m_values->val, \
                                        m_val = m_col < mtx_cols((p_mtx)) \
                                                        ? m_values->diag : 0.0,\
                                        TRUE)) \
                                     : (void) (m_col = *m_colrow++, \
                                        m_val = *m_valrow++) ) \
                        { \
                                /*@-realcompare@*/ \
                                if ( 0.0 != m_val ) /*@=realcompare@*/ {

#       define end_mtx_walk_all } \
                        } \
                        /*@-noeffect@*/ (void) m_i__all; /*@=noeffect@*/ \
                }

#endif

