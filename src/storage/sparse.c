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
*	Source description: Manage Sparse Matrix.
*	Uses: DEF: sparse.h
*	Definition of sparse matrix - sparse.h
*/

#include "sparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* Matrix iterators for internal use: iterators that change the matrix */
/*@iter mtx_change_row(sef sparse * p_mtx, sef state_index row,
                yield state_index m_col, yield double * m_p_val)@*/
#define mtx_change_row(p_mtx,row,m_col,m_p_val) \
                { \
                        state_index m_col = (row); \
                        double * m_p_val; \
                        state_count m_i__c_row; \
                        const int * m_colrow; \
                        double * m_valrow; \
                        if ( NULL == (p_mtx) || (unsigned) m_col >= \
                                                (unsigned) mtx_rows((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_change_row(" \
                                        "%p,%d,col,val)", (void *) (p_mtx), \
                                        m_col, EXIT_FAILURE)); \
                        } \
                        m_colrow = (p_mtx)->valstruc[m_col].col; \
                        m_valrow = (p_mtx)->valstruc[m_col].val; \
                        m_i__c_row = mtx_next_num((p_mtx), m_col) + 1; \
                        m_p_val = &(p_mtx)->valstruc[m_col].diag; \
                        if ( m_col >= mtx_cols((p_mtx)) ) { \
                                (void) (0 < --m_i__c_row && \
                                        (m_col=*m_colrow++, m_p_val=m_valrow++,\
                                         TRUE)); \
                        } \
                        for ( ; 0 < m_i__c_row ; \
                                (void) (0 < --m_i__c_row && \
                                        (m_col=*m_colrow++, m_p_val=m_valrow++,\
                                         TRUE)) ) \
                        {

#define end_mtx_change_row \
                        } \
                        /*@-noeffect@*/ (void) m_i__c_row; /*@=noeffect@*/ \
               }

/*@iter mtx_change_row_nodiag(sef sparse * p_mtx, sef state_index row,
                yield state_index m_col, yield double * m_p_val)@*/
#define mtx_change_row_nodiag(p_mtx,row,m_col,m_p_val) \
                { \
                        state_index m_col = (row); \
                        state_count m_i__c_row_nodiag; \
                        const int * m_colrow; \
                        double * m_valrow; \
                        if ( NULL == (p_mtx) || (unsigned) m_col >= \
                                                (unsigned) mtx_rows((p_mtx)) ) \
                        { \
                                exit(err_macro_2(err_PARAM, "mtx_change_row_" \
                                        "nodiag(%p,%d,col,p_val)", \
                                        (void*)(p_mtx), m_col, EXIT_FAILURE)); \
                        } \
                        m_colrow = (p_mtx)->valstruc[m_col].col; \
                        m_valrow = (p_mtx)->valstruc[m_col].val; \
                        for ( m_i__c_row_nodiag = mtx_next_num((p_mtx),m_col) ;\
                                        0 < m_i__c_row_nodiag-- ; \
                                                m_colrow++, m_valrow++ ) \
                        { \
                                double * const m_p_val = m_valrow; \
                                m_col = *m_colrow;

#define end_mtx_change_row_nodiag \
                        } \
                        /*@-noeffect@*/(void) m_i__c_row_nodiag;/*@=noeffect@*/\
                }

/*@iter mtx_change_all(sef sparse * p_mtx, yield state_index m_row,
        yield state_index m_col, yield double * m_p_val)@*/
#define mtx_change_all(p_mtx,m_row,m_col,m_p_val) \
                { \
                        state_index m_row, m_col; \
                        double * m_p_val; \
                        state_count m_i__c_all; \
                        struct values * m_values; \
                        const int * m_colrow; \
                        double * m_valrow; \
                        if ( NULL == (p_mtx) ) { \
                                exit(err_macro_1(err_PARAM, "mtx_change_all(" \
                                        "%p,row,col,p_val)", (void *) (p_mtx), \
                                        EXIT_FAILURE)); \
                        } \
                        m_row = mtx_rows((p_mtx)); \
                        m_values = &(p_mtx)->valstruc[m_row]; \
                        (void) (m_col = --m_row, --m_values, \
                        m_i__c_all = m_values->ncols + 1, \
                        m_colrow = m_values->col, m_valrow = m_values->val, \
                        m_p_val = m_col < mtx_cols((p_mtx)) ? &m_values->diag \
                                                            : NULL, \
                        TRUE); \
                        for ( ; 0 < m_i__c_all ; \
                                0 >= --m_i__c_all \
                                     ? (void) (0 < m_row && \
                                       (m_col = --m_row, --m_values, \
                                        m_i__c_all = m_values->ncols + 1, \
                                        m_colrow = m_values->col, \
                                        m_valrow = m_values->val, \
                                        m_p_val = m_col < mtx_cols((p_mtx)) \
                                                ? &m_values->diag : NULL, \
                                        TRUE)) \
                                     : (void) (m_col = *m_colrow++, \
                                        m_p_val = m_valrow++) ) \
                        { \
                                if ( NULL != m_p_val ) {

#define end_mtx_change_all      } \
                        } \
                        /*@-noeffect@*/ (void) m_i__c_all; /*@=noeffect@*/ \
                }


/*======================================================================*/
/************************************************************************/
/************************A trivial MATRIX ALLOCATION*********************/
/*======================================================================*/

/**
* Allocates a sparse matrix
* @param rows the number of rows in the matrix.
* @param cols the number of cols in the matrix.
* @return the pointer to the newly created sparse matrix
*/
/*@only@*/ /*@null@*/
sparse * allocate_sparse_matrix(const int rows, const int cols){
        sparse * pMatrix;

        if ( 0 >= rows || 0 >= cols ) {
                err_msg_2(err_PARAM, "allocate_sparse_matrix(%d,%d)", rows,
                                cols, NULL);
        }
        pMatrix = (sparse *) calloc(1, sizeof(sparse)-sizeof(pMatrix->valstruc)
                                + rows * sizeof(pMatrix->valstruc[0]));
        if ( NULL == pMatrix ) {
                err_msg_2(err_MEMORY, "allocate_sparse_matrix(%d,%d)", rows,
                                cols, NULL);
        }
	pMatrix->rows = rows;
	pMatrix->cols = cols;
	return pMatrix;
}

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
* NOTE: Matrices which were allocated with this method are faster to
*	work with, at the same time they should be treated with care
*	because you can not just add - reallocate rows in them
*	So there is a specific set of methods to be used for this.
*/
/*@only@*/ /*@null@*/
sparse * allocate_sparse_matrix_ncolse(const int rows, const int cols, const int *ncolse){
	int i, sum=0;
	int *col_val;
	double *val;
        sparse * pMatrix;

        if ( 0 >= rows || 0 >= cols || NULL == ncolse ) {
                err_msg_3(err_PARAM, "allocate_sparse_matrix_ncolse(%d,%d,%p)",
                                rows, cols, (const void *) ncolse, NULL);
        }

        /* In an ncolse matrix, I allocate one valstruc.col too much to
           allow for checking the space size even for the last row. */
        pMatrix = (sparse*) calloc(1, sizeof(sparse) - sizeof(pMatrix->valstruc)
                                + sizeof(pMatrix->valstruc[rows].col)
                                + rows * sizeof(pMatrix->valstruc[0]));
        if ( NULL == pMatrix ) {
                err_msg_3(err_MEMORY, "allocate_sparse_matrix_ncolse(%d,%d,%p)",
                                rows, cols, (const void *) ncolse, NULL);
        }
	pMatrix->rows = rows;
	pMatrix->cols = cols;
	for(i=0;i<rows;i++){
		sum+=ncolse[i];
	}

        /* Allocate the continuous arrays for all non-zero column indices and
           values */
        col_val = (int *) calloc((size_t) sum, sizeof(int));
        val = (double *) calloc((size_t) sum, sizeof(double));
        if ( NULL == col_val || NULL == val ) {
                err_msg_3(err_MEMORY, "allocate_sparse_matrix_ncolse(%d,%d,%p)",
                                rows, cols, (const void *) ncolse,
                                (free(val),free(col_val), free(pMatrix), NULL));
        }

	/* Initialize pointers for the matrix rows */
	sum=0;
        for ( i = 0 ; ; i++ )
	{
		pMatrix->valstruc[i].col = &col_val[sum];
                /* valstruc[rows].col is also allocated */
                if ( i >= rows )
                        break;
                /* but valstruc[rows].val or ncolse[rows] may generate a bus
                   error */
	 	pMatrix->valstruc[i].val = &val[sum];
		sum+=ncolse[i];
	}
	return pMatrix;
}

/**
* Frees the sparse matrix allocated with the allocate_sparse_matrix_ncolse(...) method.
* WARNING: This method should be used with the matrices allocated with
*	the allocate_sparse_matrix_ncolse(...) method.
* @param pM the matrix to remove
*/
err_state free_sparse_ncolse(/*@only@*/ /*@i1@*/ /*@null@*/ sparse * pM) {
	int i;
	if( pM != NULL ) {
		if(pM->rows==pM->cols){
                        for ( i = 0 ; i< mtx_rows(pM) ; i++ ) {
                                free(pM->valstruc[i].back_set);
			}
		}
                free(pM->valstruc[0].col);
                free(pM->valstruc[0].val);
		free(pM);
	}
        else {
                err_msg_3(err_PARAM, "free_sparse_ncolse(%p[%dx%d])",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, err_ERROR);
        }
        return err_OK;
}

/**
* Sets a value to a pre-allocated element in the matrix.
* WARNING: This method should be used with the matrices allocated with
*	the allocate_sparse_matrix_ncolse(...) method.
*       It assumes that the matrix be filled from left to right and from top to
*       bottom: It appends new elements to the end of the respective arrays.
* @param pM the matrix to add element to
* @param row the row number in which the value is to be set.
* @param col the column number in which the value is to be set.
* @param val the value to be set.
*/
err_state set_mtx_val_ncolse(/*@i1@*/ /*@null@*/ sparse * pM,
                int row, int col, double val)
{
        int idx, idxp = 0;

        if ( NULL == pM || row < 0 || row >= mtx_rows(pM) || col < 0
                                || col >= mtx_cols(pM) )
        {
                err_msg_6(err_PARAM, "set_mtx_val_ncolse(%p[%dx%d],%d,%d,%g)",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, row, col,
                                val, err_ERROR);
        }

		/* Here pM->ncols[row] is initially 0. */
		if(val != 0 && row != col){
                        idx = mtx_next_num(pM, row);
                        /* test if enough space for another entry has been
                           reserved */
                        if ( &pM->valstruc[row].col[idx]
                                                >= pM->valstruc[row+1].col )
                        {
                                /* Too many elements in this row */
                                err_msg_6(err_INCONSISTENT, "set_mtx_val_ncolse"
                                        "(%p[%dx%d],%d,%d,%g)", (void *) pM,
                                        mtx_rows(pM), mtx_cols(pM), row, col,
                                        val, err_ERROR);
                        }
                        pM->valstruc[row].col[idx] = col; /* Set the new nonzero
                                                        element column id */
                        pM->valstruc[row].val[idx] = val; /* Set the new nonzero
                                                        element value */
                        /* Check if it is a square matrix */
                        if ( mtx_rows(pM) == mtx_cols(pM) ) {
                                int * temp_back_set;
                                /* Get the number of states from which we can
                                   reach state col */
                                idxp = mtx_prev_num(pM, col);
                                /* increase the size of the back_set array */
                                temp_back_set = realloc(
                                                pM->valstruc[col].back_set,
                                                (idxp + 1) * sizeof(int));
                                if ( NULL == temp_back_set ) {
                                        /*@-usereleased@*/ /*@-compdef@*/
                                        err_msg_6(err_MEMORY,
                                                "set_mtx_val_ncolse(%p[%dx%d],"
                                                "%d,%d,%g)", (void *) pM,
                                                mtx_rows(pM), mtx_cols(pM), row,
                                                col, val, err_ERROR);
                                        /*@=usereleased@*/ /*@=compdef@*/
                                }
                                pM->valstruc[col].back_set = temp_back_set;
                                /* Add the current row to this set */
                                temp_back_set[idxp] = row;
                                /* Increase the counter for the number of
                                   elements in the back_set */
				pM->valstruc[col].pcols++;
			}
                        pM->valstruc[row].ncols++; /* Increase the number of
                                        nonzero entries in the current row */
		} else {
			if(val != 0 && row == col){
                                /* Handle the diagonal element */
                                mtx_set_diag_val_nt(pM, row, val);
			}
		}
        return err_OK;
}

/**
* Searches for a pre-allocated matrix element and adds a value to it.
* WARNING: This method should be used with the matrices allocated with
*	the allocate_sparse_matrix_ncolse(...) method.
* @param pM the matrix
* @param row the row number in which the value is to be added.
* @param col the column number in which the value is to be added.
* @param val the value to be added.
*/
err_state add_mtx_val_ncolse(/*@i1@*/ /*@null@*/ sparse * pM,
                int row, int col, double val)
{
        int ncols, i;
	BOOL found = FALSE;
        if ( NULL == pM || 0 > row || row >= mtx_rows(pM) || 0 > col
                        || col >= mtx_cols(pM) )
                err_msg_6(err_PARAM, "add_mtx_val_ncolse(%p[%dx%d],%d,%d,%g)",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, row, col, val,
                                err_ERROR);
        if ( /*@-realcompare@*/ 0.0 != val /*@=realcompare@*/ ) {
		if ( row != col ) {
                        ncols = mtx_next_num(pM, row);
			/* TODO: We could be smart and make at least a binary */
			/* search here, as soon as all elements in pM->val[row].col */
                        /* array are actually required to be ordered by value */
			/* (increasing with the arrayindex) */
			for( i=0; i < ncols; i++ ) {
                                if ( pM->valstruc[row].col[i] == col ) {
                                        pM->valstruc[row].val[i] += val;
					found = TRUE;
					break;
				}
			}
		} else {
                        mtx_set_diag_val_nt(pM, row,
                                        mtx_get_diag_val_nt(pM, row) + val);
			found = TRUE;
		}
		if( found == FALSE ) {
                        if ( err_state_iserror(set_mtx_val_ncolse(pM, row, col,
                                                                val)) )
                                err_msg_6(err_CALLBY, "add_mtx_val_ncolse"
                                        "(%p[%dx%d],%d,%d,%g)", (void *) pM,
                                        mtx_rows(pM), mtx_cols(pM), row, col,
                                        val, err_ERROR);
		}
        }
        return err_OK;
}

/*======================================================================*/
/************************************************************************/
/************************General Sparse matrix methods*******************/
/**********For matrices allocated with "allocate_sparse_matrix"**********/
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
*/
err_state initMatrix(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pStateSpace,
                /*@i1@*/ /*@null@*/ sparse * pQ,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pValidStates)
{
	int i;
	int valid_state;

        if ( NULL == pStateSpace || NULL == pQ
                        || NULL == pValidStates || 0 > pValidStates[0] )
        {
                err_msg_8(err_PARAM, "initMatrix(%p[%dx%d],%p[%dx%d],%p[%d])",
                                (const void *) pStateSpace,
                                NULL != pStateSpace ? mtx_rows(pStateSpace) : 0,
                                NULL != pStateSpace ? mtx_cols(pStateSpace) : 0,
                                (void *) pQ, NULL != pQ ? mtx_rows(pQ) : 0,
                                NULL != pQ ? mtx_cols(pQ) : 0,
                                (const void *) pValidStates,
                                NULL != pValidStates ? pValidStates[0] : 0,
                                err_ERROR);
        }

	/*Form the pQ matrix*/
	for( i = 1; i <= pValidStates[0] ; i++ )
	{
		valid_state = pValidStates[i];
                if ( err_state_iserror(mtx_copy_row(pQ, valid_state,
                                        pStateSpace)) )
                {
                        err_msg_8(err_CALLBY, "initMatrix(%p[%dx%d],%p[%dx%d],"
                                "%p[%d])", (const void *) pStateSpace,
                                mtx_rows(pStateSpace), mtx_cols(pStateSpace),
                                (void *) pQ, mtx_rows(pQ), mtx_cols(pQ),
                                (const void *) pValidStates, pValidStates[0],
                                err_ERROR);
                }
	}
        return err_OK;
}

/**
* This method removes all the non zero rows from the Q matrix of the BSCC
* @param pQ the BSCC nodes of the original matrix
* @param pValidStates the valid states i.e. ids of the rows which have to be cleaned
*                     this array contains the number of nodes as the first element
*/
err_state cleanMatrix(/*@i1@*/ /*@null@*/ sparse * pQ,
                /*@i1@*/ /*@null@*/ const int * pValidStates)
{
	/*Remove the non zero rows*/
        int i, length;

        if ( NULL == pQ || NULL == pValidStates || 0 > pValidStates[0] ) {
                err_msg_5(err_PARAM, "cleanMatrix(%p[%dx%d],%p[%d])",
                                (void *) pQ, NULL != pQ ? mtx_rows(pQ) : 0,
                                NULL != pQ ? mtx_cols(pQ) : 0,
                                (const void *) pValidStates,
                                NULL != pValidStates ? pValidStates[0] : 0,
                                err_ERROR);
        }

        length = pValidStates[0];
	for( i = 1; i <= length ; i++ )
	{
                struct values * valrow = &pQ->valstruc[*++pValidStates];
                valrow->diag = 0.0;
                if ( 0 != valrow->ncols ) {
                        valrow->ncols = 0;
                        valrow->col = (int *) NULL;
                        valrow->val = (double *) NULL;
                        free(valrow->back_set);
                        valrow->back_set = (int *) NULL;
                        valrow->pcols = 0;
                }
	}
        /*Note: Because of the sparse matrix functions we do not need to
	restore the pStateSpace->val diagonal values as we do not modify them*/
        return err_OK;
}

/*****************************************************************************
name            : free_mtx_sparse
role		: frees the sparse matrix.
@param		: sparse * pM the matrix to remove
remark		:
******************************************************************************/
err_state free_mtx_sparse(/*@only@*/ /*@i1@*/ /*@null@*/ sparse * pM)
{
        state_index i, rows;
	if(pM)
	{
                rows = mtx_rows(pM);
		for(i=0;i<rows;i++)
		{
                        free(pM->valstruc[i].back_set);
                        free(pM->valstruc[i].col);
                        free(pM->valstruc[i].val);
		}
		free(pM);
	}
        else {
                err_msg_3(err_PARAM, "free_mtx_sparse(%p[%dx%d])",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, err_ERROR);
        }
        return err_OK;
}

/**
* ADD a NEW value to the matrix.
* WARNING: Should NOT be used with the matrices allocated with
* allocate_sparse_matrix_ncolse(,,) method!
* @param	: sparse * pM the matrix to add element to
* @param	: int row: The row number in which the value is to be set.
* @param	: int col: The column number in which the value is to be set.
* @param	: double val: The value to be set.
* @return	: int: 0: fail, 1: success
* NOTE: return 0 indicates that row/col does not meet the bounds of matrix
*/
err_state set_mtx_val(/*@i1@*/ /*@null@*/ sparse * pM,
                int row, int col, double val)
{
        /* Correction of an MRMC bug: It is an error if row == mtx_rows(pM) or
           col == mtx_cols(pM). David N. Jansen. */
        if( NULL == pM || row < 0 || row >= mtx_rows(pM) || col < 0
                                || col >= mtx_cols(pM) )
        {
                err_msg_6(err_PARAM, "set_mtx_val(%p[%dx%d],%d,%d,%g)",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, row, col, val,
                                err_ERROR);
        }

	if( val != 0 && row != col )
	{
                int idx = mtx_next_num(pM, row);
                /*@null@*/ int *temp_col = (int*) realloc(pM->valstruc[row].col,
                                (idx + 1) * sizeof(int));
                /*@null@*/ double * temp_val = NULL;
                if ( NULL != temp_col ) {
                        pM->valstruc[row].col = temp_col;
                        temp_val = (double *) realloc(pM->valstruc[row].val,
                                (idx + 1) * sizeof(double));
                        if ( NULL != temp_val ) {
                                pM->valstruc[row].val = temp_val;
                                if ( mtx_rows(pM) == mtx_cols(pM) ) {
                                        int idxp = mtx_prev_num(pM, col);
                                        /*@null@*/ int * temp_back_set =
                                                (int *) realloc(
                                                pM->valstruc[col].back_set,
                                                (idxp + 1) * sizeof(int));
                                        if ( NULL != temp_back_set ) {
                                                pM->valstruc[col].back_set
                                                        = temp_back_set;
                                                temp_back_set[idxp] = row;
                                                pM->valstruc[col].pcols++;
                                        } else
                                                temp_col = NULL;
                                }
                        }
                }
                if ( NULL == temp_col || NULL == temp_val ) {
                        err_msg_6(err_MEMORY, "set_mtx_val(%p[%dx%d],%d,%d,%g)",
                                        (void *) pM, mtx_rows(pM), mtx_cols(pM),
                                        row, col, val, err_ERROR);
                }
                temp_col[idx] = col;
                temp_val[idx] = val;
                pM->valstruc[row].ncols++;
	}
	else if( val != 0 && row == col )
                mtx_set_diag_val_nt(pM, row, val);
        return err_OK;
}

/*****************************************************************************
name		: free_mtx_wrap
role		: frees the sparse matrix. Does not clean rows.
@param		: sparse * pM the matrix to remove
******************************************************************************/
err_state free_mtx_wrap(/*@only@*/ /*@i1@*/ /*@null@*/ sparse * pM)
{
	if(pM)
	{
                values * valrow = &pM->valstruc[0];
                state_index i = mtx_rows(pM);
                do {
                        free(valrow++->back_set);
                } while ( 0 < --i );
		free(pM);
	}
        else {
                err_msg_3(err_PARAM, "free_mtx_wrap(%p[%dx%d])",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, err_ERROR);
        }
        return err_OK;
}

/**
* Add a value to an existing one in the matrix or inserts a new value if the is none.
* WARNING: Should NOT be used with the matrices allocated with
* allocate_sparse_matrix_ncolse(,,) method! Unless you are 100% sure that
* the pM[row,col] value already exists in the matrix!
* @param	: sparse * pM the matrix
* @param        : int row: The row number in which the value is to be added.
* @param	: int col: The column number in which the value is to be added.
* @param	: double val: to be added.
*/
err_state add_mtx_val(/*@i1@*/ /*@null@*/ sparse * pM,
                int row, int col, double val)
{
        int ncols, i;
        BOOL found = FALSE;
        /* Correction of an MRMC bug: It is an error if row == mtx_rows(pM) or
           col == mtx_cols(pM). David N. Jansen. */
        if ( NULL == pM || row < 0 || row >= mtx_rows(pM) || col < 0
                                || col >= mtx_cols(pM) )
        {
                err_msg_6(err_PARAM, "add_mtx_val(%p[%dx%d],%d,%d,%g)",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, row, col, val,
                                err_ERROR);
        }

        if ( /*@-realcompare@*/ 0.0 == val /*@=realcompare@*/ ) {
                return err_OK;
        }
	if ( row != col )
	{
                ncols = mtx_next_num(pM, row);
		for( i=0; i < ncols; i++ )
		{
                        if ( pM->valstruc[row].col[i] == col )
			{
                                pM->valstruc[row].val[i] += val;
				found=1;
				break;
			}
		}
	}
	else
	{
                mtx_set_diag_val_nt(pM, row, mtx_get_diag_val_nt(pM,row) + val);
		found = 1;
	}
	if( found == 0 )
                if ( err_state_iserror(set_mtx_val(pM, row, col, val)) ) {
                        err_msg_6(err_CALLBY, "add_mtx_val(%p[%dx%d],%d,%d,%g)",
                                (void *) pM, mtx_rows(pM), mtx_cols(pM), row,
                                col, val, err_ERROR);
                }
        return err_OK;
}

/*****************************************************************************
name		: get_mtx_val
role		: get a value in the matrix. Perform binary search to get a
		non-diagonal value in the matrix when the number of
		non-diagonal elements in a row is >4.
@param		: sparse * pM the matrix to get element from
@param		: int row: The row number in which the value is to be set.
@param		: int col: The column number in which the value is to be set.
@param		: double * val (by reference): The value.
@return         : err_ERROR: fail, err_OK: success
******************************************************************************/
err_state get_mtx_val(/*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM,
                int row, int col, /*@out@*/ /*@i1@*/ /*@null@*/ double * val)
{
        int i;

        /* Correction of an MRMC bug: It is an error if row == mtx_rows(pM) or
           col == mtx_cols(pM). David N. Jansen. */
        if ( NULL == pM || row < 0 || row >= mtx_rows(pM) || col < 0
                                || col >= mtx_cols(pM) || NULL == val )
	{
                /*@-mustdefine@*/
                err_msg_6(err_PARAM, "get_mtx_val(%p[%dx%d],%d,%d,%p)",
                                (const void *)pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, row, col,
                                (void *) val, err_ERROR);
                /*@=mustdefine@*/
	}

	if ( row != col )
	{
                int high = mtx_next_num(pM, row) - 1, low = 0;
                                /* high and low are used for binary search */
                *val = 0.0;
		/* If there are <= than 4 elements in a row then we can do a regular search */
		if(high > 4){
			while(low < high){
				i = (low+high)/2;
                                if ( pM->valstruc[row].col[i] == col ) {
                                        *val = pM->valstruc[row].val[i];
					break;
                                } else if ( col > pM->valstruc[row].col[i] ) {
					low = i + 1;
				}else{
					high = i - 1;
				}
			}
			/* low == high => there may be no required element in the (row,col) cell */
                        if ( pM->valstruc[row].col[low] == col )
			{
                                *val = pM->valstruc[row].val[low];
			}
		}else{
                        const int ncols = mtx_next_num(pM, row);
			for( i=0; i < ncols; i++ )
			{
                                if ( pM->valstruc[row].col[i] == col )
				{
                                        *val = pM->valstruc[row].val[i];
					break;
				}
			}
		}
	}
        else
	{
                *val = mtx_get_diag_val_nt(pM, row);
	}
        return err_OK;
}

/*****************************************************************************
name		: get_mtx_row_sums
role		: add the elements in each row, including the diagonal element.
@param		: sparse * pM the matrix to work with
@return		: double *: row sums for each row.
remark		:
******************************************************************************/
/*@only@*/ /*@null@*/ double * get_mtx_row_sums(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM)
{
        int size;
        double * row_sum;

        if ( NULL == pM ) {
                err_msg_3(err_PARAM, "get_mtx_row_sums(%p[%dx%d])",
                                (const void *)pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, NULL);
        }

        size = pM->rows;
        row_sum = (double *) calloc((size_t) size, sizeof(double));
        if ( NULL == row_sum ) {
                err_msg_3(err_MEMORY, "get_mtx_row_sums(%p[%dx%d])",
                                (const void*)pM,mtx_rows(pM),mtx_cols(pM),NULL);
        }
        mtx_walk_all(pM, row, col, val)
	{
                row_sum[row] += val;
	}
        end_mtx_walk_all;
	return row_sum;
}

/*****************************************************************************
name		: get_mtx_cer_row_sums
role		: add elements in certain rows and return the vector of values.
@param		: sparse * pM the matrix to work with
@param		: int * pIds the ids of the rows for which the sums should be computed
			the pIds[0] contains the number of row ids
@return		: the rows' sums
******************************************************************************/
/*@only@*/ /*@null@*/ double * get_mtx_cer_row_sums(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pIds)
{
        int i, num;
        double * pE, * pE_i;

        if ( NULL == pM || NULL == pIds || 0 > pIds[0] ) {
                err_msg_5(err_PARAM, "get_mtx_cer_row_sums(%p[%dx%d],%p[%d])",
                                (const void *)pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) pIds, NULL != pIds ? pIds[0] : 0,
                                NULL);
        }

        num = pIds[0];
        pE = (double *) calloc((size_t) num, sizeof(double));
        if ( NULL == pE ) {
                err_msg_5(err_MEMORY, "get_mtx_cer_row_sums(%p[%dx%d],%p[%d])",
                                (const void *) pM, mtx_rows(pM), mtx_cols(pM),
                                (const void *) pIds, pIds[0], NULL);
        }
        pE_i = pE;
	for( i=1; i<= num; i++ )
	{
		double result=0.0;
                const state_index id = *++pIds;
                mtx_walk_row(pM, id, col, val) {
                        result += val;
                } end_mtx_walk_row;
                *pE_i++ = result;
	}
	return pE;
}

/*****************************************************************************
name		: mult_mtx_const
role		: multiply all elements in the matrix with a constant.
@param		: sparse * pM the matrix to multiply
@param		: double constant: The constant value.
******************************************************************************/
err_state mult_mtx_const(/*@i1@*/ /*@null@*/ sparse * pM,
                double constant)
{
        if ( NULL == pM ) {
                err_msg_4(err_PARAM, "mult_mtx_const(%p[%dx%d],%g)", (void *)pM,
                                NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, constant,
                                err_ERROR);
        }

        mtx_change_all(pM, row, col, pVal)
	{
                *pVal *= constant;
	}
        end_mtx_change_all;
        return err_OK;
}

/*****************************************************************************
name		: mult_mtx_cer_const
role		: multiply all elements in the matrix with a constant.
@param		: sparse * pM the matrix to multiply
@param		: double constant: The constant value.
@param		: int * pValidStates the valid rows to multiply pValidStates[0]
		  is the number of states (rows) in the array
******************************************************************************/
err_state mult_mtx_cer_const(/*@i1@*/ /*@null@*/ sparse * pM,
                double constant,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pValidStates)
{
        int i, id, length;

        if ( NULL == pM || NULL == pValidStates || 0 > pValidStates[0] ) {
                err_msg_6(err_PARAM, "mult_mtx_cer_const(%p[%dx%d],%g,%p[%d])",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, constant,
                                (const void *) pValidStates,
                                NULL != pValidStates ? pValidStates[0] : 0,
                                err_ERROR);
        }

        length = pValidStates[0];
	for( i=1; i <= length; i++ )
	{
		id = pValidStates[i];
                mtx_change_row(pM, id, col, pVal) {
                        *pVal *= constant;
                } end_mtx_change_row;
	}
        return err_OK;
}

/**
* Multiply all elements in all i'th rows of the matrix by a constant[i] (or 1/constant[i] ).
* @param pM the matrix to multiply
* @param constant The constant value.
* @param pValidStates the valid rows to multiply pValidStates[0]
*			is the number of states (rows) in the array
* @param invert if == 0 then nothing, otherwise we invert the elements of constant
*			before the multiplication.
*/
err_state mult_mtx_cer_const_array(/*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * constant,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pValidStates,
                BOOL invert)
{
        int i, id, length;
	double cons;

        if ( NULL == pM || NULL == constant
                                || NULL == pValidStates || 0 > pValidStates[0] )
        {
                err_msg_7(err_PARAM,
                        "mult_mtx_cer_const_array(%p[%dx%d],%p,%p[%d],%d)",
                        (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                        NULL != pM ? mtx_cols(pM) : 0, (const void *) constant,
                        (const void *) pValidStates,
                        NULL != pValidStates ? pValidStates[0] : 0, (int)invert,
                        err_ERROR);
        }

        length = pValidStates[0];
	for( i=1; i <= length; i++ ){
		id = pValidStates[i];
		cons = constant[id];

		if ( invert && cons ) cons = 1 / cons;
                mtx_change_row(pM, (const state_index) id, col, pVal) {
                        *pVal *= cons;
		}
                end_mtx_change_row;
	}
        return err_OK;
}

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
*	NOTE: the array dimension corresponds to the state-space dimensions
*/
/*@only@*/ /*@null@*/ double * make_embedded_dtmc_out_of_rate_mtx_vs(
                /*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pRowSums,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pValidStates)
{
        int row_idx, length;
        double * pExitRates;

        if ( NULL == pM || NULL == pRowSums
                                || NULL == pValidStates || 0 > pValidStates[0] )
        {
                err_msg_6(err_PARAM, "make_embedded_dtmc_out_of_rate_mtx_vs"
                                "(%p[%dx%d],%p,%p[%d])", (void *) pM,
                                NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) pRowSums,
                                (const void *) pValidStates,
                                NULL != pValidStates ? pValidStates[0] : 0,
                                NULL);
        }

        length = pValidStates[0];
        pExitRates = (double *) calloc((size_t) mtx_rows(pM),
                        sizeof(double));
        if ( NULL == pExitRates ) {
                err_msg_6(err_MEMORY, "make_embedded_dtmc_out_of_rate_mtx_vs"
                                "(%p[%dx%d],%p,%p[%d])", (void *) pM,
                                mtx_rows(pM), mtx_cols(pM),
                                (const void *) pRowSums,
                                (const void *) pValidStates, pValidStates[0],
                                NULL);
        }
	/* Iterate through the allowed states */
	for( row_idx = 1; row_idx <= length; row_idx++ ){
		/* Obtain the allowed state index */
		int id = pValidStates[row_idx];
		/* Get the number of non diagonal elements in the row */
                int non_zero_off_diag_columns = mtx_next_num(pM, id);

		/* if there are non-diagonal elements */
		if( non_zero_off_diag_columns > 0 ){
			/* Create the normalization factor, the diagonal element should be excluded */
			/* because we do not want to have loops in the resulting DTMC */
			/* The normalization factor is basically the exit rate of the state*/
                        double normalization_factor = pExitRates[id]
                                = pRowSums[id] - mtx_get_diag_val_nt(pM, id);

			/* NOTE: Since non_zero_off_diag_columns is positive */
			/* We shall assume that normalization_factor > 0, */
			/* but let's do a safety check just in case */
                        if ( normalization_factor <= 0.0 ) {
                                err_msg_6(err_INCONSISTENT,
                                        "make_embedded_dtmc_out_of_rate_mtx_vs"
                                        "(%p[%dx%d],%p,%p[%d])", (void *) pM,
                                        mtx_rows(pM), mtx_cols(pM),
                                        (const void *) pRowSums,
                                        (const void *) pValidStates,
                                        pValidStates[0],
                                        (free(pExitRates), NULL));
                        }
				/* Get the non-diagonal row elements */
                        mtx_change_row_nodiag(pM, id, col, pVal) {
					/* Normalize the off-diagonal elements of the row */
                                *pVal /= normalization_factor;
				}
                        end_mtx_change_row_nodiag;
		}
	}

	return pExitRates;
}

/**
* This method restores the rate matrix from the "embedded" DTMC obtained
* with the "make_embedded_dtmc_out_of_rate_mtx" method.
* WARNING: The diagonal elements remain untouched, and therefore
* we expect them to be preserved from before the matrix was turned in to an embedded DTMC.
* @param pM the matrix to operate on
* @param pValidStates the valid rows of the matrix to operate on.
*		NOTE: The element pValidStates[0] should store
*		the number of states (rows) in the array.
* @param pExitRates the array of exit rates returned by make_embedded_dtmc_out_of_rate_mtx
* WARNING: The pExitRates has to be obtained with the same pValidStates on the same pM matrix
*		using the make_embedded_dtmc_out_of_rate_mtx method
*/
err_state restore_rate_mtx_out_of_embedded_dtmc_vs(
                /*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pValidStates,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pExitRates)
{
        int row_idx, length;

        if ( NULL == pM || NULL == pValidStates
                                || 0 > pValidStates[0] || NULL == pExitRates )
        {
                err_msg_6(err_PARAM, "restore_rate_mtx_out_of_embedded_dtmc_vs"
                                "(%p[%dx%d],%p[%d],%p)", (void *) pM,
                                NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) pValidStates,
                                NULL != pValidStates ? pValidStates[0] : 0,
                                (const void *) pExitRates, err_ERROR);
        }

        length = pValidStates[0];
		/* Iterate through the allowed states */
		for( row_idx = 1; row_idx <= length; row_idx++ ){
			/* Obtain the allowed state index */
			int id = pValidStates[row_idx];
			/* Get the number of non diagonal elements in the row */
                        int non_zero_off_diag_columns = mtx_next_num(pM, id);

			/* if there are non-diagonal elements */
			if( non_zero_off_diag_columns > 0 ){
				/* Take the pre-computed exit rate */
				double normalization_factor =  pExitRates[id];
				/* WARNING: We assume that pMtxRow.diag contains the diagonal */
				/* element of the rate mtx, i.e. it was not modified */

				/* Get the non-diagonal row elements */
                                mtx_change_row_nodiag(pM, id, col, pVal) {
					/* Restore the off-diagonal elements of the row */
                                        *pVal *= normalization_factor;
				}
                                end_mtx_change_row_nodiag;
			}
		}
        return err_OK;
}

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
*	NOTE: the array dimension corresponds to the state-space dimensions
*/
/*@only@*/ /*@null@*/ double * make_embedded_dtmc_out_of_rate_mtx_all(
                /*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double *pRowSums)
{
        int row_idx, length;
        double * pExitRates;

        if ( NULL == pM || NULL == pRowSums ) {
                err_msg_4(err_PARAM, "make_embedded_dtmc_out_of_rate_mtx_all"
                                "(%p[%dx%d],%p)", (void *) pM,
                                NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) pRowSums, NULL);
        }

        pExitRates = (double *) calloc((size_t) mtx_rows(pM), sizeof(double));
        if ( NULL == pExitRates ) {
                err_msg_4(err_MEMORY, "make_embedded_dtmc_out_of_rate_mtx_all"
                                "(%p[%dx%d],%p)", (void *) pM, mtx_rows(pM),
                                mtx_cols(pM), (const void *) pRowSums, NULL);
        }

        length = mtx_rows(pM);
	/* Iterate through the allowed states */
	for( row_idx = 0; row_idx < length; row_idx++ ){
		/* Get the number of non diagonal elements in the row */
                int non_zero_off_diag_columns = mtx_next_num(pM, row_idx);

		/* if there are non-diagonal elements */
		if( non_zero_off_diag_columns > 0 ){
			/* Create the normalization factor, the diagonal element should be excluded */
			/* because we do not want to have loops in the resulting DTMC */
			/* The normalization factor is basically the exit rate of the state*/
                        double normalization_factor = pExitRates[row_idx]
                                        = pRowSums[row_idx] -
                                                mtx_get_diag_val_nt(pM,row_idx);

			/* NOTE: Since non_zero_off_diag_columns is positive */
			/* We shall assume that normalization_factor > 0, */
			/* but let's do a safety check just in case */
                        if ( normalization_factor <= 0.0 ) {
                                err_msg_4(err_INCONSISTENT,
                                        "make_embedded_dtmc_out_of_rate_mtx_all"
                                        "(%p[%dx%d],%p)", (void *) pM,
                                        mtx_rows(pM), mtx_cols(pM),
                                        (const void *) pRowSums,
                                        (free(pExitRates), NULL));
                        }
				/* Get the non-diagonal row elements */
                        mtx_change_row_nodiag(pM, (const state_index) row_idx,
                                                        column_idx, pVal)
                        {
					/* Normalize the off-diagonal elements of the row */
                                *pVal /= normalization_factor;
				}
                        end_mtx_change_row_nodiag;
		}
	}

	return pExitRates;
}

/**
* This method restores the rate matrix from the "embedded" DTMC obtained
* with the "make_embedded_dtmc_out_of_rate_mtx_all" method.
* WARNING: The diagonal elements remain untouched, and therefore
* we expect them to be preserved from before the matrix was turned in to an embedded DTMC.
* @param pM the matrix to operate on
* @param pExitRates the array of exit rates returned by make_embedded_dtmc_out_of_rate_mtx
*/
err_state restore_rate_mtx_out_of_embedded_dtmc_all(
                /*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pExitRates)
{
        int row_idx, length;

        if ( NULL == pM || NULL == pExitRates ) {
                err_msg_4(err_PARAM, "restore_rate_mtx_out_of_embedded_dtmc_all"
                                "(%p[%dx%d],%p)", (void *) pM,
                                NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) pExitRates, err_ERROR);
        }

        length = mtx_rows(pM);
		/* Iterate through the allowed states */
		for( row_idx = 0; row_idx < length; row_idx++ ){
			/* Get the number of non diagonal elements in the row */
                        int non_zero_off_diag_columns=mtx_next_num(pM, row_idx);

			/* if there are non-diagonal elements */
			if( non_zero_off_diag_columns > 0 ){
				/* Take the pre-computed exit rate */
				double normalization_factor =  pExitRates[row_idx];
				/* WARNING: We assume that pMtxRow.diag contains the diagonal */
				/* element of the rate mtx, i.e. it was not modified */

				/* Get the non-diagonal row elements */
                                mtx_change_row_nodiag(pM,
                                        (const state_index)row_idx, dummy, pVal)
                                {
					/* Restore the off-diagonal elements of the row */
                                        *pVal *= normalization_factor;
				}
                                end_mtx_change_row_nodiag;
			}
		}

        return err_OK;
}

/*****************************************************************************
name            : add_mtx_diagonal
role		: add the elements of the input array to the diagonal one-one.
@param		: sparse * pM the matrix to work with
@param		: double *diagonal: The addendum to the diagonal.
@return         : err_ERROR: fail, err_OK: success
remark          : return err_ERROR: sparse matrix is not valid.
		  size should be correct.
******************************************************************************/
err_state add_mtx_diagonal(/*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * diagonal)
{
	int i, size;

        if ( NULL == pM || NULL == diagonal ) {
                err_msg_4(err_PARAM, "add_mtx_diagonal(%p[%dx%d],%p)",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) diagonal, err_ERROR);
        }

	size = pM->rows;
	for( i=0; i < size; i++ )
                mtx_set_diag_val_nt(pM, i,
                                mtx_get_diag_val_nt(pM, i) + diagonal[i]);
        return err_OK;
}

/*****************************************************************************
name		: sub_mtx_diagonal
role            : sub the elements of the input array from the diagonal one-one.
@param		: sparse * pM the matrix to work with
@param          : double *diagonal: The subtrahend from the diagonal.
@return         : err_ERROR: fail, err_OK: success
remark          : return err_ERROR: sparse matrix is not valid.
		  size should be correct.
******************************************************************************/
err_state sub_mtx_diagonal(/*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * diagonal)
{
	int i, size;

        if ( NULL == pM || NULL == diagonal ) {
                err_msg_4(err_PARAM, "sub_mtx_diagonal(%p[%dx%d],%p)",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) diagonal, err_ERROR);
        }

	size = pM->rows;
	for( i=0; i < pM->rows; i++ )
	{
                mtx_set_diag_val_nt(pM, i,
                                mtx_get_diag_val_nt(pM, i) - diagonal[i]);
	}
        return err_OK;
}

/****************************************************************************
name            : add_cer_cons_diagonal_mtx
role            : add the elements of the input array to the diagonal one-one.
@param		: sparse * pM the matrix to work with
@param		: int * pIds the ids of the rows in which diagonal elements
                  should be increased by the corresponding pE values
		  the pIds[0] contains the number of row ids
@param		: int constant: an operand (cons)
******************************************************************************/
err_state add_cer_cons_diagonal_mtx(/*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pIds,
                int constant)
{
	int i;
        int num;

        if ( NULL == pM || NULL == pIds || 0 > pIds[0] ) {
                err_msg_6(err_PARAM, "add_cer_cons_diagonal_mtx"
                                "(%p[%dx%d],%p[%d],%d)", (void *) pM,
                                NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) pIds, NULL != pIds ? pIds[0] : 0,
                                constant, err_ERROR);
        }

        num = pIds[0];
	for( i=1; i <= num; i++ )
        {
                const state_index id = pIds[i];
                mtx_set_diag_val_nt(pM, id,
                                        mtx_get_diag_val_nt(pM, id) + constant);
        }
        return err_OK;
}

/****************************************************************************
name		: sub_cer_cons_diagonal_mtx
role		: diag elements = cons - diag elements
@param		: sparse * pM the matrix to work with
@param		: int * pIds the ids of the rows in which diagonal elements
                  should be subtracted from the corresponding pE values
		  the pIds[0] contains the number of row ids
@param		: int constant: an operand (cons)
******************************************************************************/
/* The function is never used. Therefore, I have commented it out. David N.
   Jansen.
err_state sub_cer_cons_diagonal_mtx(/ *@i1@* / / *@null@* / sparse * pM,
                / *@observer@* / / *@i1@* / / *@null@* / const int * pIds,
                int constant)
{
	int i;
	int num = pIds[0];
	for( i=1; i <= num; i++ )
	{
                mtx_set_diag_val_nt(pM, pIds[i],
                                constant - mtx_get_diag_val_nt(pM, pIds[i]));
	}
        return err_OK;
}
*/

/****************************************************************************
name		: sub_mtx_cer_diagonal
role		: Subtract elements of the pE array from the predefined diagonal
		  elements of the matrix
@param		: sparse * pM the matrix to work with
@param		: int * pIds the ids of the rows in which diagonal elements
                  should be reduced by the corresponding pE values
		  the pIds[0] contains the number of row ids
@param		: double * pE: The elements to be subtracted
******************************************************************************/
err_state sub_mtx_cer_diagonal(/*@i1@*/ /*@null@*/ sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pIds,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pE)
{
	int i;
        int num;

        if ( NULL == pM || NULL == pIds || 0 > pIds[0] || NULL == pE ) {
                err_msg_6(err_PARAM, "sub_mtx_cer_diagonal"
                                "(%p[%dx%d],%p[%d],%p)", (void *) pM,
                                NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                (const void *) pIds, NULL != pIds ? pIds[0] : 0,
                                (const void *) pE, err_ERROR);
        }

        num = pIds[0];
	for( i=1; i <= num; i++ )
	{
                const state_index id = pIds[i];

                mtx_set_diag_val_nt(pM, id,
                                mtx_get_diag_val_nt(pM, id) - pE[i-1]);
	}
        return err_OK;
}

/*****************************************************************************
name		: add_mtx_cons_diagonal
role            : add a constant to all elements of the diagonal.
@param		: sparse * pM the matrix to work with
@param		: double constant: The constant value.
@return         : err_ERROR: fail, err_OK: success
remark          : return err_ERROR: sparse matrix is not valid.
******************************************************************************/
err_state add_mtx_cons_diagonal(/*@i1@*/ /*@null@*/ sparse * pM,
                double cons)
{
	int i, size;

        if ( NULL == pM ) {
                err_msg_4(err_PARAM, "add_mtx_cons_diagonal(%p[%dx%d],%g)",
                                (void *) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, cons, err_ERROR);
        }

	size = pM->rows;
	for(i=0;i<pM->rows;i++)
                mtx_set_diag_val_nt(pM, i, mtx_get_diag_val_nt(pM, i) + cons);
        return err_OK;
}

/*****************************************************************************
name            : mtx_copy_row
role            : copy a row from one matrix to another.
                  NOTE: This function copies the values structure. That is why
                  changes of fields in the original row later on will not affect
                  the one which was set by this function. In particular,
                  changing the diagonal element
		  in the original matrix will not affect the row added to the other matrix
		  by this method.
@param          : sparse * to: the goal matrix to copy to
@param		: int row: The row number in which the value is to be set.
@param          : sparse * from: the source matrix to copy from
@return         : err_ERROR: fail, err_OK: success
remark          : return err_ERROR: indicates that row/col does not meet the
                  bounds
		  of matrix
******************************************************************************/
err_state mtx_copy_row(/*@i1@*/ /*@null@*/ sparse * to,
                int row, /*@i1@*/ /*@null@*/ const sparse * from)
{
        /* Correction of an MRMC bug: It is an error if row == mtx_rows(pM).
           David N. Jansen. */
        if ( NULL == to || (unsigned) row >= (unsigned) mtx_rows(to)
                                || NULL == from
                                || mtx_rows(to) != mtx_rows(from)
                                || mtx_cols(to) != mtx_cols(from)
                                || NULL != to->valstruc[row].col
                                || NULL != to->valstruc[row].val
                                || 0 != mtx_next_num(to, row) )
        {
                err_msg_7(err_PARAM,
                                "mtx_copy_row(%p[%dx%d],%d,%p[%dx%d])",
                                (void *) to, NULL != to ? mtx_rows(to) : 0,
                                NULL != to ? mtx_cols(to) : 0, row,
                                (const void *) from,
                                NULL != from ? mtx_rows(from) : 0,
                                NULL != from ? mtx_cols(from) : 0, err_ERROR);
        }

        free(to->valstruc[row].back_set);
        to->valstruc[row].back_set = NULL;
        to->valstruc[row].pcols = 0;
        mtx_set_diag_val_nt(to, row, mtx_get_diag_val_nt(from, row));
        if ( (to->valstruc[row].ncols = mtx_next_num(from, row)) != 0 ) {
                to->valstruc[row].col = from->valstruc[row].col;
                to->valstruc[row].val = from->valstruc[row].val;
        }
        /*@-compmempass@*/ /* pM is a wrapper matrix, but splint thinks it is
                              a normal one. */
        return err_OK; /*@=compmempass@*/
}

/*****************************************************************************
name		: multiply_mtx_MV
role		: multiply a matrix with a vector.
@param		: sparse * pM: operand matrix.
@param		: double *vec: The operand vector.
@param		: double *res: The resulting vector.
@return         : err_ERROR: fail, err_OK: success
remark		: size should be correct.
******************************************************************************/
err_state multiply_mtx_MV(/*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * vec,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * res)
{
        if ( NULL == pM || NULL == vec || NULL == res ) {
                /*@-mustdefine@*/
                err_msg_5(err_PARAM, "multiply_mtx_MV(%p[%dx%d],%p,%p)",
                                (const void *)pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, (const void*)vec,
                                (void *) res, err_ERROR);
                /*@=mustdefine@*/
        }

        memset(res, (int) '\0', sizeof(double) * mtx_rows(pM));
        mtx_walk_all(pM, row, col, val)
	{
                res[row] += vec[col] * val;
	}
        end_mtx_walk_all;
        return err_OK;
}

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
err_state multiply_mtx_cer_MV(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * vec,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * res,
                int num,/*@observer@*/ /*@i1@*/ /*@null@*/ const int*valid_rows)
{
        int i, v;
        double result;

        if ( NULL == pM || NULL == vec || NULL == res || 0 > num
                                || NULL == valid_rows )
        {
                /*@-mustdefine@*/
                err_msg_7(err_PARAM,
                                "multiply_mtx_cer_MV(%p[%dx%d],%p,%p,%d,%p)",
                                (const void *)pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, (const void*)vec,
                                (void *) res, num, (const void *) valid_rows,
                                err_ERROR);
                /*@=mustdefine@*/
        }

	for(i=0;i<num;i++, valid_rows++)
	{
		v = *valid_rows;
                result = 0.0;

                mtx_walk_row(pM, v, col, val)
		{
                        result += vec[col] * val;
		}
                end_mtx_walk_row;
		res[v]=result;
	}
        /*@-mustdefine@*/ /* If num == 0, nothing should be written */
        return err_OK; /*@=mustdefine@*/
}

/*****************************************************************************
name            : multiply_mtx_TMV
role		: multiply a vector to a matrix.
@param		: sparse * pM the matrix to be used
@param		: double *vec: The operand vector.
@param		: double *res: The resulting vector.
@return         : err_OK if everything went fine; err_ERROR otherwise.
remark		: size should be correct.
******************************************************************************/
err_state multiply_mtx_TMV(/*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * vec,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * res)
{
        int rows;

        if ( NULL == pM || NULL == vec || NULL == res ) {
                /*@-mustdefine@*/
                err_msg_5(err_PARAM, "multiply_mtx_TMV(%p[%dx%d],%p,%p)",
                                (const void*) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, (const void*)vec,
                                (void *) res, err_ERROR);
                /*@=mustdefine@*/
        }

        rows = mtx_rows(pM);
	/*Set the resulting array to zero values*/
        memset(res, (int) '\0', sizeof(double) * rows);
	/*Compute*/
        mtx_walk_all(pM, row, col, val)
	{
                res[col] += vec[row] * val;
	}
        end_mtx_walk_all;
        return err_OK;
}

/*****************************************************************************
name		: multiply_mtx_cer_TMV
role		: multiply a vector to a matrix.
@param		: sparse * pM the matrix to be used
@param          : double * vec: The operand vector.
@param          : double * res: The resulting vector.
@param		: int * pIds the valid ids, pIds[0] contains the number of ids
@return         : err_OK if everything went fine; err_ERROR otherwise.
remark		: size should be correct.
******************************************************************************/
err_state multiply_mtx_cer_TMV(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * vec,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * res,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pIds)
{
	/*ToDo: Implement*/
        int i, total_length, length, id;
        double vec_i;

        if ( NULL==pM || NULL==vec || NULL==res || NULL==pIds || 0 > pIds[0] ) {
                /*@-mustdefine@*/
                err_msg_7(err_PARAM,
                                "multiply_mtx_cer_TMV(%p[%dx%d],%p,%p,%p[%d])",
                                (const void*) pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, (const void*)vec,
                                (void *) res, (const void *) pIds,
                                NULL != pIds ? pIds[0] : 0, err_ERROR);
                /*@=mustdefine@*/
        }

        total_length = mtx_rows(pM);
        length = pIds[0];
	/*Set the resulting array to zero values*/
        memset(res, (int) '\0', sizeof(double) * total_length);

	/*Compute*/
	for( i=1; i <= length; i++ )
	{
		id = pIds[i];		/*Get the valid row id*/
		vec_i = vec[id];	/*The corresponding element in the vec vector*/
                mtx_walk_row(pM, id, col, val) {
                        res[col] += vec_i * val;
                } end_mtx_walk_row;
	}
        return err_OK;
}

/************************************************************************/
/**************SUPPLEMENTARTY METHODS USED FOR PRINTING ETC**************/
/************************************************************************/

err_state print_pattern_vec_double(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const char * pattern,
                int length, /*@observer@*/ /*@null@*/ const double * vec)
{
	int i;

        if ( NULL == pattern || 0 > length ) {
                err_msg_3(err_PARAM, "print_pattern_vec_double(\"%s\",%d,%p)",
                                NULL == pattern ? "NULL" : pattern, length,
                                (const void *) vec, err_ERROR);
        }

	if( vec != NULL){
		printf("( ");
		for( i = 0 ; i < length ; i++ )
		{
			printf( pattern, vec[i] );

			if (i != length-1) printf(", ");
		}
		printf(" )");
	}else{
		printf("NULL");
	}
        return err_OK;
}

err_state print_vec_double(int length,
                /*@observer@*/ /*@null@*/ const double * vec)
{
        if ( 0 > length ) {
                err_msg_2(err_PARAM, "print_vec_double(%d,%p)", length,
                                (const void *) vec, err_ERROR);
        }
        if ( err_state_iserror(print_pattern_vec_double("%e", length, vec)) ) {
                err_msg_2(err_CALLBY, "print_vec_double(%d,%p)", length,
                                (const void *) vec, err_ERROR);
        }
        return err_OK;
}

err_state print_vec_int(int length,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * vec)
{
	int i;

        if ( 0 > length || NULL == vec ) {
                err_msg_2(err_PARAM, "print_vec_int(%d,%p)", length,
                                (const void *) vec, err_ERROR);
        }

	printf("( ");
	for( i = 0 ; i < length ; i++ )
	{
		printf( "%d", vec[i] );
		if (i != length-1) printf(", ");
	}
	printf(" )\n");
        return err_OK;
}

/* The function is never used. Therefore, I have commented it out. David N.
   Jansen.
err_state print_vec_short(int length,
                / *@observer@* / / *@null@* / const short * vec)
{
	int i;
        if ( 0 > length || NULL == vec ) {
                err_msg_2(err_PARAM, "print_vec_short(%d,%p)", length,
                                (const void *) vec, err_ERROR);
        }
	if( vec != NULL){
		printf("( ");
		for( i = 0 ; i < length ; i++ )
		{
                        printf("%hd", vec[i]);
			if (i != length-1) printf(", ");
		}
		printf(" )\n");
	}else{
		printf("NULL\n");
	}
        return err_OK;
}
*/

/**
* Write the matrix to a file in tra-format
* @param pM the matrix to work with
* @param fname: the name of the output file
* @return 0: if writing failed, 1: for success
* NOTE: This function returns 0 if a sparse matrix is not valid, or the file
*	could not be opened for writing. This is also indicated via a printf statement.
*/
/* The function is never used. Therefore, I have commented it out. David N.
   Jansen.
err_state write_tra_file(/ *@observer@* // *@i1@* // *@null@* /const sparse *pM,
                / *@observer@* / / *@i1@* / / *@null@* / const char * fname)
{
        int i, tn;
	FILE *fp;

        if ( NULL == pM || NULL == fname || '\0' == fname[0] ) {
                err_msg_4(err_PARAM, "write_tra_file(%p[%dx%d],\"%s\")",
                                (const void *)pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0,
                                NULL == fname ? "NULL" : fname, err_ERROR);
        }

	/ * First we want to calculate the number of transitions. * /
	tn = 0;
	for(i=0;i<pM->rows;i++) {
                tn += mtx_next_num(pM, i);
                / *@-realcompare@* /
                if ( 0.0 != mtx_get_diag_val_nt(pM, i) ) / *@=realcompare@* /
			tn++;
	}

	printf("Writing transition file to %s.\n", fname);
        errno = 0;
	fp = fopen(fname, "w");
        if( NULL == fp || 0 != errno ) {
                err_msg_4(err_FILE, "write_tra_file(%p[%dx%d],\"%s\")",
                                (const void *) pM, mtx_rows(pM), mtx_cols(pM),
                                fname, err_ERROR);
	}

        fprintf(fp, "STATES %d\n", mtx_rows(pM));
	fprintf(fp, "TRANSITIONS %d\n", tn);

	/ * Now, print all the transitions. There is one catch, namely,
	 * we have to make sure we insert diag at the right time, if
	 * needed.
         * /
        i = 0;
        do {
                mtx_walk_row_sorted(pM, (const state_index) i, target, val) {
                        fprintf(fp, "%d %d %.10f\n", i + 1, target + 1, val);
		}
                end_mtx_walk_row_sorted;
	}
        while ( ++i < mtx_rows(pM) );
        if ( 0 != fclose(fp) || 0 != errno ) {
                err_msg_4(err_FILE, "write_tra_file(%p[%dx%d],\"%s\")",
                                (const void *) pM, mtx_rows(pM), mtx_cols(pM),
                                fname, err_ERROR);
        }
        return err_OK;
}
*/

/*****************************************************************************
name            : print_mtx_sparse
role		: print the matrix.
@param		: sparse * pM the matrix to work with
@return	: int: 0: fail, 1: success
remark		: return 0: sparse matrix is not valid.
******************************************************************************/
err_state print_mtx_sparse(/*@observer@*/ /*@i1@*/ /*@null@*/ const sparse * pM)
{
        int i;

        if ( NULL == pM ) {
                err_msg_3(err_PARAM, "print_mtx_sparse(%p[%dx%d])",
                                (const void *)pM, NULL != pM ? mtx_rows(pM) : 0,
                                NULL != pM ? mtx_cols(pM) : 0, err_ERROR);
        }

        i = 0;
        do
	{
                int col = mtx_next_num(pM, i);
		printf("ncols[%d]:: %d\n", i, col);
                mtx_walk_row_nodiag(pM, (const int) i, coll, val) {
                        printf("[%d][%d]=%f\n", i, coll, val);
                } end_mtx_walk_row_nodiag;
	}
        while ( ++i < mtx_rows(pM) );

        i = 0;
        do {
                printf("diag:: [%d][%d]=%f\n", i, i, mtx_get_diag_val_nt(pM,i));
        } while ( ++i < mtx_rows(pM) );
        if ( mtx_rows(pM) == mtx_cols(pM) )
	{
		printf("********************* Printing Back Sets *************************************\n");
                i = 0;
                do
		{
			printf("state: %d : BS : [ ", i+1);
                        mtx_walk_column_nodiag_noval(pM, row,
                                                (const state_index) i)
                        {
                                printf("%d, ", row + 1);
                        } end_mtx_walk_column_nodiag_noval;
			printf("]\n");
		}
                while ( ++i < mtx_rows(pM) );
		printf("******************************************************************************\n");
	}
        return err_OK;
}

/*****************************************************************************
name		: split_A_into_DI_LU
role		: this method splits the given A matrix into two matrixes DI and L+U
		: where DI is an inverted diagonal matrix and L+U are all the rest elements of
		  the initial matrix A.
@param		: sparse *pA the initial matrix A
@param		: sparse *pDI the empty matrix which will be filled with the D elements
@param		: sparse *pLU the empty matrix which will be filled with the LU elements
*****************************************************************************/
err_state split_A_into_DI_LU(/*@observer@*/ /*@i1@*/ /*@null@*/ const sparse*pA,
                /*@i1@*/ /*@null@*/ sparse * pDI,
                /*@i1@*/ /*@null@*/ sparse * pLU)
{
	int i;

        if ( NULL == pA || NULL == pDI || NULL == pLU
                || mtx_rows(pA)!=mtx_rows(pDI)||mtx_cols(pA)!=mtx_cols(pDI)
                || mtx_rows(pA)!=mtx_rows(pLU)||mtx_cols(pA)!=mtx_cols(pLU))
        {
                err_msg_9(err_PARAM, "split_A_into_DI_LU(%p[%dx%d],%p[%dx%d],"
                                "%p[%dx%d])", (const void *) pA,
                                NULL != pA ? mtx_rows(pA) : 0,
                                NULL != pA ? mtx_cols(pA) : 0, (void *) pDI,
                                NULL != pDI ? mtx_rows(pDI) : 0,
                                NULL != pDI ? mtx_cols(pDI) : 0, (void *) pLU,
                                NULL != pLU ? mtx_rows(pLU) : 0,
                                NULL != pLU ? mtx_cols(pLU) : 0, err_ERROR);
        }

        i = mtx_rows(pA);
        do
	{
                double diag;
                i--;
		/*Copy values to the D matrix*/
                diag = mtx_get_diag_val_nt(pA, i);
                if ( /*@-realcompare@*/ 0.0 != diag /*@=realcompare@*/ ) {
                        mtx_set_diag_val_nt(pDI, i, 1 / diag);
                }

		/*Copy values to the LU matrix*/
                if ( err_state_iserror(mtx_copy_row(pLU, i, pA)) ) {
                        err_msg_9(err_CALLBY, "split_A_into_DI_LU(%p[%dx%d],"
                                "%p[%dx%d],%p[%dx%d])", (const void *) pA,
                                mtx_rows(pA), mtx_cols(pA), (void *) pDI,
                                mtx_rows(pDI), mtx_cols(pDI), (void *) pLU,
                                mtx_rows(pLU), mtx_cols(pLU), err_ERROR);
                }
                mtx_set_diag_val_nt(pLU, i, 0.0);
	}
        while ( 0 < i );
        return err_OK;
}

/*****************************************************************************
name		: split_A_into_D_LU
role		: this method splits the given A matrix into two matrixes DI and L+U
		: where D is a diagonal matrix and L+U are all the rest elements of
		  the initial matrix A.
@param		: sparse *pA the initial matrix A
@param		: double *pD the empty vector which will be filled with the diagonal elements
@param		: sparse *pLU the empty matrix which will be filled with the LU elements
*****************************************************************************/
err_state split_A_into_D_LU(/*@observer@*/ /*@i1@*/ /*@null@*/ const sparse *pA,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * pD,
                /*@i1@*/ /*@null@*/ sparse * pLU)
{
	int i;

        if ( NULL==pA || NULL==pD || NULL==pLU || mtx_rows(pA)!=mtx_rows(pLU)
                                        || mtx_cols(pA)!=mtx_cols(pLU) )
        {
                /*@-mustdefine@*/
                err_msg_7(err_PARAM,"split_A_into_D_LU(%p[%dx%d],%p,%p[%dx%d])",
                                (const void *)pA, NULL != pA ? mtx_rows(pA) : 0,
                                NULL != pA ? mtx_cols(pA) : 0, (void *) pD,
                                (void *) pLU, NULL != pLU ? mtx_rows(pLU) : 0,
                                NULL != pLU ? mtx_cols(pLU) : 0, err_ERROR);
                /*@=mustdefine@*/
        }

        i = mtx_rows(pA);
        do
	{
                --i;
		/*Copy values to the D vector*/
                pD[i] = mtx_get_diag_val_nt(pA, i);

		/*Copy values to the LU matrix*/
                if ( err_state_iserror(mtx_copy_row(pLU, i, pA)) ) {
                        err_msg_7(err_PARAM, "split_A_into_D_LU(%p[%dx%d],%p,"
                                "%p[%dx%d])", (const void *) pA, mtx_rows(pA),
                                mtx_cols(pA), (void *) pD, (void *) pLU,
                                mtx_rows(pLU), mtx_cols(pLU), err_ERROR);
                }
                mtx_set_diag_val_nt(pLU, i, 0.0);
	}
        while ( i > 0 );
        return err_OK;
}

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
*****************************************************************************/
err_state split_A_into_D_L_U(/*@observer@*/ /*@i1@*/ /*@null@*/ const sparse*pA,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * pD,
                /*@i1@*/ /*@null@*/ sparse * pL,
                /*@i1@*/ /*@null@*/ sparse * pU)
{
        BOOL isRFirst;
	int i, j, theNColsL, theNColsU, col;

        if ( NULL == pA || NULL == pD || NULL == pL || NULL == pU
                      || mtx_rows(pA)!=mtx_rows(pL)||mtx_cols(pA)!=mtx_cols(pL)
                      || mtx_rows(pA)!=mtx_rows(pU)||mtx_cols(pA)!=mtx_cols(pU))
        {
                /*@-mustdefine@*/
                err_msg_10(err_PARAM, "split_A_into_D_L_U"
                                "(%p[%dx%d],%p,%p[%dx%d],%p[%dx%d])",
                                (const void *)pA, NULL != pA ? mtx_rows(pA) : 0,
                                NULL != pA ? mtx_cols(pA) : 0, (void *) pD,
                                (void *) pL, NULL != pL ? mtx_rows(pL) : 0,
                                NULL != pL ? mtx_cols(pL) : 0, (void *) pU,
                                NULL != pU ? mtx_rows(pU) : 0,
                                NULL != pU ? mtx_cols(pU) : 0, err_ERROR);
                /*@=mustdefine@*/
        }

        i = mtx_rows(pA);
        do
	{
                --i;
		/*Copy values to the pD vector*/
                pD[i] = mtx_get_diag_val_nt(pA, i);

		/*Initialize the counter for the number of row elements of the pL and pU matrixes*/
		theNColsL = 0;
		theNColsU = 0;

                isRFirst = TRUE;

                /*Store possible Low diagonal values. Even if there are no
                values, then the number of columns will be zero and thus this
                will not
		affect the entire process*/
                pL->valstruc[i].col = pA->valstruc[i].col;
                pL->valstruc[i].val = pA->valstruc[i].val;

		/*Copy values to the pL, pU matrixes*/
                for ( j = 0 ; j < /*@-compmempass@*/ mtx_next_num(pA, i)
                                                /*@=compmempass@*/ ; j++ )
		{
			/*Get the column index of the next element in a row*/
                        col = pA->valstruc[i].col[j];
			/*Check and store element in pL or in pU*/
			if( col < i )
			{
				theNColsL++;
			}
			else
			{
				/* Case of col > i*/
				if( isRFirst )
				{
					/*Copy to the pU matrix*/
                                        pU->valstruc[i].col =
                                                        &pA->valstruc[i].col[j];
                                        pU->valstruc[i].val =
                                                        &pA->valstruc[i].val[j];
					isRFirst = 0;
				}
				theNColsU++;
			}
		}
		pL->valstruc[i].ncols = theNColsL;
		pU->valstruc[i].ncols = theNColsU;
	}
        while ( 0 < i );
        /*@-compmempass@*/
        return err_OK;
        /*@=compmempass@*/
}

/*****************************************************************************
name		: sum_vec_vec
role		: This method is used to compute the sum of two vectors A+B=Res.
@param		: int length the length of vectors
@param		: double *pA vector A
@param		: double *pB vector B
@param		: double *pRes vector Res
*****************************************************************************/
err_state sum_vec_vec(int length,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pA,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pB,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * pRes)
{
	int i;

        if ( 0 > length || NULL == pA || NULL == pB || NULL == pRes ) {
                /*@-mustdefine@*/
                err_msg_4(err_PARAM, "sum_vec_vec(%d,%p,%p,%p)", length,
                                (const void *) pA, (const void *) pB,
                                (void *) pRes, err_ERROR);
                /*@=mustdefine@*/
        }

	for( i = 0; i < length; i++ )
	{
		pRes[i] = pB[i] + pA[i];
	}
        /*@-mustdefine@*/ /* If length == 0, nothing should be written */
        return err_OK; /*@=mustdefine@*/
}

/*****************************************************************************
name            : sum_cer_vec_vec
role		: This method is used to compute the sum of two vectors A+B=Res.
@param		: double *pA vector A
@param		: double *pB vector B
@param		: double *pRes vector Res
@param		: int * pIds the valid ids, pIds[0] contains the number of ids
*****************************************************************************/
err_state sum_cer_vec_vec(/*@observer@*/ /*@i1@*/ /*@null@*/ const double * pA,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pB,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * pRes,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pIds)
{
        int i, length, id;

        if ( NULL==pA || NULL==pB || NULL==pRes || NULL==pIds || 0 > pIds[0] ) {
                /*@-mustdefine@*/
                err_msg_5(err_PARAM, "sum_cer_vec_vec(%p,%p,%p,%p[%d])",
                                (const void *) pA, (const void *) pB,
                                (void *) pRes, (const void *) pIds,
                                NULL != pIds ? pIds[0] : 0, err_ERROR);
                /*@=mustdefine@*/
        }

        length = pIds[0];
	for( i = 1; i <= length ; i++ )
	{
		id = pIds[i];
		pRes[ id ] = pB[ id ] + pA[ id ];
	}
        /*@-mustdefine@*/ /* If length == 0, nothing should be written */
        return err_OK; /*@=mustdefine@*/
}

/*****************************************************************************
name		: multiply_inv_diag_D_V
role		: multiply v*diag(d)I. I.e. take vector of diagonal elements
		  create matrix diag(d) invert it and multiply by the vector v
		  from the left side.
@param		: double * pD: The diagonal elements.
@param		: double * pV: The vector to multiply with.
@param		: double * pR: The resulting vector.
@param		: int length : The length of vectors
******************************************************************************/
err_state multiply_inv_diag_D_V(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pD,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pV,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * pR,
                int length)
{
	int i;
	double diag;

        if ( NULL == pD || NULL == pV || NULL == pR || 0 > length ) {
                /*@-mustdefine@*/
                err_msg_4(err_PARAM, "multiply_inv_diag_D_V(%p,%p,%p,%d)",
                                (const void *) pD, (const void *) pV,
                                (void *) pR, length, err_ERROR);
                /*@=mustdefine@*/
        }

	for( i = 0; i < length; i++ )
	{
		diag = pD[i];
		if( diag != 0 )
		{
			pR[i] = pV[i]/diag;
		}
	}
        /*@-mustdefine@*/ /* If length == 0, nothing should be written */
        return err_OK; /*@=mustdefine@*/
}

/*****************************************************************************
name		: multiply_cer_inv_diag_D_V
role		: multiply v*diag(d)I. I.e. take vector of diagonal elements
		  create matrix diag(d) invert it and multiply by the vector v
		  from the left side.
@param		: double * pD: The diagonal elements.
@param		: double * pV: The vector to multiply with.
@param		: double * pR: The resulting vector.
@param		: int * pValidStates : The set of valid states
******************************************************************************/
err_state multiply_cer_inv_diag_D_V(
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pD,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const double * pV,
                /*@out@*/ /*@i1@*/ /*@null@*/ double * pR,
                /*@observer@*/ /*@i1@*/ /*@null@*/ const int * pValidStates)
{
	int i,idx;
	double diag;
        int length;

        if ( NULL == pD || NULL == pV || NULL == pR || NULL == pValidStates
                                                || 0 > pValidStates[0] )
        {
                /*@-mustdefine@*/
                err_msg_5(err_PARAM,
                                "multiply_cer_inv_diag_D_V(%p,%p,%p,%p[%d])",
                                (const void *) pD, (const void *) pV,
                                (void *) pR, (const void *) pValidStates,
                                NULL != pValidStates ? pValidStates[0] : 0,
                                err_ERROR);
                /*@=mustdefine@*/
        }

        length = pValidStates[0];
	for( i = 1; i <= length; i++ )
	{
		idx = pValidStates[i];
		diag = pD[idx];
		if( diag != 0 )
		{
			pR[idx] = pV[idx]/diag;
		}
	}
        /*@-mustdefine@*/ /* If length == 0, nothing should be written */
        return err_OK; /*@=mustdefine@*/
}
