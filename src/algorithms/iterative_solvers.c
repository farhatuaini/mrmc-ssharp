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
*	Source description: This file contains functions for the iterative
*		solving of the systems of linear equations.
*/

#include "iterative_solvers.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
* This method normalizes solution i.e. makes the sum of its elements to be 1
* @param length the length of vector pSolution if pValidStates == NULL
* @param pSolution the solution vector
* @param pValidStates the valid states of the solution,
*                     pValidStates[0] == the number of valid states
*                     pValidStates == NULL if all states are valid
* @param in_any_case if TRUE then the normalization is performed even if sum < 1
*/
void normalizeSolution( const int length, double * pSolution,
                        const
			int * pValidStates, BOOL in_any_case )
{
	int i;
	double sum = 0;
        if ( NULL != pValidStates )
	{
		/*If not all states are valid*/
		int size = pValidStates[0];
		for( i = 1; i <= size ; i++ )
		{
			sum += pSolution[ pValidStates[i] ];
		}
		if( ( sum > 0 ) && (in_any_case || sum > 1) )
		{
			/*If we can make a solution smaller*/
			for( i = 1; i <= size ; i++ )
			{
				pSolution[ pValidStates[i] ] /= sum;
			}
		}
	}
	else
	{
		/*If all the states are valid*/
		for( i = 0; i < length ; i++ )
		{
			sum += pSolution[i];
		}
		if( ( sum > 0 ) && (in_any_case || sum > 1) )
		{
			/*If we can make a solution smaller*/
			for( i = 0; i < length ; i++ )
			{
				pSolution[i] /= sum;
			}
		}
	}
}

/**
* This method checks the convergence of solution
* @param pX1 the first vector
* @param pX2 the second vector
* @param err the allowed difference between corresponding vector elements
* @param pValidStates this array contains the number of nodes as the first
*		element
*               all the other elements are the node ids
* @return TRUE if vectors are close enough otherwise FALSE;
*/
static BOOL stopGJIC(const double * pX1, const double * pX2, double err,
                const int * pValidStates)
{
	BOOL result = TRUE;
	int i, id;
	const int size = pValidStates[0];
	for( i = 1; i <= size; i++ )
	{
		id = pValidStates[i];
		if( fabs( pX1[ id ]- pX2[ id ] ) > err )
		{
			result = FALSE;
			break;
		}
	}
	return result;
}

/**
* This method checks the convergence of solution
* @param size the length of vectors
* @param pX1 the first vector
* @param pX2 the second vector
* @param err the allowed difference between corresponding vector elements
* @return TRUE if vectors are close enough otherwise FALSE;
*/
static BOOL stopGJI(const int size, const double * pX1, const double * pX2,
                double err)
{
	BOOL result = TRUE;
	int i;
	for( i = 0; i < size; i++ )
	{
		if( fabs( pX1[ i ]- pX2[ i ] ) > err )
		{
			result = FALSE;
			break;
		}
	}
	return result;
}

/**
* If the pB vector is non zero result=DI(b+(-LU)x), and pValidStates are available
*/
static double * solveGaussJacobi_V_B(const sparse * pDI, const sparse * pLU,
                                        double * pIntermediate, double * pX,
                                        const double * pB, double err,
                                        int max_iterations,
                                        const int * pValidStates,
                                        const int N_STATES)
{
        double * pResult = (double *) calloc((size_t) N_STATES, sizeof(double));
	int i=0;
	double *pTmp;

	/* Copy the initial vector into the results. */
	/* pX is initialized to 1 for states that already fulfill the unbounded until formula or reach only  */
	/* states which fulfill the formula with probability 1, therefore no further probability has to be calculated */
	/* and for this states result = initial value */
	/* analog for states which reach only states fulfilling the formula with probability 0 */
	memcpy(pResult, pX, sizeof(double) * N_STATES);

	while(TRUE)
	{
		i++;
		/* intermediate = (-LU)*x */
                if ( err_state_iserror(multiply_mtx_cer_MV(pLU,pX,pIntermediate,
                                pValidStates[0], &pValidStates[1]))
		/* intermediate = b + intermediate */
                                || err_state_iserror(sum_vec_vec(N_STATES, pB,
                                                pIntermediate, pIntermediate))
		/* result = DI * intermediate */
                                || err_state_iserror(multiply_mtx_cer_MV(pDI,
                                                pIntermediate, pResult,
                                                pValidStates[0],
                                                &pValidStates[1])) )
                {
                        err_msg_13(err_CALLBY, "solveGaussJacobi_V_B(%p[%dx%d],"
                                "%p[%dx%d],%p,%p,%p,%g,%d,%p,%d)",
                                (const void *) pDI, mtx_rows(pDI),mtx_cols(pDI),
                                (const void *) pLU, mtx_rows(pLU),mtx_cols(pLU),
                                (void *) pIntermediate, (void *) pX,
                                (const void *) pB, err, max_iterations,
                                (const void *) pValidStates, N_STATES,
                                (free(pX), free(pResult), NULL));
                }

		/* Stop if we need or can(b+(-LU)x) */
		if( stopGJIC( pResult, pX, err, pValidStates ) || ( i > max_iterations ) )
			break;

		/* Switch values between pX and pResult */
		pTmp = pX;
		pX = pResult;
		pResult = pTmp;
	}
	free(pX);

	printf("Gauss Jacobi V_B: The number of Gauss-Jacobi iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

	return pResult;
}

/**
* If the pB vector is zero result=DI((-LU)x), and pValidStates are available
*/
static double * solveGaussJacobi_V_0(const sparse * pDI, const sparse * pLU,
                                        double * pIntermediate, double * pX,
                                        double err, int max_iterations,
                                        const int * pValidStates,
                                        const int N_STATES)
{
        double * pResult = (double *) calloc((size_t) N_STATES, sizeof(double));
	int i=0;
	double *pTmp;

	/* Copy the initial vector into the results. */
	/* pX is initialized to 1 for states that already fulfill the unbounded until formula or reach only  */
	/* states which fulfill the formula with probability 1, therefore no further probability has to be calculated */
	/* and for this states result = initial value */
	/* analog for states which reach only states fulfilling the formula with probability 0 */
	memcpy(pResult, pX, sizeof(double) * N_STATES);

	while(TRUE)
	{
		i++;
		/* intermediate = (-LU)*x */
                if ( err_state_iserror(multiply_mtx_cer_MV(pLU,pX,pIntermediate,
                                pValidStates[0], &pValidStates[1]))
		/* result =  DI * intermediate */
                                || err_state_iserror(multiply_mtx_cer_MV(pDI,
                                                pIntermediate, pResult,
                                                pValidStates[0],
                                                &pValidStates[1])) )
                {
                        err_msg_12(err_CALLBY, "solveGaussJacobi_V_0(%p[%dx%d],"
                                "%p[%dx%d],%p,%p,%g,%d,%p,%d)",
                                (const void *) pDI, mtx_rows(pDI),mtx_cols(pDI),
                                (const void *) pLU, mtx_rows(pLU),mtx_cols(pLU),
                                (void *) pIntermediate, (void *) pX, err,
                                max_iterations, (const void *) pValidStates,
                                N_STATES, (free(pX), free(pResult), NULL));
                }

		/* Stop if we need or can */
		if( stopGJIC( pResult, pX, err, pValidStates) || ( i > max_iterations ) )
			break;

		/* Switch values between pX and pResult */
		pTmp = pX;
		pX = pResult;
		pResult = pTmp;
	}
	free(pX);

	printf("Gauss Jacobi V_0: The number of Gauss-Jacobi iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

	return pResult;
}

/**
* Solves the system of linear equations Ax=b using the Gauss-Jacobi method
* @param pA the A matrix
* @param pX the initial x vector
*           NOTE: Is possibly freed inside by solveGaussJacobi method
* @param pB the b vector
* @param err the difference between two successive x vector values which
*	     indicates when we can stop iterations
* @param max_iterations the max number of iterations
* @param pValidStates this array contains the number of nodes as the first
*         element all the other elements are the node ids, if it is NULL then
*                     all the nodes from the A matrix are valid
* @return the solution of the system
*/
static double * solveGaussJacobi(const sparse *pA, double *pX, const double *pB,
                                  double err, int max_iterations,
                                  const int * pValidStates)
{
        sparse * pDI = NULL, * pLU = NULL;
        const int N_STATES = mtx_rows(pA);
	double *pResult = NULL;
        double * pIntermediate = (double *) calloc((size_t) N_STATES,
                        sizeof(double));

        const char * error_str = err_MEMORY;
        if ( NULL == pIntermediate
	/* Init the pDI matrix */
                        || (error_str = err_CALLBY,
                            pDI = allocate_sparse_matrix(N_STATES,mtx_cols(pA)))
                                        == NULL
	/* Init the pLU matrix */
                        || (pLU = allocate_sparse_matrix(N_STATES,mtx_cols(pA)))
                                        == NULL
	/* Get D inverted and L+U */
                        || err_state_iserror(split_A_into_DI_LU(pA, pDI, pLU))
	/* Turn LU into -LU */
                        || err_state_iserror(mult_mtx_const(pLU, -1.0)) )
        {
                err_msg_8(error_str, "solveGaussJacobi(%p[%dx%d],%p,%p,%g,%d,"
                        "%p)", (const void *) pA, mtx_rows(pA), mtx_cols(pA),
                        (void *) pX, (const void *) pB, err, max_iterations,
                        (const void *) pValidStates,
                        ((void) (NULL == pIntermediate || (free(pIntermediate),
                                NULL == pDI || (free_mtx_wrap(pDI),
                                        NULL == pLU || (free_mtx_wrap(pLU),
                                                FALSE)))),
                         NULL));
        }

	/*Solve system iteratively*/
        if ( NULL != pValidStates )
	{
                if ( NULL != pB )
		{
			/* If there is a non zero pB vector result=DI(b+(-LU)x) */
			pResult = solveGaussJacobi_V_B(pDI, pLU, pIntermediate, pX, pB, err,
                                                        max_iterations,
                                                        pValidStates, N_STATES);
		}
		else
		{
			/* If the pB vector is zero result=DI((-LU)x) */
			pResult = solveGaussJacobi_V_0(pDI, pLU, pIntermediate, pX, err, max_iterations,
                                                        pValidStates, N_STATES);
		}
	}
	else
	{
                if ( NULL != pB )
		{
			/*ToDo: Implement*/
			printf("ERROR: solveGaussJacobi_NULL_B is not implemented.\n");
		}
		else
		{
			/*ToDo: Implement*/
			printf("ERROR: solveGaussJacobi_NULL_0 is not implemented.\n");
		}
	}

	free(pIntermediate);
        if ( NULL == pResult
	/* Turn -LU into LU */
                        || err_state_iserror(mult_mtx_const(pLU, -1.0))

                        || (err_state_iserror(free_mtx_wrap(pDI))
                                        && (pDI = NULL, TRUE))
                        || (pDI = NULL, err_state_iserror(free_mtx_wrap(pLU))
                                        && (pLU = NULL, TRUE)) )
        {
                err_msg_8(err_CALLBY, "solveGaussJacobi(%p[%dx%d],%p,%p,%g,%d,"
                        "%p)", (const void *) pA, mtx_rows(pA), mtx_cols(pA),
                        (void *) pX, (const void*) pB, err, max_iterations,
                        (const void *) pValidStates,
                        ((void) (NULL == pLU || (free_mtx_wrap(pLU),
                                NULL == pDI || (free_mtx_wrap(pDI),
                                        FALSE))),
                        NULL));
        }

	return pResult;
}

/**
* Multiply nRow row of pL by the elem vector and add result to the pLUx[nRow]
*/
static void multiplyLRowByVecAndAddToLUx(const int nRow, const sparse * pL,
                const double * elem, double * pLUx)
{
        double result = 0.0;

        mtx_walk_row_nodiag(pL, (const int) nRow, i, val)
	{
                result += val * elem[i];
	}
        end_mtx_walk_row_nodiag;
	pLUx[nRow] += result;
}

/**
* This function is used to order an array of valid states
* @param pValidStates this array contains the number of nodes as the first element
*                     all the other elements are the node ids, if it is NULL then
*                     all the nodes from the A matrix are valid
*/
static int compare_int(const void *e1, const void *e2 )
{
        return * (const int *) e1 - * (const int *) e2;
}

static void orderValidStates(int *pValidStates)
{
        qsort(&pValidStates[1], (size_t) pValidStates[0], sizeof(int),
                        compare_int);
}

/**
* If the pB vector is non zero result=(b+(-LU)x))DI, and pValidStates are available
*/
static void solveGaussSeidel_V_B(const sparse * pL, const sparse * pU,
                        const double * pD, double * pX, const double * pB,
                        double err, int max_iterations, const int *pValidStates,
                        double * pLUx)
{
	int i=0, j=0, cur_id;
	double tmp;
	BOOL endOfIterations;
	const int LENGTH = pValidStates[0];

		/*If no state is valid, stop computation, initial vector is result*/
    if ( 0 == LENGTH ) {
        i = 1;
    } else {
	while(TRUE)
	{
		i++;

		/*Compute the non changeable part of each iteration*/
                if ( err_state_iserror(multiply_mtx_cer_MV(pU, pX, pLUx, LENGTH,
                                                        &pValidStates[1])) )
                {
                        exit(err_macro_13(err_CALLBY, "solveGaussSeidel_V_B(%p"
                                "[%dx%d],%p[%dx%d],%p,%p,%p,%g,%d,%p,%p)",
                                (const void *) pL, mtx_rows(pL), mtx_cols(pL),
                                (const void *) pU, mtx_rows(pU), mtx_cols(pU),
                                (const void *) pD, (void*) pX, (const void*) pB,
                                err, max_iterations, (const void*) pValidStates,
                                (void *) pLUx, EXIT_FAILURE));
                }
		endOfIterations = TRUE;

                /*Start computing the (Xi+1)th solution element by element
                Perform the first iteration here as it is simpler than others*/
		cur_id = pValidStates[1];
		if( cur_id != 0 ) multiplyLRowByVecAndAddToLUx(cur_id, pL, pX, pLUx);
		/*xi+1 = (b+(-LU)x)DI*/
		tmp = (pB[cur_id] - pLUx[cur_id]) / pD[cur_id];
		endOfIterations = endOfIterations && fabs( tmp - pX[ cur_id ] ) <= err;
		pX[cur_id] = tmp;
                /*Perform all the remaining iterations*/
		for( j = 2; j <= LENGTH; j++)
		{
			cur_id = pValidStates[j];
			multiplyLRowByVecAndAddToLUx(cur_id, pL, pX, pLUx);
			/*xi+1 = (b+(-LU)x)DI*/
			tmp = (pB[cur_id] - pLUx[cur_id]) / pD[cur_id];
			endOfIterations = endOfIterations && fabs( tmp - pX[ cur_id ] ) <= err ;
			/* printf("new = %1.15le - old = %1.15le == %1.15le\n",tmp,pX[ cur_id ],tmp - pX[ cur_id ]); */
			pX[cur_id] = tmp;
		}
		/*Stop if we need or can*/
		if( endOfIterations || ( i > max_iterations ) ) break;
	}
    }

	printf("Gauss Seidel V_B: The number of Gauss-Seidel iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);
}

/**
* If the pB vector is zero result=((-LU)x)DI, and pValidStates are available
*/
static void solveGaussSeidel_V_0(const sparse * pL, const sparse * pU,
                        const double * pD,
			double * pX, double err, int max_iterations,
                        const int * pValidStates, double * pLUx)
{
	int i=0, j=0, cur_id;
	double tmp;
	BOOL endOfIterations;
	const int LENGTH = pValidStates[0];

		/*If no state is valid, stop computation, initial vector is result*/
    if ( 0 == LENGTH ) {
        i = 1;
    } else {
	while(TRUE)
	{
		i++;

		/*Compute the non changeable part of each iteration*/
                if ( err_state_iserror(multiply_mtx_cer_MV(pU, pX, pLUx, LENGTH,
                                                        &pValidStates[1])) )
                {
                        exit(err_macro_12(err_CALLBY, "solveGaussSeidel_V_0(%p"
                                "[%dx%d],%p[%dx%d],%p,%p,%g,%d,%p,%p)",
                                (const void *) pL, mtx_rows(pL), mtx_cols(pL),
                                (const void *) pU, mtx_rows(pU), mtx_cols(pU),
                                (const void *) pD, (void *) pX, err,
                                max_iterations, (const void *) pValidStates,
                                (void *) pLUx, EXIT_FAILURE));
                }
		endOfIterations = TRUE;

                /*Start computing the (Xi+1)th solution element by element
                Perform the first iteration here as it is simpler than others*/
		cur_id = pValidStates[1];
		if( cur_id != 0 ) multiplyLRowByVecAndAddToLUx(cur_id, pL, pX, pLUx);
		/*xi+1 = ((-LU)x)DI*/
		tmp =  - pLUx[cur_id] / pD[cur_id];
		endOfIterations = endOfIterations && fabs( tmp - pX[ cur_id ] ) <= err;
		/* printf("new = %1.15le - old = %1.15le == %1.15le\n",tmp,pX[ cur_id ],tmp - pX[ cur_id ]); */
		pX[cur_id] = tmp;
                /*Perform all the remaining iterations*/
		for( j = 2; j <= LENGTH; j++)
		{
			cur_id = pValidStates[j];
			multiplyLRowByVecAndAddToLUx(cur_id, pL, pX, pLUx);
			/*xi+1 = ((-LU)x)DI*/
			tmp =  - pLUx[cur_id] / pD[cur_id];
			endOfIterations = endOfIterations && fabs( tmp - pX[ cur_id ] ) <= err ;
			/* printf("new = %1.15le - old = %1.15le == %1.15le\n",tmp,pX[ cur_id ],tmp - pX[ cur_id ]); */
			pX[cur_id] = tmp;
		}
		/*Stop if we need or can*/
		if( endOfIterations || ( i > max_iterations ) ) break;
	}
    }

	printf("Gauss Seidel V_0: The number of Gauss-Seidel iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);
}

/**
* Solves the system of linear equations Ax=b using the Gauss-Seidel method
* @param pA the A matrix
* @param pX the initial x vector
* @param pB the b vector
* @param err the difference between two successive x vector values which
*            indicates
*	     when we can stop iterations
* @param max_iterations the max number of iterations
* @param pValidStates this array contains the number of nodes as the first element
*                     all the other elements are the node ids, if it is NULL then
*                     all the nodes from the A matrix are valid
*                     It is not const int * because solveGaussSeidel may sort
*                     it.
* @return the solution of the system
*/
static double * solveGaussSeidel(const sparse *pA, double *pX, const double *pB,
			  double err, int max_iterations, int * pValidStates)
{
        int N_STATES = mtx_rows(pA);
        sparse * pL = NULL, * pU = NULL;
        double * pD = (double *) calloc((size_t) N_STATES,sizeof(double));
        double * pLUx = NULL;

        const char * error_str = err_MEMORY;
        if ( NULL == pD
                        || (pLUx = (double *) calloc((size_t) N_STATES,
                                        sizeof(double))) == NULL
	/* Init the pL matrix */
                        || (error_str = err_CALLBY,
                            pL = allocate_sparse_matrix(N_STATES, mtx_cols(pA)))
                                        == NULL
	/* Init the pU matrix */
                        || (pU = allocate_sparse_matrix(N_STATES, mtx_cols(pA)))
                                        == NULL
	/*Get D, L and U, A = L + diag(D) + U*/
                        || err_state_iserror(split_A_into_D_L_U(pA, pD,pL,pU)) )
        {
                err_msg_8(error_str, "solveGaussSeidel(%p[%dx%d],%p,%p,%g,%d,"
                        "%p)", (const void *) pA, mtx_rows(pA), mtx_cols(pA),
                        (void *) pX, (const void *) pB, err, max_iterations,
                        (void *) pValidStates,
                        ((void) (NULL == pD || (free(pD),
                                NULL == pLUx || (free(pLUx),
                                        NULL == pL || (free_mtx_wrap(pL),
                                                NULL==pU || (free_mtx_wrap(pU),
                                                        FALSE))))),
                        NULL));
        }

	/*solve system iteratively*/
        if ( NULL != pValidStates )
	{
		/*NOTE: The pValidStates states should be ordered, otherwise we will get a wrong solution
		I.e. for example instead of {5,9,3} there should be {3,5,9}*/
		orderValidStates(pValidStates);
                if ( NULL != pB )
		{
			/*If the pB vector is non zero result=(b+(-LU)x)DI*/
			solveGaussSeidel_V_B(pL, pU, pD, pX, pB, err, max_iterations,
                                                pValidStates, pLUx);
		}
		else
		{
			/*If the pB vector is zero result=(-LU)xDI*/
			solveGaussSeidel_V_0(pL, pU, pD, pX, err, max_iterations, pValidStates,
                                                pLUx);
		}
	}
	else
	{
                if ( NULL != pB )
		{
			/*If the pB vector is non zero result=(b+(-LU)x)DI*/
			printf("ERROR: solveGaussSeidel_NULL_B is not implemented.\n");
		}
		else
		{
			/*If the pB vector is zero result=(-LU)xDI*/
			printf("ERROR: solveGaussSeidel_NULL_0 is not implemented.\n");
		}
	}

	free(pD);
	free(pLUx);
        if ( err_state_iserror(free_mtx_wrap(pL))
                        || (err_state_iserror(free_mtx_wrap(pU))
                                        && (pU = NULL, TRUE)) )
        {
                err_msg_8(err_CALLBY, "solveGaussSeidel(%p[%dx%d],%p,%p,%g,%d,"
                        "%p)", (const void *) pA, mtx_rows(pA), mtx_cols(pA),
                        (void *) pX, (const void *) pB, err, max_iterations,
                        (void *) pValidStates,
                        ((void) (NULL == pU || (free_mtx_wrap(pU), FALSE)),
                         NULL));
        }

        /*Return the result as it is modified internally and finally contains
          the resulting solution*/
	return pX;
}

/**
* If there is a non zero pB vector result=(b+x(-LU))DI, and pValidStates are available
*/
static double *solveGaussJacobiInverted_V_B(const sparse *pLU, const double *pD,
                                           double * pIntermediate, double * pX,
                                           const double * pB, double err,
                                           int max_iterations, const
					   int * pValidStates, const int N_STATES)
{
        double * pResult = (double *) calloc((size_t) N_STATES, sizeof(double));
	int i=0;
	double * pTmp;
	while( TRUE )
	{
		i++;
		/*intermediate = x*(-LU)*/
                if ( err_state_iserror(multiply_mtx_cer_TMV( pLU, pX,
                                                pIntermediate, pValidStates))
		/*intermediate = b + intermediate*/
                                || err_state_iserror(sum_cer_vec_vec(pB,
                                                pIntermediate, pIntermediate,
                                                pValidStates))
		/*result = intermediate * DI*/
                                || err_state_iserror(multiply_cer_inv_diag_D_V(
                                                pD, pIntermediate, pResult,
                                                pValidStates)) )
                {
                        err_msg_11(err_CALLBY, "solveGaussJacobiInverted_V_B(%p"
                                "[%dx%d],%p,%p,%p,%p,%g,%d,%p,%d)",
                                (const void *) pLU, mtx_rows(pLU),mtx_cols(pLU),
                                (const void *) pD, (void *) pIntermediate,
                                (void *) pX, (const void *) pB, err,
                                max_iterations, (const void *) pValidStates,
                                N_STATES, (free(pX), free(pResult), NULL));
                }
		/*Stop if we need or can*/
		if(stopGJIC(pResult,pX,err,pValidStates)||(i>max_iterations))break;

		/*Improve the convergence
		if ( i % 10 == 0) normalizeSolution( 0, pResult, pValidStates, FALSE );

		Switch values between pX and pResult*/
		pTmp = pX;
		pX = pResult;
		pResult = pTmp;
	}
	free(pX);

	printf("Gauss Jacobi Inverted V_B: The number of Gauss-Jacobi iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

	return pResult;
}

/**
* If the pB vector is zero result=x(-LU)DI, and pValidStates are available
*/
static double *solveGaussJacobiInverted_V_0(const sparse *pLU, const double *pD,
                                           double * pIntermediate,
					   double * pX, double err, int max_iterations,
                                           const
					   int * pValidStates, const int N_STATES)
{
        double * pResult = (double *) calloc((size_t) N_STATES, sizeof(double));
	int i=0;
	double * pTmp;
	while(TRUE)
	{
		i++;
		/*intermediate = x*(-LU)*/
                if ( err_state_iserror(multiply_mtx_cer_TMV(pLU, pX,
                                                pIntermediate, pValidStates))
		/*result = intermediate * DI*/
                                || err_state_iserror(multiply_cer_inv_diag_D_V(
                                                pD, pIntermediate, pResult,
                                                pValidStates)) )
                {
                        err_msg_10(err_CALLBY, "solveGaussJacobiInverted_V_0(%p"
                                "[%dx%d],%p,%p,%p,%g,%d,%p,%d)",
                                (const void *) pLU, mtx_rows(pLU),mtx_cols(pLU),
                                (const void *) pD, (void *) pIntermediate,
                                (void *) pX, err, max_iterations,
                                (const void *) pValidStates, N_STATES,
                                (free(pX), free(pResult), NULL));
                }
		/*Stop if we need or can*/
		if(stopGJIC(pResult,pX,err,pValidStates)||(i>max_iterations))break;
		/*Improve the convergence
		if ( i % 10 == 0) normalizeSolution( 0, pResult, pValidStates, FALSE );
		Switch values between pX and pResult */
		pTmp = pX;
		pX = pResult;
		pResult = pTmp;
	}
	free(pX);

	printf("Gauss Jacobi Inverted V_0: The number of Gauss-Jacobi iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

	return pResult;
}

/**
* If there is a non zero pB vector result=(b+x(-LU))DI and pValidStates are unavailable
*/
static double * solveGaussJacobiInverted_NULL_B(const sparse * pLU,
                                              const double * pD,
                                              double *pIntermediate, double *pX,
                                              const double * pB, double err,
                                              int max_iterations,
					      const int N_STATES)
{
        double * pResult = (double *) calloc((size_t) N_STATES, sizeof(double));
	int i=0;
	double * pTmp;
	while( TRUE )
	{
		i++;
		/*intermediate = x*(-LU)*/
                if ( err_state_iserror(multiply_mtx_TMV(pLU, pX, pIntermediate))
		/*intermediate = b + intermediate*/
                                || err_state_iserror(sum_vec_vec(N_STATES, pB,
                                                pIntermediate, pIntermediate))
		/*result = intermediate * DI*/
                                || err_state_iserror(multiply_inv_diag_D_V(pD,
                                                pIntermediate, pResult,
                                                N_STATES)) )
                {
                        err_msg_10(err_CALLBY, "solveGaussJacobiInverted_NULL_B"
                                "(%p[%dx%d],%p,%p,%p,%p,%g,%d,%d)",
                                (const void *) pLU, mtx_rows(pLU),mtx_cols(pLU),
                                (const void *) pD, (void *) pIntermediate,
                                (void *) pX, (const void *) pB, err,
                                max_iterations, N_STATES,
                                (free(pX), free(pResult), NULL));
                }
		/*Stop if we need or can*/
		if( stopGJI( N_STATES, pResult, pX, err ) || ( i > max_iterations ) ) break;
		/*Improve the convergence
		if ( i % 10 == 0) normalizeSolution( N_STATES, pResult, NULL, FALSE );
		Switch values between pX and pResult*/
		pTmp = pX;
		pX = pResult;
		pResult = pTmp;
	}
	free(pX);

	printf("Gauss Jacobi Inverted NULL_B: The number of Gauss-Jacobi iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

	return pResult;
}

/**
* If the pB vector is zero result=x(-LU)DI and pValidStates are unavailable
*/
static double * solveGaussJacobiInverted_NULL_0(const sparse * pLU,
                                              const double * pD,
                                              double * pIntermediate,
					      double * pX, double err, int max_iterations,
					      const int N_STATES)
{
        double * pResult = (double *) calloc((size_t) N_STATES, sizeof(double));
	int i=0;
	double * pTmp;
	while( TRUE )
	{
		i++;
		/*intermediate = x*(-LU)*/
                if ( err_state_iserror(multiply_mtx_TMV(pLU, pX, pIntermediate))
		/*result = intermediate * DI*/
                                || err_state_iserror(multiply_inv_diag_D_V(pD,
                                                pIntermediate, pResult,
                                                N_STATES)) )
                {
                        err_msg_9(err_CALLBY, "solveGaussJacobiInverted_NULL_0"
                                "(%p[%dx%d],%p,%p,%p,%g,%d,%d)",
                                (const void *) pLU, mtx_rows(pLU),mtx_cols(pLU),
                                (const void *) pD, (void *) pIntermediate,
                                (void *) pX, err, max_iterations, N_STATES,
                                (free(pX), free(pResult), NULL));
                }
		/*Stop if we need or can*/
		if( stopGJI( N_STATES, pResult, pX, err ) || ( i > max_iterations )) break;
		/*Improve the convergence
		if ( i % 10 == 0) normalizeSolution( N_STATES, pResult, NULL, FALSE );
		Switch values between pX and pResult*/
		pTmp = pX;
		pX = pResult;
		pResult = pTmp;
	}
	free(pX);

	printf("Gauss Jacobi Inverted NULL_0: The number of Gauss-Jacobi iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

	return pResult;
}

/**
* Solves the system of linear equations xA=b using the Inverted Gauss-Jacobi method
* @param pA the A matrix
* @param pX the initial x vector
* @param pB the b vector
* @param err the difference between two successive x vector values which
*            indicates
*	     when we can stop iterations
* @param max_iterations the max number of iterations
* @param pValidStates this array contains the number of nodes as the first element
*                     all the other elements are the node ids, if it is NULL then
*                     all the nodes from the A matrix are valid
* @return the solution of the system
*/
static double * solveGaussJacobiInverted(const sparse * pA, double * pX,
                                  const double * pB, double err,
                                  int max_iterations, const int * pValidStates)
{
        const int N_STATES = mtx_rows(pA);
        sparse * pLU = NULL;
        double * pD = (double *) calloc((size_t) N_STATES, sizeof(double));
	double * pResult = NULL;
        double * pIntermediate = NULL;

        const char * error_str = err_MEMORY;
        if ( NULL == pD
                        || (pIntermediate = (double *) calloc((size_t) N_STATES,
                                        sizeof(double))) == NULL
	/*Init the pLU matrix*/
                        || (error_str = err_CALLBY,
                            pLU = allocate_sparse_matrix(N_STATES,mtx_cols(pA)))
                                        == NULL
	/*Get D and L+U, A = L + diag(D) + U*/
                        || err_state_iserror(split_A_into_D_LU(pA, pD, pLU)) )
        {
                err_msg_8(error_str, "solveGaussJacobiInverted(%p[%dx%d],%p,%p,"
                        "%g,%d,%p)", (const void*)pA, mtx_rows(pA),mtx_cols(pA),
                        (void *) pX, (const void *) pB, err, max_iterations,
                        (const void *) pValidStates,
                        ((void) (NULL == pD || (free(pD),
                                NULL == pIntermediate || (free(pIntermediate),
                                        NULL == pLU || (free_mtx_wrap(pLU),
                                                FALSE)))),
                        NULL));
        }
	/*Solve system iteratively*/
        if ( NULL != pValidStates )
	{
		/*If not all states are valid
		Turn LU into -LU*/
                if ( err_state_iserror(mult_mtx_cer_const(pLU, -1.0,
                                                        pValidStates)) )
                {
                        err_msg_8(err_CALLBY, "solveGaussJacobiInverted(%p[%dx"
                                "%d],%p,%p,%g,%d,%p)", (const void *) pA,
                                mtx_rows(pA), mtx_cols(pA), (void *) pX,
                                (const void *) pB, err, max_iterations,
                                (const void*) pValidStates, (free_mtx_wrap(pLU),
                                free(pIntermediate), free(pD), NULL));
                }
                if ( NULL != pB )
		{
			/*If there is a non zero pB vector result=(b+x(-LU))DI*/
			pResult = solveGaussJacobiInverted_V_B(pLU, pD, pIntermediate,
						     pX, pB, err, max_iterations, pValidStates, N_STATES);
		}
		else
		{
			/*If the pB vector is zero result=x(-LU)DI*/
			pResult = solveGaussJacobiInverted_V_0(pLU, pD, pIntermediate,
						     pX, err, max_iterations, pValidStates, N_STATES);
		}
                if ( NULL == pResult
		/*Turn -LU into LU*/
                                || err_state_iserror(mult_mtx_cer_const(pLU,
                                                -1.0, pValidStates)) )
                {
                        err_msg_8(err_CALLBY, "solveGaussJacobiInverted(%p[%dx"
                                "%d],%p,%p,%g,%d,%p)", (const void *) pA,
                                mtx_rows(pA), mtx_cols(pA), (void *) pX,
                                (const void *) pB, err, max_iterations,
                                (const void*) pValidStates,
                                (free(pResult), free_mtx_wrap(pLU),
                                free(pIntermediate), free(pD), NULL));
                }
	}
	else
	{
		/*If all states are valid
		Turn LU into -LU*/
                if ( err_state_iserror(mult_mtx_const(pLU, -1.0)) ) {
                        err_msg_8(err_CALLBY, "solveGaussJacobiInverted(%p[%dx"
                                "%d],%p,%p,%g,%d,%p)", (const void *) pA,
                                mtx_rows(pA), mtx_cols(pA), (void *) pX,
                                (const void *) pB, err, max_iterations,
                                (const void*) pValidStates, (free_mtx_wrap(pLU),
                                free(pIntermediate), free(pD), NULL));
                }
                if ( NULL != pB )
		{
			/*If there is a non zero pB vector result=(b+x(-LU))DI*/
			pResult = solveGaussJacobiInverted_NULL_B(pLU, pD, pIntermediate,
							pX, pB, err, max_iterations, N_STATES);
		}
		else
		{
			/*If the pB vector is zero result=x(-LU)DI*/
			pResult = solveGaussJacobiInverted_NULL_0(pLU, pD, pIntermediate,
							pX, err, max_iterations, N_STATES);
		}
                if ( NULL == pResult
		/*Turn -LU into LU*/
                                || err_state_iserror(mult_mtx_const(pLU,-1.0)) )
                {
                        err_msg_8(err_CALLBY, "solveGaussJacobiInverted(%p[%dx"
                                "%d],%p,%p,%g,%d,%p)", (const void *) pA,
                                mtx_rows(pA), mtx_cols(pA), (void *) pX,
                                (const void *) pB, err, max_iterations,
                                (const void*) pValidStates,
                                (free(pResult), free_mtx_wrap(pLU),
                                free(pIntermediate), free(pD), NULL));
                }
	}

	free(pIntermediate);
	free(pD);
        if ( err_state_iserror(free_mtx_wrap(pLU)) ) {
                err_msg_8(err_CALLBY, "solveGaussJacobiInverted(%p[%dx%d],%p,"
                                "%p,%g,%d,%p)", (const void *) pA,
                                mtx_rows(pA), mtx_cols(pA), (void *) pX,
                                (const void *) pB, err, max_iterations,
                                (const void *) pValidStates,
                                (free(pResult), NULL));
        }

	return pResult;
}

/*
* Multiply nRow row of pU by the elem value and add resulting vector to the pLUx
*/
static void multiplyUrowByConstAndAddToLUx(const int nRow, const sparse * pU,
                double elem, double * pLUx)
{
        mtx_walk_row_nodiag(pU, (const int) nRow, i, val)
	{
                pLUx[i] += elem * val;
	}
        end_mtx_walk_row_nodiag;
}

/**
* If the pB vector is zero result=x(-LU)DI, and pValidStates are available
* NOTE: The pValidStates states should be ordered, otherwise we will get a wrong solution
*       I.e. for example instead of {5,9,3} there should be {3,5,9}
*/
static void solveGaussSeidelInverted_V_0(const sparse * pL, const sparse * pU,
                                                const double * pD,
						double * pX, double err, int max_iterations,
                                                const
						int * pValidStates, double * pLUx)
{
	int i=0, j=0, prev_id, cur_id;
	double tmp;
	BOOL endOfIterations;
	const int LENGTH = pValidStates[0];
	/* print_mtx_sparse(pU); */
	/* print_mtx_sparse(pL); */
        /* print_vec_double(mtx_rows(pU), pD); */
        /* print_vec_double(mtx_rows(pU), pX); */
	while(TRUE)
	{
		i++;
		/*Compute the non changeable part of each iteration*/
                if ( err_state_iserror(multiply_mtx_cer_TMV(pL, pX, pLUx,
                                                pValidStates)) )
                {
                        exit(err_macro_12(err_CALLBY,
                                "solveGaussSeidelInverted_V_0(%p[%dx%d],%p[%dx"
                                "%d],%p,%p,%g,%d,%p,%p)", (const void *) pL,
                                mtx_rows(pL), mtx_cols(pL), (const void *) pU,
                                mtx_rows(pU), mtx_cols(pU), (const void *) pD,
                                (void *) pX, err, max_iterations,
                                (const void *) pValidStates, (void *) pLUx,
                                EXIT_FAILURE));
                }

		endOfIterations = TRUE;

                /*Start computing the (Xi+1)th solution element by element
                Perform the first iteration here as it is simpler than others*/
		cur_id = pValidStates[1];
		tmp =  - pLUx[cur_id] / pD[cur_id];
		endOfIterations = endOfIterations && fabs( tmp - pX[ cur_id ] ) <= err;
		/* printf("new = %1.15le - old = %1.15le == %1.15le\n",tmp,pX[ cur_id ],tmp - pX[ cur_id ]); */
		pX[cur_id] = tmp;
		prev_id = cur_id;
                /*Perform all the remaining iterations*/
		for( j = 2; j <= LENGTH; j++)
		{
			cur_id = pValidStates[j];
			multiplyUrowByConstAndAddToLUx(prev_id, pU, pX[prev_id], pLUx);
			tmp =  - pLUx[cur_id] / pD[cur_id];
			endOfIterations = endOfIterations && fabs( tmp - pX[ cur_id ] ) <= err ;
			/* printf("new = %1.15le - old = %1.15le == %1.15le\n",tmp,pX[ cur_id ],tmp - pX[ cur_id ]); */
			pX[cur_id] = tmp;
			prev_id = cur_id;
		}
		/*Stop if we need or can*/
		if( endOfIterations || ( i > max_iterations ) ) break;
		/*Improve the convergence*/
		if ( i % 10 == 0) normalizeSolution( 0, pX, pValidStates, FALSE );
                /* print_vec_double(mtx_rows(pU), pX); */
	}

	printf("Gauss Seidel Inverted V_0: The number of Gauss-Seidel iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

}

/**
* If the pB vector is zero result=x(-LU)DI and pValidStates are unavailable
*/
static void solveGaussSeidelInverted_NULL_0(const sparse * pL, const sparse *pU,
                                                const double * pD,
						double * pX, double err, int max_iterations,
						const int N_STATES, double * pLUx)
{
	int i=0, j=0;
	double tmp;
	BOOL endOfIterations;
	while(TRUE)
	{
		i++;
		/*Compute the non changeable part of each iteration*/
                if ( err_state_iserror(multiply_mtx_TMV(pL, pX, pLUx)) ) {
                        exit(err_macro_12(err_CALLBY,
                                "solveGaussSeidelInverted_NULL_0(%p[%dx%d],%p"
                                "[%dx%d],%p,%p,%g,%d,%d,%p)", (const void *) pL,
                                mtx_rows(pL), mtx_cols(pL), (const void *) pU,
                                mtx_rows(pU), mtx_cols(pU), (const void *) pD,
                                (void *) pX, err, max_iterations,
                                N_STATES, (void *) pLUx, EXIT_FAILURE));
                }
		endOfIterations = TRUE;
                /*Start computing the (Xi+1)th solution element by element
                Perform the first iteration here as it is simpler than others*/
		tmp =  - pLUx[0] / pD[0];
		endOfIterations = endOfIterations && fabs( tmp - pX[ 0 ] ) <= err;
		pX[0] = tmp;
                /*Perform all the remaining iterations*/
		for( j = 1; j < N_STATES; j++)
		{
			multiplyUrowByConstAndAddToLUx(j-1, pU, pX[j-1], pLUx);
			tmp =  - pLUx[j] / pD[j];
			endOfIterations = endOfIterations && fabs( tmp - pX[ j ] ) <= err ;
			pX[j] = tmp;
		}
		/*Stop if we need or can*/
		if( endOfIterations || ( i > max_iterations ) ) break;
		/*Improve the convergence*/
		if ( i % 10 == 0) normalizeSolution( N_STATES, pX, NULL, FALSE );
	}

	printf("Gauss Seidel Inverted NULL_0: The number of Gauss-Seidel iterations %d\n",i);
	if( i > max_iterations ) printf("ERROR: %s\n", METHOD_DIVERGENCE_MSG);

}

/**
* Solves the system of linear equations xA=b using the Inverted Gauss-Seidel method
* @param pA the A matrix
* @param pX the initial x vector
* @param pB the b vector
* @param err the difference between two successive x vector values which
*            indicates
*	     when we can stop iterations
* @param max_iterations the max number of iterations
* @param pValidStates this array contains the number of nodes as the first element
*                     all the other elements are the node ids, if it is NULL then
*                     all the nodes from the A matrix are valid
*                     It is not const int * because solveGaussSeidelInverted may
*                     sort it.
* @return the solution of the system
*/
static double * solveGaussSeidelInverted(const sparse * pA, double * pX,
                                  const double * pB,
				  double err, int max_iterations, int * pValidStates)
{
        const int N_STATES = mtx_rows(pA);
        sparse * pL = NULL, * pU = NULL;
        double * pD = (double *) calloc((size_t) N_STATES, sizeof(double));
        double * pLUx = NULL;

        const char * error_str = err_MEMORY;
        if ( NULL == pD
                        || (pLUx = (double *) calloc((size_t) N_STATES,
                                        sizeof(double))) == NULL
	/*Init the pL matrix*/
                        || (error_str = err_CALLBY,
                            pL = allocate_sparse_matrix(N_STATES, mtx_cols(pA)))
                                        == NULL
	/*Init the pU matrix*/
                        || (pU = allocate_sparse_matrix(N_STATES, mtx_cols(pA)))
                                        == NULL
	/*Get D, L and U, A = L + diag(D) + U*/
                        || err_state_iserror(split_A_into_D_L_U(pA, pD,pL,pU)) )
        {
                err_msg_8(error_str, "solveGaussSeidelInverted(%p[%dx%d],%p,%p,"
                        "%g,%d,%p)", (const void*)pA, mtx_rows(pA),mtx_cols(pA),
                        (void *) pX, (const void *) pB, err, max_iterations,
                        (void *) pValidStates,
                        ((void) (NULL == pD || (free(pD),
                                NULL == pLUx || (free(pLUx),
                                        NULL == pL || (free_mtx_wrap(pL),
                                                NULL==pU || (free_mtx_wrap(pU),
                                                        FALSE))))),
                        NULL));
        }

	/*Solve system iteratively*/
        if ( NULL != pValidStates )
	{
		/*NOTE: The pValidStates states should be ordered, otherwise we will get a wrong solution
		I.e. for example instead of {5,9,3} there should be {3,5,9}*/
		orderValidStates(pValidStates);
                if ( NULL != pB )
		{
			/*ToDo: Implement*/
			printf("ERROR: solveGaussSeidelInverted_V_B is not implemented.\n");
		}
		else
		{
			/*If the pB vector is zero result=x(-LU)DI*/
			solveGaussSeidelInverted_V_0(pL, pU, pD, pX, err,
						max_iterations, pValidStates, pLUx);
		}
	}
	else
	{
                if ( NULL != pB )
		{
			/*ToDo: Implement*/
			printf("ERROR: solveGaussSeidelInverted_NULL_B is not implemented.\n");
		}
		else
		{
			/*If the pB vector is zero result=x(-LU)DI*/
			solveGaussSeidelInverted_NULL_0(pL, pU, pD, pX, err,
							max_iterations, N_STATES, pLUx);
		}
	}

	free(pD);
	free(pLUx);
        if ( err_state_iserror(free_mtx_wrap(pL))
                        || (err_state_iserror(free_mtx_wrap(pU))
                                        && (pU = NULL, TRUE)) )
        {
                err_msg_8(err_CALLBY, "solveGaussSeidelInverted(%p[%dx%d],%p,"
                        "%p,%g,%d,%p)", (const void *) pA, mtx_rows(pA),
                        mtx_cols(pA), (void *) pX, (const void *) pB, err,
                        max_iterations, (void *) pValidStates,
                        (free(pX),
                         (void) (NULL == pU || (free_mtx_wrap(pU), FALSE)),
                         NULL));
        }

        /*Return the result as it is modified internally and finally contains
          the resulting solution*/
	return pX;
}

/**
* Simple printing info method
*/
static void print_info(const char * method, double err, int max_iterations)
{
	printf("%s\n", method);
	printf("\tError: %e\n", err);
	printf("\tMax iter: %d\n", max_iterations);
}

/**
* Solves the system of linear equations xA=b using one of the Inverted methods
* @param type the type of the method to use (GAUSS_JACOBI,GAUSS_JACOBI_INV,GAUSS_SEIDEL,...)
* @param pA the A matrix
* @param pX the initial x vector
*           NOTE: Is possibly freed inside by a solver method (like solveGaussJacobi)
* @param pB the b vector if NULL then b = (0,...,0)
* @param err the difference between two successive x vector values which
*            indicates
*	     when we can stop iterations
* @param max_iterations the max number of iterations
* @param pValidStates this array contains the number of nodes as the first element
*                     all the other elements are the node ids, if it is NULL then
*                     all the nodes from the A matrix are valid
*                     It is not const int * because solve may sort it.
* @return the solution of the system
*/
static double * solve(int type, const sparse * pA, double *pX, const double *pB,
		double err, int max_iterations, int * pValidStates)
{
	double * result = NULL;
	printf("Solving the system of linear equations:\n\tMethod: ");
	switch(type)
	{
		case GAUSS_JACOBI:
			print_info("GAUSS-JACOBI", err, max_iterations);
			result = solveGaussJacobi(pA, pX, pB, err, max_iterations, pValidStates);
			break;
		case GAUSS_SEIDEL:
			print_info("GAUSS-SEIDEL", err, max_iterations);
			result = solveGaussSeidel(pA, pX, pB, err, max_iterations, pValidStates);
			break;
		case GAUSS_JACOBI_INV:
			print_info("GAUSS-JACOBI-INVERTED", err, max_iterations);
			result = solveGaussJacobiInverted(pA, pX, pB, err, max_iterations, pValidStates);
			break;
		case GAUSS_SEIDEL_INV:
			print_info("GAUSS-SEIDEL-INVERTED", err, max_iterations);
			result = solveGaussSeidelInverted(pA, pX, pB, err, max_iterations, pValidStates);
			break;
		default:
			printf("Bug: The method to solve a system of linear equations is not defined.\n");
                        exit(EXIT_FAILURE);
	}
	return result;
}

/**
* This method uses one of the method to iteratively solve the system of linear equations.
* The initial X value is computed inside this method
* @param pQ the matrix to be solved
* @param pX an optional initial vector, if not provided the one is created internally.
* @param pB the b vector if NULL then b = (0,...,0)
* @param err the difference between two successive x vector values which
*            indicates
*            when we can stop iterations
* @param max_iterations the max number of iterations
* @param pValidStates the valid states i.e. ids of the rows which are valid and for which
*                     the system is solved, this array contains the number of nodes as the
*                     first element
*                     It is not const int * because solve_initial may sort it.
*/
double * solve_initial(int type, const sparse * pQ, double * pX,
                                const double * pB, double err,
				int max_iterations, int * pValidStates)
{
	/* If a specific initial vector is not provided. */
	/* For example it is provided for the unbounded until operator (PCTL, CSL) */
	if(pX == NULL){
		double init_value;
		int i, length;

		/*Create the initial vector*/
                pX = (double *) calloc((size_t) mtx_rows(pQ), sizeof(double));
		init_value = 0.0;
		/*Initialize the initial X vector*/
                if ( NULL != pValidStates )
		{
			length = pValidStates[0];
			init_value = 1.0/length;
			for( i = 1; i <= length; i++)
			{
				pX[ pValidStates[i] ] = init_value;
			}
		}
		else
		{
                        length = mtx_rows(pQ);
			init_value = 1.0/length;
			for( i = 0; i < length; i++)
			{
				pX[i] = init_value;
			}
		}
	}

	/*Solve system of linear equations*/
	return solve(type, pQ, pX, pB, err, max_iterations, pValidStates);
}
