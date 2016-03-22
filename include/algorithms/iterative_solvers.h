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
*	Source description: This file contains functions for the iterative
*		solving of the systems of linear equations.
*/

#ifndef ITERATIVE_SOLVERS_H
#define ITERATIVE_SOLVERS_H

#include "sparse.h"

#       define GAUSS_JACOBI 1
#       define GAUSS_SEIDEL 2
#       define GAUSS_JACOBI_INV 3
#       define GAUSS_SEIDEL_INV 4

#       define METHOD_DIVERGENCE_MSG "The results are UNRELIABLE!  The " \
                "numerical method did not converge. Please use another " \
                "method (see 'set method_path M' and 'set method_steady M') " \
                "or increase the maximum number of iterations (see 'set " \
                "max_iter I')."

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
        * @param pValidStates the valid states, i.e. ids of the rows which are
        *                     valid and for which the system is solved. This
        *                     array contains the number of nodes as the
	*                     first element
        *                     It is not const int * because solve_initial may
        *                     sort it.
	*/
        extern double * solve_initial(int type, const sparse * pQ, double * pX,
                                        const double * pB, double err,
					int max_iterations, int * pValidStates);

	/**
	* This method normalizes solution i.e. makes the sum of its elements to be 1
	* @param length the length of vector pSolution if pValidStates == NULL
	* @param pSolution the solution vector
	* @param pValidStates the valid states of the solution,
	*                     pValidStates[0] == the number of valid states
	*                     pValidStates == NULL if all states are valid
	* @param in_any_case if TRUE then the normalization is performed even if sum < 1
	*/
        extern void normalizeSolution(int length, double * pSolution,
                        const int * pValidStates, BOOL in_any_case);

#endif
