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
*	Authors: Maneesh Khattri, Ivan Zapreev, Tim Kemna
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
*	Source description: Contains some common methods for transient
*				 analysis of PCTL/CSL/PRCTL/CSRL - X, U.
*/

#ifndef TRANSIENT_COMMON_H
#define TRANSIENT_COMMON_H

#include "bitset.h"
#include "sparse.h"

/**
* Solve E(phi U psi) until formula.
* @param: sparse *state_space: the state space
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: bitset *psi: satisfaction relation for psi formula.
* @return: bitset *: result of E(SAT(phi) U SAT(psi)) for all states.
* NOTE: adapted for t-bounded-until from
*         1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
extern
bitset * get_exist_until(const sparse *state_space, const bitset *phi, const bitset *psi);

/**
* Solve A(phi U psi) until formula.
* @param: sparse *state_space: the state space
* @param: bitset *phi: satisfaction relation for phi formula.
* @param: bitset *psi: satisfaction relation for psi formula.
* @param: bitset *e_phi_psi: The indicator set for the formula E(Phi U Psi)
* @return: bitset *: result of A(SAT(phi) U SAT(psi)) for all states.
* NOTE: adapted for t-bounded-until from
*         1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
extern
bitset * get_always_until(const sparse *state_space, const bitset *phi, const bitset *psi, const bitset *e_phi_psi);

/**
* Universal part of PCTL and CSL unbounded until
* Solves the system of linear equations Ax=b
* @param: pM: the A matrix
* @param: pValidStates: this array contains the number of nodes as the
*	    first element all the other elements are the node ids, if it is
*	    NULL then all the nodes from the A matrix are valid
* @param: pX: the initial x vector.
*	NOTE: Is possibly freed inside by solveGaussJacobi method
* @param: pB: the b vector
* @param: BOOL revert_matrix: revert the matrix afterwards, if true
* @return: the solution of the system
*/
extern double * unbounded_until_universal(sparse * pM, int * pValidStates,
                double * pX, const double * pB, BOOL revert_matrix);

/**
* Do optimization, we exclude states from Phi which are Psi states and
* also states from which you always go to bad i.e. not Phi ^ not Psi states
* NOTE: This optimization is applicable for until if only it has lower time bound 'subi'
* (and reward bound 'subj' if any) equal to 0. The latter is because Psi states should
* not be made absorbing before the lower time bound 'subi' is reached!
* So it should be used with time (and reward) bounded until, NOT interval until,
* though there it may be applied in a smart way as well.
* @param phi the Phi states
* @param psi the Psi states
* @param state_space this is the sparse matrix
* @return phi_and_not_psi without phi states from which you never reach psi states
*/
extern bitset * get_good_phi_states(const bitset * phi, const bitset * psi,
                const sparse * state_space);

/**
* Make certain states (not in n_absorbing) absorbing.
* Creates a new empty sparse matrix, then assighs non absorbing rows via pointers
* and computes uniformization rates for them.
* @param: sparse *state_space: the state space
* @param: bitset *n_absorbing: not absorbing states.
* @param: double *lambda(ref.): return value of uniformizing rate(also q).
* @param: double *abse: return row_sums
* @param: int    *non_absorbing(ref.): return value for the number of
*                       non-absorbing states in the model
* @return: sparse *: new sparse matrix with some states made absorbing.
* NOTE: makes states absorbing (in ab_state_space) by assigning pointers
*	  to rows (in state_space) only for those rows which are not to be
*	  made absorbing.
* NOTE: initially *lambda should be equal to 0!!!!
******************************************************************************/
extern
sparse * ab_state_space(const sparse * state_space, const bitset *n_absorbing, double *lambda, double *abse, int *non_absorbing);

/**
* Free absorbing sparse matrix.
* @param                : sparse *abs_local: the state space
* NOTE: Its row elements were copied as structures,
*	 so the internal row (values structure) pointers refer to the original matrix
*	 and thus should not be cleaned !!!
******************************************************************************/
extern void free_abs(sparse * abs_local);
#define free_abs(abs_local) (free_mtx_wrap((abs_local)))

#endif
