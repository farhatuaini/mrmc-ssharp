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
*	Source description: Perform Transient Analysis for CSL - X, U.
*/

#ifndef TRANSIENT_CTMC_H
#define TRANSIENT_CTMC_H

#include "bitset.h"

/**
* Solve the unbounded until operator for all states.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @return	: double *: result of the unbounded until operator for all states.
* NOTE: adapted from
*        1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*        In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*        eds.: Validation of stochastic systems.
*	 LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
extern
double * unbounded_until(const bitset *phi, const bitset *psi);

/**
* Solve the unbounded until operator for all states with lumping.
* @param: bitset *phi: SAT(phi).
* @param: bitset *psi: SAT(psi).
* @return: double *: result of the unbounded until operator for all states.
*/
extern
double * unbounded_until_lumping(const bitset *phi, const bitset *psi);

/**
* Solve the bounded until operator.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @param		: double supi: sup I
* @return	: double *: result of the unbounded until operator for all states.
* NOTE: 1. J.-P. Katoen, M. Kwiatkowska, G. Norman, D. Parker.
*         Faster and symbolic CTMC model checking. In: L. de Alfaro, S. Gilmore,
*         eds., Process algebra and probabilistic methods. LNCS Vol. 2165,
*         Springer, Berlin, pp. 23-38, 2001.
*/
extern
double * bounded_until(const bitset *phi, const bitset *psi, double supi);

/**
* Solve the bounded until operator with lumping.
* @param: bitset *phi: SAT(phi).
* @param: bitset *psi: SAT(psi).
* @param: double supi: sup I
* @return: double *: result of the unbounded until operator for all states.
*/
extern
double * bounded_until_lumping(const bitset *phi, const bitset *psi, double supi);

/**
* @param: bitset *phi: SAT(phi).
* @param: bitset *psi: SAT(psi).
* @param: double subi: sub I
* @param: double supi: sup I
* @return: double *: result of the interval until operator for all states.
* NOTE: 1. C.Baier, B.R. Haverkort, H. Hermanns and J.-P. Katoen.
*	Model Checking Algorithms for Continuous-Time Markov Chains.
*	IEEE Transactions on Software Engineering, Vol. 29, No. 6,
*	pp. 299-316, 2000.
*/
extern
double * interval_until(const bitset *phi, const bitset *psi, double subi, double supi);

/**
* Solve the interval until operator with lumping.
* @param: bitset *phi: SAT(phi).
* @param: bitset *psi: SAT(psi).
* @param: double subi: sub I
* @param: double supi: sup I
* @return: double *: result of the interval until operator for all states.
*/
extern
double * interval_until_lumping(const bitset *phi, const bitset *psi, double subi, double supi);

/**
* Solve the unbounded next operator.
* @param: bitset *phi: SAT(phi).
* @return: double *: result of the unbounded next operator for states.
*/
extern double * unbounded_next(const bitset * phi);

/**
* Solve the bounded next operator.
* @param: bitset *phi: SAT(phi).
* @param: double supi: sup I
* @return: double *: result of the bounded next operator for all states.
*/
extern double * bounded_next(const bitset * phi, double supi);

/**
* Solve the interval next operator.
* @param: bitset *phi: SAT(phi).
* @param: double subi: sub I
* @param: double supi: sup I
* @return: double *: result of the interval next operator for all states.
*/
extern double * interval_next(const bitset * phi, double subi, double supi);

#endif
