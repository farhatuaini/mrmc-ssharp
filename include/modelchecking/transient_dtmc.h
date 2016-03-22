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
*	Source description: Perform Transient Analysis for PCTL - X, U.
*/

#ifndef TRANSIENT_DTMC_H
#define TRANSIENT_DTMC_H

#include "bitset.h"

/**
* Solve the unbounded until operator for all states for DTMC.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @return	: double *: result of unbounded until operator for all states.
* NOTE: 1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
extern
double * dtmc_unbounded_until(const bitset *phi, const bitset *psi);

/**
* Solve the unbounded until operator for all states for DTMC with lumping.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @return		: double *: result of unbounded until operator for all states.
*/
extern
double * dtmc_unbounded_until_lumping(const bitset *phi, const bitset *psi);

/**
* Solve the bounded until operator for DTMC.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @param		: double supi: sup I should contain the Natural number
* @return	: double *: result of the unbounded until operator for all states.
* NOTE: 1. F. Ciesinski, M. Größer. On Probabilistic Computation Tree Logic.
*         In: C. Baier, B.R. Haverkort, H. Hermanns, J.-P. Katoen, M. Siegle,
*         eds.: Validation of stochastic systems.
*	  LNCS, Vol. 2925, Springer, pp. 147-188, 2004.
*/
extern
double * dtmc_bounded_until(const bitset *phi, const bitset *psi, double supi);

/**
* Solve the bounded until operator for DTMC with lumping.
* @param		: bitset *phi: SAT(phi).
* @param		: bitset *psi: SAT(psi).
* @param		: double supi: sup I should contain the Natural number
* @return		: double *: result of the unbounded until operator for all states.
*/
extern
double * dtmc_bounded_until_lumping(const bitset *phi, const bitset *psi, double supi);

/**
* This method modelchecks the Xphi formula for the DTMC
* @param the phi states
* @return the probabilities
*/
extern double * dtmc_unbounded_next(const bitset * phi);

#endif
