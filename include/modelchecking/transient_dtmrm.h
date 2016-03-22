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
*	Source description: Perform Transient Analysis for PRCTL - X, U.
*/

#ifndef TRANSIENT_DTMRM_H
#define TRANSIENT_DTMRM_H

#include "bitset.h"

/**
* Solve the N-J-bounded until operator for DTMRM.
* Adds optimization for the case subi == 0 & subj == 0
* @param phi SAT(phi).
* @param psi SAT(psi).
* @param subi lower rime bound
* @param supi upper time bound
* @param subj lower reward bound
* @param supj upper reward bound
* @return result of the N-J-bounded until operator for all states.
* NOTE: see (slightly modified path_graph):
*	1. S. Andova, H. Hermanns and J.-P. Katoen.
*	Discrete-time rewards model-checked.
*	FORMATS 2003, LNCS, Vol. 2791, Springer, pp. 88 - 103, 2003.
*/
extern double * dtmrm_bounded_until(const bitset * phi, const bitset * psi,
                double subi, double supi, double subj, double supj);

/**
* Solve the N-J-bounded until operator for DTMRM with formula dependent lumping.
* Adds optimization for the case subi == 0 & subj == 0
* NOTE: NO support for impulse rewards!!!
* @param phi SAT(phi).
* @param psi SAT(psi).
* @param subi lower rime bound
* @param supi upper time bound
* @param subj lower reward bound
* @param supj upper reward bound
* @return result of the N-J-bounded until operator for all states.
* NOTE: see (slightly modified path_graph):
*	1. S. Andova, H. Hermanns and J.-P. Katoen.
*	Discrete-time rewards model-checked.
*	FORMATS 2003, LNCS, Vol. 2791, Springer, pp. 88 - 103, 2003.
*/
extern double *dtmrm_bounded_until_lumping(const bitset *phi, const bitset *psi,
                double subi, double supi, double subj, double supj);

#endif
