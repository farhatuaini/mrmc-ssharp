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
*	Source description: Perform Transient Analysis for CSRL - X, U.
*	I. Discretization for CSRL - U.
*		Uses: DEF: bitset.h, sparse.h, label.h, runtime.h
*			LIB: bitset.c, sparse.c, label.c, runtime.c
*		Remarks: Discretization for CSRL - U.
*			SLOW setting FTSKN
*		For discussion see:
*		  1. B.R. Haverkort, L. Cloth, H.Hermanns, J.-P. Katoen,
*		     C. Baier. Model Checking Performability Properties.
*		     DSN-2002, IEEE CS Press, pp. 103-112, 2002.
*		  2. M. Khattri and R. Pulungan. Model Checking Markov
*		     Reward Models with Impulse Rewards. M.S. Thesis.
*		     University of Twente. 2004.
*		  3. L. Cloth, J.-P. Katoen, M. Khattri, R. Pulungan,
*		     Model Checking Markov Reward Models with Impulse
*		     Rewards, DSN-PDS-05.
*
*		NOTE: This is too slow though. haven't had time to think
*			of something better.
*	II. Uniformization based algorithm by Bruno Sericola for CSRL - U.
*		Uses: DEF: bitset.h, sparse.h, label.h, runtime.h
*			LIB: bitset.c, sparse.c, label.c, runtime.c
*		Remarks: Uniformization based algorithm by Bruno Sericola
*			for CSRL - U.
*		For discussion see:
*		  1. Bruno Sericola. Occupation Times in Markov Processes.
*                    Stochastic Models: 16(5), 2000. pp. 479-510.
*		  2. B.R. Haverkort, L. Cloth, H.Hermanns, J.-P. Katoen,
*		     C. Baier. Model Checking Performability Properties.
*		     DSN-2002, IEEE CS Press, pp. 103-112, 2002.
*		  Only valid for MRMs without impulse rewards.
*		NOTE: INCOMPLETE - Needs to be implemented
*	III. Uniformization based algorithm by Qureshi&Sanders for CSRL
*		- U.
*		Uses: DEF: bitset.h, sparse.h, label.h, runtime.h
*			LIB: bitset.c, sparse.c, label.c, runtime.c
*		Remarks: Uniformization based algorithm by
*			Qureshi&Sanders for CSRL - U.
*
*		For discussion see:
*                    1. M.A. Qureshi and W.H.Sanders. Reward Model
*			Solution Methods with Impulse and Rate Rewards:
*			An Algorithm and Numerical Results.
*			Performance Evaluation, 20(4), 1994.
*                    2. M.A. Qureshi and W.H. Sanders. A new
*			methodology for calculating distributions of
*			reward accumulated during a finite interval.
*			Fault-Tolerant Computing Symposium, IEEE CS Press,
*			pp. 116-125, 1996.
*
*		FOR OMEGA SEE
*                    3. M.C. Diniz, E. de Souza e Silva and H.R. Gail.
*			Calculating the distribution of a linear
*			combination of uniform of order statistics.
*			INFORMS Journal on Computing: 14(2), pp. 124-131,
*			2002.
*
*		FOR Model Checking Markov Reward Models with Impulse
*		Rewards see:
*                    4. M. Khattri, R. Pulungan, Model Checking Markov
*			Reward Models with Impulse Rewards, M.Sc. Thesis,
*			UT-FMG-2004.
*                    5. L. Cloth, J.-P. Katoen, M. Khattri, R. Pulungan,
*			Model Checking Markov Reward Models with Impulse
*			Rewards, DSN-PDS-05.
*		NOTE: THIS DOES NOT SEEM TO WORK CORRECTLY OR MAY BE IS
*			NUMERICALLY UNSTABLE
*/

#ifndef TRANSIENT_CTMRM_H
#define TRANSIENT_CTMRM_H

#include "bitset.h"

/**
* Solve bounded until formula with rewards.
* Does optimization: excludes states from phi_and_not_psi
* from which you always go to bad i.e. not Phi and not Psi states
* Also removes Phi & Psi states
* NOTE: The lower bound for time and reward is zero that is
* why the optimization is possible!
* @param phi the phi states.
* @param psi the psi states.
* @param supi the upper time bound
* @param supj the upper reward bound
* @param ppResultError  the pointer to the array of doubles. This array will
*                       store the errors
*			for the resulting probabilities. (This is the return variable)
* @return the result of the until formula.
*/
extern double * ctmrm_bounded_until(const bitset * phi, const bitset * psi,
                double supi, double supj, double ** ppResultError);

/**
* Solve bounded until formula with rewards with lumping.
* NOTE: The lower bound for time and reward is zero
* @param phi the phi states.
* @param psi the psi states.
* @param supi the upper time bound
* @param supj the upper reward bound
* @param ppResultError  the pointer to the array of doubles. This array will
*                       store the errors
*			for the resulting probabilities. (This is the return variable)
* @return the result of the until formula for all states.
*/
extern double *ctmrm_bounded_until_lumping(const bitset *phi, const bitset *psi,
                double supi, double supj, double ** ppResultError);

/**
* Solve next_rewards formula.
* @param phi satisfaction relation for phi formula.
* @param subi rime bound
* @param supi time bound
* @param subj reward bound
* @param supj reward bound
* @return result of the next formula for all states.
* NOTE: ********without impulse rewards at present***********
*		  ********to be tested***********
*/
extern double * ctmrm_interval_next(const bitset *phi, double subi, double supi,
                double subj, double supj);

/* Just for the test cases, should not be used otherwise outside the transient_ctmrm.c file */
extern double omega(const int * tmp_k, double r_dash, int pivot, int ndsr,
                const double * dsr);

#endif
