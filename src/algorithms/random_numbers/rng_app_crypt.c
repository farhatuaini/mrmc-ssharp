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
*	Authors: Ivan Zapreev, Christina Jansen
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
*	Source description:
*		A C library for generating random numbers.
*		The method used for generating random numbers is a combined
*		linear congruential one and returns uniformly distributed
*		random numbers between 0 and 1.
*
*		For more information see: "Applied Cryptography",
*		Bruce Schneier, John Wiley & Sons, 1996
*/

#include "rng_app_crypt.h"

#include "macro.h"

#include <stdlib.h>
#include <time.h>

/* Macro for calculating parameters for choosing a uniformly distributed random
 number,
 approach taken from Applied Cryptography", Bruce Schneier */
#define MODMULT(a, b, c, m, s) \
                do { \
                        long m_q = (s) / (a); \
                        (s) = (b) * ((s)-(a)*m_q) - (c)*m_q; \
                        if ( (s) < 0 ) \
                                (s) += (m); \
                } while ( FALSE )

static long s1, s2;

/**
* Choose uniformly distributed random number between 0 and 1, approach taken
* from
* "Applied Cryptography" by Bruce Schneier
* @param pRandGenerator the rand number generator structure, needed for GSL.
* @return a random value between 0 and 1
*/
double generateRandUnifAC(void * UNUSED(pRandGenerator)){
        long z;

        MODMULT(53668, 40014, 12211, 2147483563L, s1);
        MODMULT(52774, 40692, 3791, 2147483399L, s2);
	z = s1 - s2;
	if (z < 1){
		z += 2147483562;
	}
	return z * 4.656613e-10;
}

/**
* This function seeds the random number generator.
* @param pRandGenerator the rand number generator structure, needed for GSL.
* @param _s1 the value the stream of random numbers is started with
*/
void generateSeedAC(void * UNUSED(pRandGenerator), unsigned long _s1){
	srand( (unsigned)time( NULL) );
	s1 = _s1;
	s2 = rand();
}
