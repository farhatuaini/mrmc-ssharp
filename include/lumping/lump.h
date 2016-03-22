/**
*	WARNING: Do Not Remove This Section
*
*	$LastChangedRevision: 415 $
*	$LastChangedDate: 2010-12-18 17:21:05 +0100 (Sa, 18. Dez 2010) $
*	$LastChangedBy: davidjansen $
*
*	MRMC is a model checker for discrete-time and continuous-time Markov
*	reward models. It supports reward extensions of PCTL and CSL (PRCTL
*	and CSRL), and allows for the automated verification of properties
*	concerning long-run and instantaneous rewards as well as cumulative
*	rewards.
*
*	Copyright (C) The University of Twente, 2004-2008.
*	Copyright (C) RWTH Aachen, 2008-2009.
*       Author: David N. Jansen
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
*	Source description: Lumping of Markov chains.
*	Uses: DEF: lump.h, partition.h, splay.h, sparse.h, bitset.h,
*		label.h, runtime.h, macro.h,
*		LIB: partition.c, splay.c, sparse.c, bitset.c, label.c,
*/


#ifndef LUMP_H
#define LUMP_H

#include "sparse.h"
#include "partition.h"

/**
* This method computes the optimal lumped matrix of the input matrix.
* @param	: partition *P: The initial partition.
* @param	: sparse *Q: The matrix to lump.
* @return	: The lumped matrix.
*/
extern /*@only@*/ /*@null@*/ sparse * lump(partition * P,
                /*@observer@*/ const sparse * Q)
                /*@requires notnull P->blocks, P->pos@*/
                /*@requires isnull P->sum@*/
                /*@ensures isnull P->first_Sp, P->first_PredCl, P->sum,
                        P->pos@*/
                /*@modifies *P@*/;

/**
* This method changes the labelling structure such that it corresponds to the partition.
* @param	: labelling *labellin: the labelling structure.
* @param	: partition *P: the partition.
*/
extern err_state change_labelling(labelling * labellin,
                /*@observer@*/ const partition * P)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl, P->sum,
                        P->pos@*/
                /*@modifies *labellin@*/;

/**
* This method changes the state rewards structure such that it
* corresponds to the partition. The actual change is done only
* in the case when the provided p_old_rew pointer is not NULL.
* NOTE:		We do not free the old_rewards structure because
*		it should be done on the outer level.
* @param	: double * p_old_rew: the the rewards structure.
* @param	: partition *P: the partition.
* @return	: The rewards structure after lumping or NULL if the original
*		  rewards structure is NULL
*/
extern /*@only@*/ /*@null@*/ double * change_state_rewards(
                /*@observer@*/ const double * p_old_rew,
                /*@observer@*/ const partition * P)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl, P->sum,
                        P->pos@*/
                /*@modifies nothing@*/;

/**
* Unlump a probability vector with respect to the partition.
* @param	: partition *P: the partition of the original state space.
* @param	: int size: the number of states in the original state space.
* @param	: double *result_lumped: the probability vector of the lumped state space.
* @return	: the probability vector of the original state space.
*/
extern /*@only@*/ /*@null@*/ double * unlump_vector(
                /*@observer@*/ const partition * P, state_count size,
                /*@observer@*/ const double * result_lumped)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl, P->sum,
                        P->pos@*/
                /*@modifies nothing@*/;

/**
* Lump a bitset with respect to the partition.
* @param	: partition *P: the partition of the original state space.
* @param	: bitset *phi: the original bitset.
* @return	: the lumped bitset.
*/
extern /*@only@*/ /*@null@*/ bitset * lump_bitset(
                /*@observer@*/ const partition * P,
                /*@observer@*/ const bitset * phi)
                /*@requires notnull P->blocks@*/
                /*@requires isnull P->first_Sp, P->first_PredCl, P->sum,
                        P->pos@*/
                /*@modifies nothing@*/;

#endif
