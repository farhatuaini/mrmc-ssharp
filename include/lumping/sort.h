/**
*
*	MRMC is a model checker for discrete-time and continuous-time Markov
*	reward models. It supports reward extensions of PCTL and CSL (PRCTL
*	and CSRL), and allows for the automated verification of properties
*	concerning long-run and instantaneous rewards as well as cumulative
*	rewards.
*
* 	Copyright (C) The University of Twente, 2004-2008.
* 	Copyright (C) RWTH Aachen, 2008-2009.
*	Authors: David N. Jansen
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
*	Source description: Sort partition structure for lumping.
*	Uses: DEF: partition.h
*/

#ifndef SORT_H
#define SORT_H

#include "partition.h"

/**
* The function refines a block according to the floating point numbers in key.
* It splits the block into subsets; each subset consists of states that have
* (almost) the same key value.
* Parameters: P = partition
*	B = block to be refined;
*	key = array of doubles (key[s] = value for state s).
*	set_b = TRUE if the array P->states[...].b has to be set to the correct
*		value;
*	largest = largest block found until now. may be NULL.
* Result: the largest block found (either the parameter largest or a larger
*	block). NULL if some error has happened.
*/
extern err_state sort_and_split_block(partition * P, /*@dependent@*/ block * B,
		/*@observer@*/ const double * key, BOOL set_b)
		/*@requires notnull P->blocks, P->first_Sp, P->pos@*/
		/*@modifies P->num_blocks, P->ib[], B->u.next_Sp, B->next,
			*P->pos@*/;

#endif
