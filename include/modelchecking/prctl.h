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
*	Authors: Maneesh Khattri
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
*	Source description: Solve formulas in PRCTL.
*/

#ifndef PRCTL_H
#define PRCTL_H

#include "bitset.h"

/*****************************************************************************
name		: ef
role		: solve E formula in DMRM
@param		: int steps: n.
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of E formula in DMRM.
remark		:
******************************************************************************/
extern double * ef(int, const bitset *);

/*****************************************************************************
name		: cf
role		: solve C formula in DMRM
@param		: int steps: n.
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of N-R bounded C formula in DMRM.
remark		:
******************************************************************************/
extern double * cf(int, const bitset *);

/*****************************************************************************
name		: yf
role		: solve Y formula in DMRM
@param		: int steps: n.
@param		: bitset *phi: satisfaction relation for phi formula.
@return		: double *: result of N-R bounded C formula in DMRM.
remark		:
******************************************************************************/
extern double * yf(int, const bitset *);

#endif
