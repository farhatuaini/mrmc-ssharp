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
*	Authors: Maneesh Khattri, Ivan Zapreev
*
*       Copyright (C) (this file) Saarland University, 2007
*       Authors: Moritz Hahn
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
*       Our contacts (Saarland University).
*               Sven Johr, Dependable Systems and Software
*               Stuhlsatzenhausweg 45
*               66123 Saarbruecken
*               Phone: +49 681 302 5607
*               Fax: +49 681 302 5636
*               Location:   Bldg. E1 3, Room 536
*               E-Mail: johr@cs.uni-sb.de
*
*       Source description: Helper functions to handle MDP data structures
*/
#ifndef MDP_SPARSE_H
#define MDP_SPARSE_H

#include "mdp_labelset.h"

/* datatypes and functions for MDPs with internal non-determinism
 */

/* taken from PRISM */

/* nondeterministic (mdp) sparse matrix */

typedef struct NDSparseMatrix NDSparseMatrix;

/**
* The NDsparseMatrix structure.
* @member n		: The number of states.
* @member nc		: The number of choices.
* @member nnz		: The number of non zeros in the matrix.
*
* @member k		: The maximum number of choices in a state.
* @member use_counts	: Indicates whether counts are stored.
* @member mem		: The memory used for storing.
*
* @member non_zeros	: The pointer to non zero values of the matrix.
* @member cols		: The pointer to the columns of the matrix.
* @member row_counts	:
* @member choice_counts	:
* @member choice_action	:
* @member label		: The pointer to the transition labels.
* @member rate		:
* @member uniform       : true iff CTMDP is uniform
*/
struct NDSparseMatrix
{
	int n;				/* num states */
	int nc;				/* num choices */
	int nnz;			/* num non zeros */
/*	int nm;				num matrices (upper bound on max num choices in a state) */
	int k;				/* max num choices in a state */
        BOOL use_counts;                /*store counts? (as opposed to starts)*/
	double mem;			/* memory used */

	double *non_zeros;
	unsigned int *cols;
	unsigned char *row_counts;
	unsigned char *choice_counts;
	long *choice_action;
        mdp_labelset *label;
        double rate;
        BOOL uniform;
};

/**
* This function frees the sparse matrix.
* @param sparse the matrix to be freed
*/
extern void NDSparseMatrix_free(NDSparseMatrix *);

#endif
