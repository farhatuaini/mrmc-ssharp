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
*	Source description: Read impulse rewards (.rewi) file.
*	Uses: DEF: sparse.h
*		LIB: sparse.c
*		Definition of read_impulse_rewards - read_impulse_rewards.h
*/

# include "read_impulse_rewards.h"

#include <stdlib.h>
#include <string.h>

static const char STATES[] = "STATES";

static int * make_first_pass_impulse(const char * filename, int size)
{
        FILE *p;
        char  s[1024],transitions[11];
	int row, col, nnz, *ncolse = NULL;
	double val;
	p = fopen(filename, "r");
	if(p==NULL) return NULL;
	if(fgets( s, 1024, p )!=NULL)
	{
		sscanf( s, "%s%d", transitions, &nnz);
		/* If there is a number of states present */
		if(strncmp(transitions, STATES, strlen(STATES)) == 0)
		{
			if(fgets( s, 1024, p )!=NULL){
				sscanf( s, "%s%d", transitions, &nnz);
			}
		}

                ncolse = (int *) calloc((size_t) size, sizeof(int));
		while (NULL != fgets( s, 1024, p ))
		{
			sscanf( s, "%d%d%lf", &row, &col, &val );
			if(row!=col)
				++ncolse[row-1];
		}
	}
	(void)fclose(p);
	return ncolse;
}

/*****************************************************************************
name		: read_impulse_rewards
role		: reads a .rewi file. puts the result in a sparse matrix (sparse.h).
@param		: char *filename: input .rewi file's name.
@return		: sparse *: returns a pointer to a sparse matrix.
remark		:
******************************************************************************/
sparse * read_impulse_rewards(const char * filename, int size)
{
	FILE *p;
	char  s[1024],transitions[11];
	int row, col, nnz, *ncolse;
	double val = 0.0;
	sparse *sp = NULL;

	ncolse = make_first_pass_impulse(filename, size);

	p = fopen(filename, "r");
	if(p==NULL) return NULL;
	if(fgets( s, 1024, p )!=NULL)
	{
		sscanf( s, "%s%d", transitions, &nnz);
		/* If there is a number of states present */
		if(strncmp(transitions, STATES, strlen(STATES)) == 0)
		{
			if(fgets( s, 1024, p )!=NULL){
				sscanf( s, "%s%d", transitions, &nnz);
			}
		}

		sp = allocate_sparse_matrix_ncolse(size, size, ncolse);
                if ( NULL == sp ) {
                        err_msg_2(err_CALLBY, "read_impulse_rewards(\"%s\",%d)",
                                filename, size, (free(ncolse), NULL));
                }
		while (NULL != fgets( s, 1024, p ))
		{
			sscanf( s, "%d%d%lf", &row, &col, &val );
                        if ( err_state_iserror(set_mtx_val_ncolse(sp, row - 1,
                                        col - 1, val)) )
                        {
                                err_msg_2(err_CALLBY, "read_impulse_rewards("
                                        "\"%s\",%d)", filename, size,
                                        (free_sparse_ncolse(sp), free(ncolse),
                                        NULL));
                        }
		}
	}

	(void)fclose(p);
	free(ncolse);

	return sp;
}
