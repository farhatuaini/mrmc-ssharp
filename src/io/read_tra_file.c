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
*	Authors: Maneesh Khattri, Ivan Zapreev
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
*	Source description: Read transition (.tra) file.
*	Uses: DEF: sparse.h
*		LIB: sparse.c
*		Definition of read_tra_file - read_tra_file.h
*/

#include "read_tra_file.h"

#include <stdlib.h>

/**
* This method does the first pass through the .tra file and computes
* the number of non zero elements in each row. Then it stores these
* values into the ncolse array.
*/
static int * make_first_pass(const char * filename)
{
	FILE *p;
	char  s[1024], states[7], transitions[12];
	int size=0, row=0, col=0, nnz=0, *ncolse=NULL;
	double val=0.0;
	p = fopen(filename, "r");
	if(p==NULL) return NULL;

	if(fgets( s, 1024, p )!=NULL)
	{
		sscanf( s, "%s%d", states, &size);
		if(fgets( s, 1024, p )!=NULL)
		{
			sscanf( s, "%s%d", transitions, &nnz);
                        ncolse = (int *) calloc((size_t) size, sizeof(int));
			while (NULL != fgets( s, 1024, p ))
			{
				sscanf( s, "%d%d%lf", &row, &col, &val );
				if(row != col)
					++ncolse[row-1];
			}
		}
	}
	(void)fclose(p);
	return ncolse;
}

void print_read_mtx(const sparse * sp) {
	FILE *of;
        int nnz, i;

	if(sp == NULL){
		printf("ERROR: The state space matrix seems to be unallocated!\n");
                exit(EXIT_FAILURE);
	}

	of = fopen("matrix_test_read.tra", "w");
        fprintf(of, "STATES %d\n", mtx_rows(sp));

	nnz = 0; i = 0;
        do {
                double val;

                nnz += mtx_next_num(sp, i);
                if ( err_state_iserror(get_mtx_diag_val(sp, i, &val)) )
                        exit(err_macro_3(err_CALLBY, "print_read_mtx(%p[%dx"
                                "%d])", (const void *) sp, mtx_rows(sp),
                                mtx_cols(sp), EXIT_FAILURE));
                if ( 0.0 != val )
			nnz ++;
	}
        while ( ++i < mtx_rows(sp) );

	fprintf(of, "TRANSITIONS %d\n", nnz);
        i = 0;
        do {
                mtx_walk_row_sorted(sp, i+0, j, val) {
                        fprintf(of, "%d %d %.12g\n", i + 1, j + 1, val);
		}
                end_mtx_walk_row_sorted;
	}
        while ( ++i < mtx_rows(sp) );

	(void)fclose(of);
}

/*****************************************************************************
name		: read_tra_file
role		: reads a .tra file. puts the result in a sparse matrix (sparse.h).
@param		: char *filename: input .tra file's name.
@return		: sparse *: returns a pointer to a sparse matrix.
remark		:
******************************************************************************/
sparse * read_tra_file(const char * filename)
{
	FILE *p;
	char  s[1024],states[7],transitions[12];
	int size, row, col, nnz, *ncolse;
	double val = 0.0;
	sparse *sp = NULL;

	ncolse = make_first_pass(filename);

	p = fopen(filename, "r");
	if(p==NULL) return NULL;
	if(fgets( s, 1024, p )!=NULL)
	{
		sscanf( s, "%s%d", states, &size);
		if(fgets( s, 1024, p )!=NULL)
		{
			sscanf( s, "%s%d", transitions, &nnz);
			printf("States=%d, Transitions=%d\n", size, nnz);
			sp = allocate_sparse_matrix_ncolse(size, size, ncolse);
                        if ( NULL == sp ) {
                                err_msg_1(err_CALLBY, "read_tra_file(\"%s\")",
                                        filename, (free(ncolse), NULL));
                        }
			while (NULL != fgets( s, 1024, p ))
			{
				sscanf( s, "%d%d%lf", &row, &col, &val );
                                if ( err_state_iserror(set_mtx_val_ncolse(sp,
                                                        row - 1,col - 1, val)) )
                                        err_msg_1(err_CALLBY, "read_tra_file("
                                                "\"%s\")", filename,
                                                (free_sparse_ncolse(sp),
                                                free(ncolse), NULL));
			}
		}
	}

	(void)fclose(p);
	free(ncolse);

	return sp;
}
