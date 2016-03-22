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
*       Copyright (C) (this file) Saarland University, 2007-2009
*       Author: Ernst Moritz Hahn
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
*       Our contacts (Saarland University).
*               Ernst Moritz Hahn, Dependable Systems and Software
*               Campus Saarbruecken
*               66123 Saarbruecken
*               Phone: +49 681 302 5607
*               Fax: +49 681 302 5636
*               Location:   Bldg. E1 3, Room 536
*               E-Mail: emh@cs.uni-sb.de
*
*       Source description: Reads an MDP
*/

#include "read_mdpi_file.h"

#include "token.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Defining this variable as a pointer means that gcc reserves extra space to
   keep that pointer; what we actually want, however, is only a string constant,
   not a constant pointer to a string constant. David N. Jansen. */
static const char END_DECL[] = "#END";

/**
* Checks for uniformity and sets flag in CTMDP.
*
* @param sparse CTMDP to check uniformity of
*/
static void check_uniformity(NDSparseMatrix *sparse)
{
	double E = 0.0;
	double E_primed;
	double delta;
	BOOL result = TRUE;
	unsigned i;
	int state_nr;
	unsigned choice_nr;
	unsigned *row_starts = (unsigned *) sparse->row_counts;
	unsigned *choice_starts = (unsigned *) sparse->choice_counts;
	double *non_zeros = sparse->non_zeros;
	const double epsilon = 1.0E-10;

	/* calculate reference rate */
	/* Take the first non-zero rate we find.
	 * As soon as we have a situation where this is not the case, we
	 * leave the loops.
	 */
	for (state_nr = 0; state_nr < sparse->n; state_nr++) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			/* Add up all outgoing rates of the distribution */
			unsigned i_start = choice_starts[choice_nr];
			unsigned i_end = choice_starts[choice_nr + 1];
			E = 0.0;
			for (i = i_start; i < i_end; i++) {
				E = E + non_zeros[i];
			}
			if (E > 0.0) {
				break;
			}
		}
		if (E > 0.0) {
			break;
		}
	}

	for (state_nr = 0; state_nr < sparse->n; state_nr++) {
		unsigned state_start = row_starts[state_nr];
		unsigned state_end = row_starts[state_nr + 1];
		for (choice_nr = state_start; choice_nr < state_end; choice_nr++) {
			/* Add up all outgoing rates of the distribution */
			unsigned i_start = choice_starts[choice_nr];
			unsigned i_end = choice_starts[choice_nr + 1];
			E_primed = 0.0;
			for (i = i_start; i < i_end; i++) {
				E_primed = E_primed + non_zeros[i];
			}
			/* if non-zero rate differs more than epsilon
			 * from reference rate this is incorrect */
			delta = fabs(E - E_primed);
			if ((E_primed > 0) && (delta > epsilon)) {
				result = FALSE;
				break;
			}
		}
		if (FALSE == result) {
			break;
		}
	}
	sparse->uniform = result;
}

/**
* Checks whether a line of a MDP file is of correct form.
*
* @param line the line itself
* @param line_nr number of the line
* @param filename name of file the line is of
*/
static BOOL
check_line_NDSparseMatrix
(const char *line, const NDSparseMatrix *mdp, int line_nr, const char *filename)
{
	BOOL error = 0;
	int state;
	char label[MAX_LINE_LENGTH];
	int label_nr;
	double rate;
	char star[MAX_LINE_LENGTH];

	int num_states = mdp->n;
	if (NULL == line) {
		/* this should not happen */
		fprintf(stderr, "Internal error.\n");
		error = TRUE;
	}
	if (!error &&
	((strcmp(line, "") == 0) || (strcmp(line, "\n") == 0))) {
		fprintf(stderr, "Line %d of file \"%s\" is empty.\n", line_nr, filename);
		error = TRUE;
	}
	if (!error && (line[0] != '*')) {
		if (2 != sscanf(line, "%d%s", &state, label)) {
			fprintf(stderr, "Invalid line %d of file \"%s\"\n",
			line_nr, filename);
			error = TRUE;
		}else if (state > num_states) {
			fprintf(stderr, "In line %d of file \"%s\" 'from' state exceeds "
					"number of states in model.\n", line_nr, filename);
			error = TRUE;
		}else if (state < 1) {
			fprintf(stderr, "In line %d of file \"%s\" non-positive 'from'-state "
			"is given.\n", line_nr, filename);
			error = TRUE;
		}else {
			label_nr = mdp_labelset_get_label_no(mdp->label, label);
			if (-1 == label_nr) {
				fprintf(stderr, "In line %d of file \"%s\" label \"%s\" not "
				"specified before.\n", line_nr, filename, label);
				error = TRUE;
			}
		}
	}
	if (!error && (line[0] == '*')) {
		if ((3 != sscanf(line, "%s%d%lf", star, &state, &rate))
		|| (strcmp(star, "*") != 0)) {
			fprintf(stderr, "Invalid line %d of file \"%s\"\n",
			line_nr, filename);
			error = TRUE;
		}else if (state > num_states) {
			fprintf(stderr, "In line %d of file \"%s\" 'to' state exceeds "
			"number of states in model.\n", line_nr, filename);
			error = TRUE;
		}else if (state < 1) {
			fprintf(stderr, "In line %d of file \"%s\" non-positive 'to'-state "
			"is given.\n", line_nr, filename);
			error = TRUE;
		}
	}

	return error;
}

/**
* Creates a new CTMDPI with given number of states and transition labelset.
*
* @param num_states number of states of new CTMDPI
* @param labels transition labels
* @return new CTMDPI
*/
static
NDSparseMatrix * NDSparseMatrix_new(unsigned num_states, mdp_labelset *labels)
{
        unsigned * row_starts = (unsigned *) calloc((size_t) (num_states + 1),
                        sizeof(unsigned));
        NDSparseMatrix*result = (NDSparseMatrix*)malloc(sizeof(NDSparseMatrix));
	result->label = labels;
	result->n = num_states;
	result->use_counts = 0;
	result->row_counts = (unsigned char *) row_starts;
	result->non_zeros = NULL;
	result->cols = NULL;
	result->choice_counts = NULL;
	result->choice_action = NULL;

	return result;
}

/**
* Read number of states from CTMDPI file.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p CTMDPI file
* @param filename filename of CTMDPI file
* @param num_states will store number of states here
*/
static void read_num_states(int *line_no, BOOL *error, FILE *p,
			    const char *filename, int *num_states) {
	char s[MAX_LINE_LENGTH];
	char states[MAX_LINE_LENGTH];

	if (!*error) {
		if (NULL == fgets(s, MAX_LINE_LENGTH, p)) {
			fprintf(stderr, "Reading line %d of file \"%s\" failed.\n",
			  *line_no, filename);
			*error = TRUE;
		}
	}
	if (!*error) {
		if ((2 != sscanf(s, "%s%d", states, num_states))
		  || (0 != strcmp(states, "STATES"))) {
			fprintf(stderr, "Line %d of file \"%s\" should be\n"
			  "STATES <number-of-states>\n", *line_no, filename);
			*error = TRUE;
		} else if (*num_states <= 0) {
			fprintf(stderr, "Non-positive number of states not allowed\n");
			*error = TRUE;
		}
	}
	(*line_no)++;
}

/**
* Read number of transition labels from CTMDPI file.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p CTMDPI file
* @param filename filename of CTMDPI file
* @param num_labels number of transition labels will be stored here
*/
static void read_num_labels(int *line_no, BOOL *error, FILE *p,
			   const char *filename, int *num_labels) {
	char s[MAX_LINE_LENGTH];
	char token[MAX_LINE_LENGTH];
	int c;
	int key;

	if (!*error) {
		if (NULL == fgets(s, MAX_LINE_LENGTH, p)) {
			fprintf(stderr, "Reading line %d of file \"%s\" failed.\n",
			  *line_no, filename);
			*error = TRUE;
		}
	}
	(*line_no)++;

	if (!*error) {
		*num_labels = 0;
		while (NULL != fgets(s, MAX_LINE_LENGTH, p)) {
			c = 0;
			while ((key = get_next_token(s,&c,token)) != EOL
			  && (0 != strcmp(token,END_DECL))){
				++(*num_labels);
			}
			if (0 == strcmp(token,END_DECL)) {
				break;
			}
		}
	}
}

/**
* Read transition labels used in CTMDPI.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p CTMDPI file (rewind before calling function!)
* @param labels labelset structure to be created
*/
static void read_labels
(int *line_no, BOOL *error, FILE *p, const int num_labels,
 mdp_labelset **labels) {
	int c;
	int key;
	char token[MAX_LINE_LENGTH];
	char s[MAX_LINE_LENGTH];

	if (!*error) {
		*labels = mdp_labelset_new(num_labels);
		if (NULL == fgets(s, MAX_LINE_LENGTH, p)) {
			fprintf(stderr, "CTMDP file has incorrect header\n");
			*error = TRUE;
		} else if (NULL == fgets(s, MAX_LINE_LENGTH, p)) {
			fprintf(stderr, "CTMDP file has incorrect header\n");
			*error = TRUE;
		} else {
			while (NULL != fgets(s, MAX_LINE_LENGTH, p)) {
				c = 0;
				while ((key = get_next_token(s,&c,token)) != EOL
				  && (0 != strcmp(token,END_DECL))){
					mdp_labelset_add(*labels, token);
				}
				if (0 == strcmp(token, END_DECL)) {
					break;
				}
				(*line_no)++;
			}
			(*line_no)++;
		}
	}
}

/**
* Allocate memory needed for CTMDPI data structures.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p CTMDPI file
* @param ctmdpi CTMDPI for which allocations will be done
*/
static void reserve_transition_memory
(int *line_no, BOOL *error, FILE *p, const char *filename,
NDSparseMatrix *ctmdpi) {
	char s[MAX_LINE_LENGTH];
	int last_from = 0;
	int from;
	char label[MAX_LINE_LENGTH];
	unsigned num_choice = 0;
	unsigned num_non_zeros = 0;
	unsigned *row_starts;

	if (!*error) {
		row_starts = (unsigned *) ctmdpi->row_counts;
	}

	while (!*error && (NULL != fgets(s, MAX_LINE_LENGTH, p))) {
		*error = check_line_NDSparseMatrix(s, ctmdpi, *line_no, filename);
		if (!*error) {
			if (s[0] != '*') {
				/* state distribution starts in line read */
				if (2 != sscanf(s, "%d%s", &from, label)) {
					*error = TRUE;
				} else {
					/* we just had a line "state" -> "label", so increment the
					* number of internal non-deterministic mappings by one */
					from--;
					/* state transitions must be given in ascending order */
					if (from < last_from) {
                                                fprintf(stderr, "State "
                                                        "transitions must be "
                                                        "given in ascending"
							" order\n");
						*error = TRUE;
					} else {
						num_choice++;
						if (from == last_from) {
						} else {
							/* new state rows start after last one */
							last_from = from;
						}
					}
				}
			} else {
				num_non_zeros++;
			}
			++(*line_no);
		} else {
			*error = TRUE;
		}
	}

	/* now allocate the memory needed */
	if (!*error) {
                unsigned * choice_starts = (unsigned *)
                        calloc((size_t) (num_choice + 1), sizeof(unsigned));
                double * non_zeros =
                        (double *) malloc(num_non_zeros * sizeof(double));
                unsigned * cols =
                        (unsigned *) malloc(num_non_zeros * sizeof(unsigned));
		ctmdpi->choice_counts = (unsigned char *) choice_starts;
		ctmdpi->non_zeros = non_zeros;
		ctmdpi->cols = cols;
	}
}

/**
* Reads the transitions of a CTMDPI.
*
* @param line_no line number in file (should be 1 before)
* @param error error flag
* @param p CTMDPI file
* @param ctmdpi CTMDPI the transitions shall be added to
*/
static void read_transitions
(int *line_no, BOOL *error, FILE *p, NDSparseMatrix *ctmdpi) {
	unsigned choice_index = 0;
	unsigned choice_size = 0;
	unsigned nz_index = 0;
	int from;
	int last_from = 0;
	char s[MAX_LINE_LENGTH];
	int c;
	int key;
	char token[MAX_LINE_LENGTH];
	int to;
	char label[MAX_LINE_LENGTH];

	if (!*error) {
		double *non_zeros = ctmdpi->non_zeros;
		unsigned *cols = ctmdpi->cols;
		unsigned *choice_starts = (unsigned *) ctmdpi->choice_counts;
		unsigned *row_starts = (unsigned *) ctmdpi->row_counts;

		/* ignore first two lines as we already have their values */
		if (NULL == fgets(s, MAX_LINE_LENGTH, p)) {
			*error = TRUE;
		} else if (NULL == fgets(s, MAX_LINE_LENGTH, p)) {
			*error = TRUE;
		}
		while ((!*error) && (NULL != fgets(s, MAX_LINE_LENGTH, p))) {
			c = 0;
			while ((key = get_next_token(s,&c,token)) != EOL
			  && (0 != strcmp(token,END_DECL))){
				/* don't do anything */
			}
			if (0 == strcmp(token,END_DECL)) {
				break;
			}
			(*line_no)++;
		}

		while (NULL != fgets(s, MAX_LINE_LENGTH, p)) {
			if (s[0] == '*') {
				double rate;
				char star[MAX_LINE_LENGTH];
				sscanf(s, "%s%d%lf", star, &to, &rate);
				to--;
				/* just read a line belonging to internal non-determinism map */
				non_zeros[nz_index] = rate;
				cols[nz_index] = to;
				nz_index++;
				choice_size++;
			} else {
				/* just read a line "state" -> "label" */
				sscanf(s, "%d%s", &from, label);
				from--;
				if (from == last_from) {
					row_starts[from + 1]++;
				} else {
					/* new state rows start after last one */
					for (; last_from < from; last_from++) {
						row_starts[last_from + 2] = row_starts[last_from + 1];
					}
					row_starts[from + 1]++;
					last_from = from;
				}
				if (choice_index > 0) {
					choice_starts[choice_index] =
					choice_starts[choice_index - 1] + choice_size;
				}
				choice_size = 0;
				choice_index++;;
			}
		}
		choice_starts[choice_index] = choice_starts[choice_index - 1] + choice_size;
		for (; from+1 < ctmdpi->n; from++) {
			row_starts[from+2] = row_starts[from+1];
		}
	}
}

/**
* Reads CTMDPI from file @a filename.
*
* @param filename file to read CTMDPI from
* @return CTMDPI read from file
*/
NDSparseMatrix *read_NDSparseMatrix_file(const char *filename)
{
	BOOL error = FALSE;
	int line_no;
	NDSparseMatrix *result = NULL;
	mdp_labelset *labels = NULL;
	FILE *p = NULL;
	int num_states = 0;
	int num_labels = 0;

	/* To minimize the memory usage we read the file in several passes,
	* which allows us to reserve the amount of memory needed before
	* inserting in the data structure as we get this information from
	* the former pass.
	* While reading the file, the size fields of the data structure are
	* misused as indices (where to place sub-structure next), so we
	* don't need to store this information in another place. At the end
        * it is ensured that these size fields fulfil their original purpose.
	*/

	if (NULL == filename) {
		fprintf(stderr, "Called with filename == NULL\n");
		error = TRUE;
	}

	/* first pass on file: get number of states and number of labels. */
	if (!error) {
		p = fopen(filename, "r");
		if (NULL == p) {
			fprintf(stderr, "Could not open file \"%s\"\n", filename);
			error = TRUE;
		}
	}

	line_no = 1;
	read_num_states(&line_no, &error, p, filename, &num_states);
	read_num_labels(&line_no, &error, p, filename, &num_labels);

	if (NULL != p) {
		rewind(p);
	}

	/* second pass: read set of labels, reserve transition memory */
	read_labels(&line_no, &error, p, num_labels, &labels);
	if (!error) {
                result = NDSparseMatrix_new((unsigned) num_states, labels);
	}
	reserve_transition_memory(&line_no, &error, p, filename, result);
	if (!error) {
		rewind(p);
	}

	/* third pass: read model transitions */
	read_transitions(&line_no, &error, p, result);

	if (NULL != p) {
		fclose(p);
	}
	if (error) {
		/* free the halfly-complete MDP structure if an error has occured */
		NDSparseMatrix_free(result);
		result = NULL;
	} else {
	  check_uniformity(result);
	}

	return result;
}
