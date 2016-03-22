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
*       Our contacts (Saarland University)
*               Sven Johr, Dependable Systems and Software
*               Stuhlsatzenhausweg 45
*               66123 Saarbruecken
*               Phone: +49 681 302 5607
*               Fax: +49 681 302 5636
*               Location:   Bldg. E1 3, Room 536
*               E-Mail: johr@cs.uni-sb.de
*
*       Source description: Transition labels for MDPs
*/

#include "mdp_labelset.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
* This function creates a new labelset for an MDP.
* @param size the number of labels to be stored.
*/
mdp_labelset *mdp_labelset_new(int size)
{
	mdp_labelset *result;

        result = (mdp_labelset *) malloc(sizeof(mdp_labelset));
	result->reserved = size;
	result->next = 0;
        result->label = (char **) calloc((size_t) size, sizeof(char *));

	return result;
}

/**
* Frees a MDP transition labelling structure. It should also be
* possible to free partly build structure without crashing if
* possible.
*
* @param labels the transition labelling structure to be freed.
*/
void mdp_labelset_free(mdp_labelset *labels)
{
	int label_no;

	if (NULL != labels) {
		if (NULL != labels->label) {
			for (label_no = 0; label_no < labels->next; label_no++) {
				if (NULL != labels->label[label_no]) {
					free(labels->label[label_no]);
				}
			}
			free(labels->label);
		}
		free(labels);
	}
}

/**
* Adds a label to a set of MDP transition labels. The memory of the
* labelset must have been reserved in advance, the function will
* return 0 if a label shall be added but there is no space left for
* it.
*
* @param labels labelset the transition is to be added to
* @param name label to be added to the labelset.
* @return 1 for success, 0 else
*/
BOOL mdp_labelset_add(mdp_labelset * labels, const char * name)
{
        char * new_name;

	if (labels->reserved <= labels->next) {
                return TRUE;
	}

        new_name = (char *) malloc(strlen(name) + 1);
        if ( NULL == new_name ) {
                fprintf(stderr, "mdp_labelset_add: out of memory error\n");
                return FALSE;
        }
        strcpy(new_name, name);
        labels->label[labels->next] = new_name;
	labels->next++;

        return TRUE;
}

/**
* Returns number of a certain label. That is, the function maps the
* label in its textual form to the index it has in the internal
* representation of the labelset.
*
* @param labels labelset to look up label
* @param name label to get number for
* @return the number of the label in the internal representation
*/
int
mdp_labelset_get_label_no(const mdp_labelset *labels, const char *name)
{
	/**
	* @todo use hashmap or something (If we don't find a library
	* somewhere, we can at last use something learned in "Algorithms
	* and Data Structures"). However, iff the number of labels is
	* rather small, it should not make a difference.
	*/
	int result = -1;
	int i;

	for (i = 0; i < labels->next; i++) {
		if (strcmp(labels->label[i], name) == 0) {
			result = i;
			break;
		}
	}
	return result;
}
