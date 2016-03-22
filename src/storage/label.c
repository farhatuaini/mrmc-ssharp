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
*	Source description: Manage Storage for Labelling.
*/

# include "label.h"

#include <string.h>

/*****************************************************************************
name		: get_new_label
role		: get a new labelling structure.
@param		: int n: the number of labels.
@param		: int ns: the number of states.
@return		: labelling *: returns a pointer to a new labelling structure.
remark		:
******************************************************************************/
labelling * get_new_label(int n, int ns)
{
        labelling *new_label = (labelling*)calloc((size_t)1, sizeof(labelling));
	new_label->n = n; new_label->temp_n=0; new_label->ns = ns;
        new_label->label = (char **) calloc((size_t) n, sizeof(char *));
        new_label->b = (bitset **) calloc((size_t) n, sizeof(bitset *));
	return new_label;
}

/*****************************************************************************
name		: add_label
role		: used to add a new label in ascending lexicographic order to the
		  given labelling structure.
@param		: labelling *labellin: the labelling structure.
@param		: char *label: the label to be added.
@return		: int: 0-fail 1-success. Fails when labellin->n number of labels
		  are already added.
remark		: This method also initializes the appropriate bitset with a new
		  bitset.
******************************************************************************/
BOOL add_label(labelling * labellin, const char * newlabel)
{
	int i = labellin->temp_n, c;
        if ( i == labellin->n ) {
                return FALSE;
        }
	for(c=i;c>0;c--)
	{
		if(labellin->label[c-1]&&strcmp(labellin->label[c-1], newlabel)<=0) break;
		labellin->label[c] = labellin->label[c-1];
		labellin->b[c]=labellin->b[c-1];
	}
	labellin->label[c]=(char *)calloc(strlen(newlabel)+1,sizeof(char));
	strcpy(labellin->label[c], newlabel);
	labellin->b[c]=get_new_bitset(labellin->ns);
	++labellin->temp_n;
        return TRUE;
}

/*****************************************************************************
name		: find.
role		: Find the specified label in the given labelling structure.
@param		: labelling *labellin: the labelling structure.
@param		: char *label: the label to be found.
@return		: int: -1-fail other-index of the label.
remark		: Performs binary search.
******************************************************************************/
static
int find( const labelling *labellin, const char *label)
{
	int l=0,u=labellin->n-1, m = 0, res=0;
	while(l<=u)
	{
		m=(l+u)/2;
		res = strcmp(labellin->label[m], label);
		if(res==0) break;
		else if(res>0) u=m-1;
		else if(res<0) l=m+1;
	}
	if(res!=0) return -1;
	return m;
}

/*****************************************************************************
name		: add_label_bitset
role		: set the bitset of the given labelling structure indexed by the
		  given label to a desired value.
@param		: labelling *labellin: the labelling structure.
@param		: char *label: the label whose bitset is to be changed.
@param		: bitset *b: the new bitset.
@return		: int: 0-fail 1-success.
remark		: Fails when the given label cannot be found in the given labelling
		  structure.
******************************************************************************/
BOOL add_label_bitset(labelling * labellin, const char * label, bitset * b)
{
	int res = find(labellin, label);
        if ( -1 == res ) {
                return FALSE;
        }
	labellin->b[res]=b;
        return TRUE;
}

/**
* Set the bit of the given labelling structure indexed by the given
* label and position to 1. Fails when the given label cannot be found
* in the given labelling structure.
* @param labellin the labelling structure.
* @param label the label whose bitset is to be changed.
* @param pos the position of the bit.
*/
void set_label_bit(labelling * labellin, const char * label, int pos)
{
	int res = find(labellin, label);
	if( res == -1 ){
		printf("ERROR: The label '%s' is unknown, check the .lab file.\n", label);
                exit(EXIT_FAILURE);
	}
	set_bit_val(labellin->b[res], pos, BIT_ON);
}

/**
* Get the bitset of the given labelling structure indexed by the
* given label.
* @param labellin: the labelling structure.
* @param label: the label whose bitset is to be changed.
* @return the bitset of states labelled with 'label', if label is
*         not known returns NULL.
*/
bitset * get_label_bitset(const labelling *labellin, const char *label)
{
	int res = find(labellin, label);
	if( res == -1 ){
		printf("WARNING: The given atomic proposition '%s' is unknown, check the .lab file.\n", label);
		return NULL;
	}
	return labellin->b[res];
}

/*****************************************************************************
name		: print_labelling
role		: prints the given labelling structure.
@param		: labelling *labellin: the labelling structure to be printed.
@return		:
remark		:
******************************************************************************/
void print_labelling(const labelling * a)
{
	int i, n;
        if ( NULL != a )
	{
		n = a->n;
		for(i=0;i<n;i++)
		{
			printf("Label[%d]=%s = ",i,a->label[i]);
			if(a->b[i]) print_bitset_states(a->b[i]);
		}
	}
}

/**
* Writes the given labelling structure to a file
* @param labellin the labelling structure to be printed.
* @param fname the name of the output file
* @return 0 if writing failed, 1 otherwise.
* NOTE: This function returns 0 if labelling is not valid, or the file could not
*	be opened for writing. This is also indicated via a printf statement.
*/
BOOL write_lab_file(const labelling * a, const char * fname) {
	int i, j;
	FILE *fp;

	/* No valid labelling */
        if ( NULL == a ) {
                return FALSE;
        }

	printf("Writing labelling file to %s.\n", fname);
	fp = fopen(fname, "w");
        if ( NULL == fp ) {
		printf("Could not open file for writing.\n");
                return FALSE;
	}

	fprintf(fp, "#DECLARATION\n");
	for (i = 0; i < a->n; i++) {
		if( i>0 )
			fprintf(fp, " ");
		fprintf(fp, "%s", a->label[i]);
	}
	fprintf(fp, "\n#END\n");

	/* For every state, we see whether some label is defined. That is, we
	 * check the bit vector to see whether the bit is set for that state.
	 * If so, we about the label for that state. The first time, we also
	 * output the statenumber first.
	 */
	for(i=0;i<a->ns;i++) {
                BOOL printed = FALSE;
		for(j=0;j<a->n;j++) {
			if( a->b[j] && get_bit_val(a->b[j], i) ) {
				if( !printed ) {
                                        printed = TRUE;
					/* The .lab & .tra file start with 1,
					 * but internally we start at 0. */
					fprintf(fp, "%d", i + 1);
				}
				fprintf(fp, " %s", a->label[j]);
			}
		}
		/* Print a final newline, if the state had any labels. */
		if(printed){
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	return TRUE;
}

/*****************************************************************************
name		: free_labelling
role		: frees the given labelling structure.
@param		: labelling *labellin: the labelling structure to be freed.
@return		:
remark		:
******************************************************************************/
void free_labelling(labelling *pLabelling)
{
	int i, n = pLabelling->n;
	for(i=0;i<n;i++)
	{
		if(pLabelling->label[i]) free(pLabelling->label[i]);
		if(pLabelling->b[i]) free_bitset(pLabelling->b[i]);
	}
	if(pLabelling->label) free(pLabelling->label);
	if(pLabelling->b) free(pLabelling->b);
	free(pLabelling);
}
