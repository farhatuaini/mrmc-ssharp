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
*	Authors: Maneesh Khattri, Ivan S. Zapreev
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
*	Source description: A simple tokenizer.
*	Uses: DEF: token.h
*/

# include "token.h"

#include "macro.h"

#include <string.h>

/* Some required constant arrays definitions with the special symbols */

/* Defining these variables as pointers means that gcc reserves extra space to
   keep the pointers; what we actually want, however, are only string constants,
   not constant pointers to string constants. David N. Jansen. */

/* NOTE: The 3d symbol is '\0' which is present implicitly. */
static const char END_OF_LINE_SYMBOLS[] = "\n\r";
/* NOTE: Here the last symbol i.e. '\0' does not count */
static const char WITHIN_LINE_DELIMITER_SYMBOLS[] = " \t";

extern BOOL isEndOfLineSymbol(/*@sef@*/ char src);
	/* NOTE: We also test for the last (invisible) */
        /* symbol '\0' */
#define isEndOfLineSymbol(src) ('\0' == (src) \
                        || NULL != strchr(END_OF_LINE_SYMBOLS, (src)))

extern BOOL isWithinLineDelimiter(char src);
#define isWithinLineDelimiter(src) \
                        (NULL != strchr(WITHIN_LINE_DELIMITER_SYMBOLS, (src)))

/*****************************************************************************
name		: scan_number
role		: get the next token: scan for a number.
@param		: char *src: the source string to be tokenized.
@param		: int *start: the index from which to commence search for a token in
			 the source string.
@param		: char *token: the token in which the new token is to be put.
@param		: int *idx: the index from which to commence adding in the token.
remark          : If while scanning a non-numeric character is encountered then
                  scan for a
		  string subsequently.
******************************************************************************/
static int scan_number(const char * src, int * start, char * token, int * idx)
{
	int key = NUMBER;

        while ( *idx < MAXTOKENSIZE && isNumeric(src[*start])
                                && ! isWithinLineDelimiter(src[*start])
                                && ! isEndOfLineSymbol(src[*start]) )
        {
		token[(*idx)++] = src[(*start)++];
	}

        if ( ! isNumeric(src[*start]) && ! isWithinLineDelimiter(src[*start])
                                && ! isEndOfLineSymbol(src[*start]) )
        {
		key = ALPHANUMERIC;
                while ( *idx < MAXTOKENSIZE
                                        && ! isWithinLineDelimiter(src[*start])
                                        && ! isEndOfLineSymbol(src[*start]) )
                {
			token[(*idx)++] = src[(*start)++];
		}
	}
	return key;
}

/*****************************************************************************
name		: scan_string
role		: get the next token: scan for a string.
@param		: char *src: the source string to be tokenized.
@param		: int *start: the index from which to commence search for a token in
			 the source string.
@param		: char *token: the token in which the new token is to be put.
@param		: int *idx: the index from which to commence adding in the token.
remark		:
******************************************************************************/
static int scan_string(const char * src, int * start, char * token, int * idx)
{
	int key = ALPHANUMERIC;
	/* printf("ENTER: \\%d\n",(int) '\n'); */
	/* printf("Reading: '"); */
        while ( *idx < MAXTOKENSIZE && ! isWithinLineDelimiter(src[*start])
                                && ! isEndOfLineSymbol(src[*start]) )
        {
		/* printf("\\%d",(int) src[*start]); */
		token[(*idx)++] = src[(*start)++];
	}
	/* printf("'\n"); */
	return key;
}

/*****************************************************************************
name		: get_next_token
role		: get the next token.
@param		: char *src: the source string to be tokenized.
@param		: int *start: the index from which to commence search for a token in
			 the source string.
@param		: char *token: the token in which the new token is to be put.
remark		:
******************************************************************************/
int get_next_token(const char * src, int * start, char * token)
{
	int idx=0, key;
	/* printf("Skipping: '"); */
        while ( isWithinLineDelimiter(src[*start])
                                && ! isEndOfLineSymbol(src[*start]) )
        {
		/* printf("%c",src[*start]); */
		(*start)++;
	}
	/* printf("'\n"); */
	if( isNumeric( src[*start] ) ){
		key = scan_number(src,start,token,&idx);
	}else{
		key = scan_string(src,start,token,&idx);
	}
	token[idx] = '\0';
	if(idx==0){
		key=EOL;
	}

	return key;
}
