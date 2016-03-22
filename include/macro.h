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
*	Authors: Ivan Zapreev
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
*	Source description: This file contains main macro definitions.
*/

#ifndef MACRO_H
#define MACRO_H

	/**************************************/
	/*******THE SAFETY-IF ASSERTION********/
	/**************************************/
	/* WARNING: SHOULD BE USED ONLY FOR THE PERFORMANCE-CRYTICAL PLACES */
	/* WARNING: MUST NOT BE USED TO ASSERT ON BOOLEAN EXPRESSIONS WITH */
	/*         SIDEEFFECTS! BECAUSE IN CASE SAFETY_FIRST IS UNDEFINED */
	/*         THE BOOLEAN EXPRESSION WILL NOT BE EVALUATED AT ALL!!! */

	/* Should ONLY be defined for the development phase and testing phase. */
	/*#       define SAFETY_FIRST*/

#       ifdef SAFETY_FIRST
		/*@notfunction@*/
#               define IF_SAFETY( condition ) if( condition ){
		/*@notfunction@*/
#               define ELSE_SAFETY } else {
#       else
		/*@notfunction@*/
#               define IF_SAFETY( condition )
		/*@notfunction@*/
#               define ELSE_SAFETY if( FALSE ){
#       endif
	/*@notfunction@*/
#       define ENDIF_SAFETY }

	/**************************************/
	/*******THE SIMPLE BOOLEAN TYPE********/
	/**************************************/
#       ifndef BOOL
		typedef short BOOL;
#       endif

#       ifndef TRUE
#               define TRUE ((BOOL) 1)
#       endif

#       ifndef FALSE
#               define FALSE ((BOOL) 0)
#       endif

	/**************************************/
	/******THE THREE-VALUED LOGIC TYPE*****/
	/**************************************/
	/* For the three valued logic */
#       ifndef TV_LOGIC
                typedef short TV_LOGIC;
#       endif

	/* The YES (TRUE) answer ">0" */
#       ifndef TVL_TT
#               define TVL_TT ((TV_LOGIC) TRUE)
#       endif

	/* The NO (FALSE) answer "=0" */
#       ifndef TVL_FF
#               define TVL_FF ((TV_LOGIC) FALSE)
#       endif

	/* The I do not know answer "<0" */
#       ifndef TVL_NN
#               define TVL_NN ((TV_LOGIC) -1)
#       endif

	/**************************************/
	/******THE MIN OPERATOR DEFINITION*****/
	/**************************************/
#       ifndef MIN
		extern void MIN(/*@sef@*/ int a, /*@sef@*/ int b) /*@*/;
#               define MIN(a,b) ((a) < (b) ? (a) : (b))
#       endif

	/**************************************/
	/******THE BINARY-FLAGS MANAGEMENT*****/
	/**************************************/
	/* These macro definitions are used for working with bit flags */
	extern void set_flag(unsigned flags, unsigned bit) /*@modifies flags@*/;
	extern void clear_flag(unsigned flags, unsigned bit)
			/*@modifies flags@*/;
	extern unsigned test_flag(unsigned flags, unsigned bit) /*@*/;
	/*@-macroassign@*/
#       define set_flag(flags,bit) ((void) ((flags) |= (bit)))
#       define clear_flag(flags,bit) ((void) ((flags) &= ~(bit)))
	/*@=macroassign@*/
#       define test_flag(flags,bit) ((flags) & (bit))

        /**************************************/
        /*MACROS TO GENERATE/SUPPRESS WARNINGS*/
        /**************************************/

#ifndef UNUSED
#       ifdef __GNUC__
#               define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#       else
		/*@notfunction@*/
#               define UNUSED(x) /*@unused@*/ x
#       endif
#endif

	/***********************************************************/
	/*****A GENERAL MACRO DEFINITION FOR THE INLINE FUNCTIONS***/
	/***WE USE IT TO PREFIX THE FUNCTION NAME IN THE .H FILES***/
	/***********************************************************/
#       ifndef INLINE
#               ifdef __GNUC__
#                       define INLINE static inline
#               else
#                       define INLINE inline
#               endif
#       endif
#endif
