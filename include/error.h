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
*       Author: David N. Jansen
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
*       Source description: Generate error messages.
*/

#ifndef ERROR_H
#define ERROR_H

#include "macro.h"

#include <stdio.h>

        /*************************************/
        /* Special Types for Use with SPLINT */
        /*************************************/

	/* Splint checks whether type usage is consistent. This reduces the
	danger to mix several uses of ``int''.
	The types below should be used in a way that they could be either
	unsigned or signed; for example, the value ``-1'' should be avoided. */

	typedef int state_index; /* index of a state */
	typedef int state_count; /* number of states */
#       define state_index_NONE ((state_index) -1) /* indicates `not found' */
#       define state_index_ERROR ((state_index) -2) /* indicates an error */
#       define state_count_NONE ((state_count) -1)
#       define state_count_ERROR ((state_count) -2)

	/**************************************/
	/*********** Error checking ***********/
	/**************************************/

	typedef enum { err_OK, err_ERROR } err_state;
#       define err_PARAM "\nError: illegal parameter in "
#       define err_MEMORY "\nError: out of memory in "
#       define err_INCONSISTENT "\nError: inconsistent data in "
#       define err_FILE "\nError: file access failed in "
#       define err_CALLBY "called by "

	extern BOOL err_state_iserror(err_state a) /*@*/;
#       define err_state_iserror(a) (err_OK != (a))

	/**
        * These macros handle errors. They print an error message and return.
	* err_macro_<number> is meant for a macro with <number> parameters;
	* the retval is the result of the macro.
	* err_msg_<number> is meant to be used as a statement in a function;
	* the macro includes an explicit ``return'' statement.
	* @param type type of the error, e. g. err_PARAM, err_MEMORY
	* @param func name (with printf-spec of the parameters) of the function
	* @param p1,p2,... parameters of the function. The number of parameters
	*		must fit the printf-spec
	* @param retval value to be returned to the caller
	*/
	/*@notfunction@*/
#       define err_macro_0(type,func,retval) \
			(fprintf(stderr, "%s%s(%u) " func "\n", (type), \
				__FILE__, __LINE__), (retval))

	/*@notfunction@*/
#       define err_macro_1(type,func,p1,retval) \
			(fprintf(stderr, "%s%s(%u) " func "\n", (type), \
				__FILE__, __LINE__, /*@+voidabstract@*/ (p1) \
                                /*@=voidabstract@*/), (retval))

	/*@notfunction@*/
#       define err_macro_2(type,func,p1,p2,retval) \
			(fprintf(stderr, "%s%s(%u) " func "\n", (type), \
				__FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2) /*@=voidabstract@*/), (retval))

	/*@notfunction@*/
#       define err_macro_3(type,func,p1,p2,p3,retval) \
			(fprintf(stderr, "%s%s(%u) " func "\n", (type), \
				__FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2), (p3) /*@=voidabstract@*/), (retval))

	/*@notfunction@*/
#       define err_macro_4(type,func,p1,p2,p3,p4,retval) \
			(fprintf(stderr, "%s%s(%u) " func "\n", (type), \
				__FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2), (p3), (p4) /*@=voidabstract@*/), (retval))

	/*@notfunction@*/
#       define err_macro_5(type,func,p1,p2,p3,p4,p5,retval) \
			(fprintf(stderr, "%s%s(%u) " func "\n", (type), \
				__FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2), (p3), (p4), (p5) /*@=voidabstract@*/), \
                                (retval))

        /*@notfunction@*/
#       define err_macro_6(type,func,p1,p2,p3,p4,p5,p6,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6) /*@=voidabstract@*/), \
                                (retval))

        /*@notfunction@*/
#       define err_macro_7(type,func,p1,p2,p3,p4,p5,p6,p7,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2), (p3), (p4), (p5), (p6), (p7) \
                                /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#       define err_macro_8(type,func,p1,p2,p3,p4,p5,p6,p7,p8,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2), (p3), (p4), (p5), (p6), (p7), (p8) \
                                /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#       define err_macro_9(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2), (p3), (p4), (p5), (p6), (p7), (p8), (p9) \
                                /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#       define err_macro_10(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa) \
                                /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_11(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb) /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_12(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb), (pc) /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_13(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb), (pc), (pd) /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_17(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,pg,ph,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb),(pc),(pd),(pe),(pf),(pg),(ph) \
                                /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_18(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,pg,ph,pi,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb),(pc),(pd),(pe),(pf),(pg),(ph),(pi) \
                                /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_19(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,pg,ph,pi,pj,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb),(pc),(pd),(pe),(pf),(pg),(ph),(pi),(pj) \
                                /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_20(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,pg,ph,pi,pj,pk,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb),(pc),(pd),(pe),(pf),(pg),(ph),(pi),(pj), \
                                (pk) /*@=voidabstract@*/), (retval))

        /*@notfunction@*/
#define err_macro_21(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,pg,ph,pi,pj,pk,pl,retval) \
                        (fprintf(stderr, "%s%s(%u) " func "\n", (type), \
                                __FILE__, __LINE__, /*@+voidabstract@*/ (p1), \
                                (p2),(p3),(p4),(p5),(p6),(p7),(p8),(p9),(pa), \
                                (pb),(pc),(pd),(pe),(pf),(pg),(ph),(pi),(pj), \
                                (pk), (pl) /*@=voidabstract@*/), (retval))

#ifdef S_SPLINT_S
#       define err__start
#       define err__finish
#else
#       define err__start do{
#       define err__finish ;}while(0)
#endif

	/*@notfunction@*/
#       define err_msg_0(type,func,retval) \
                err__start \
                        return err_macro_0((type), func, (retval)) \
                err__finish

	/*@notfunction@*/
#       define err_msg_1(type,func,p1,retval) \
                err__start \
                        return err_macro_1((type), func, (p1), (retval)) \
                err__finish

	/*@notfunction@*/
#       define err_msg_2(type,func,p1,p2,retval) \
                err__start \
                        return err_macro_2((type), func, (p1), (p2), (retval)) \
                err__finish

	/*@notfunction@*/
#       define err_msg_3(type,func,p1,p2,p3,retval) \
                err__start \
                        return err_macro_3((type),func,(p1),(p2),(p3),(retval))\
                err__finish

	/*@notfunction@*/
#       define err_msg_4(type,func,p1,p2,p3,p4,retval) \
                err__start \
                        return err_macro_4((type), func, (p1), (p2), (p3), \
                                (p4), (retval)) \
                err__finish

	/*@notfunction@*/
#       define err_msg_5(type,func,p1,p2,p3,p4,p5,retval) \
                err__start \
                        return err_macro_5((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (retval)) \
                err__finish

        /*@notfunction@*/
#       define err_msg_6(type,func,p1,p2,p3,p4,p5,p6,retval) \
                err__start \
                        return err_macro_6((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (p6), (retval)) \
                err__finish

        /*@notfunction@*/
#       define err_msg_7(type,func,p1,p2,p3,p4,p5,p6,p7,retval) \
                err__start \
                        return err_macro_7((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (p6), (p7), (retval)) \
                err__finish

        /*@notfunction@*/
#       define err_msg_8(type,func,p1,p2,p3,p4,p5,p6,p7,p8,retval) \
                err__start \
                        return err_macro_8((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (p6), (p7), (p8), (retval)) \
                err__finish

        /*@notfunction@*/
#       define err_msg_9(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,retval) \
                err__start \
                        return err_macro_9((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (p6), (p7), (p8), (p9), (retval)) \
                err__finish

        /*@notfunction@*/
#       define err_msg_10(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,retval) \
                err__start \
                        return err_macro_10((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (p6), (p7), (p8), (p9), (pa), \
                                (retval)) \
                err__finish

        /*@notfunction@*/
#       define err_msg_11(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,retval) \
                err__start \
                        return err_macro_11((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (p6), (p7), (p8), (p9), (pa), (pb),\
                                (retval)) \
                err__finish

        /*@notfunction@*/
#       define err_msg_12(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,retval)\
                err__start \
                        return err_macro_12((type), func, (p1), (p2), (p3), \
                                (p4),(p5),(p6),(p7),(p8),(p9),(pa),(pb),(pc), \
                                (retval)) \
                err__finish

        /*@notfunction@*/
#define err_msg_13(type,func,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,retval) \
                err__start \
                        return err_macro_13((type), func, (p1), (p2), (p3), \
                                (p4), (p5), (p6), (p7), (p8), (p9), (pa), (pb),\
                                (pc), (pd), (retval)) \
                err__finish

#endif
