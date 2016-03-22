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
*	Source description: Perform Transient Analysis for PCTL/CSL/PRCTL/CSRL
*		- X, U.
*/

#include "transient.h"

#include "transient_dtmc.h"
#include "transient_ctmc.h"
#include "transient_dtmrm.h"
#include "transient_ctmrm.h"
#include "transient_ctmdpi_hd_uni.h"
#include "transient_ctmdpi_hd_non_uni.h"

#include "runtime.h"

/**
* Solve until formula.
* @param type one of: UNBOUNDED_FORM, INTERVAL_FORM
* @param phi satisfaction relation for phi formula.
* @param psi satisfaction relation for psi formula.
* @param subi the left time bound.
* @param supi tyhe right time bound.
* @param isIgnoreFLumping if we want the formula dependent lumping to be ignored
* @return result of the until formula for all states.
*/
double * until(const int type, const bitset *phi, const bitset *psi, double subi, double supi, BOOL isIgnoreFLumping)
{
	const BOOL isLumping = isRunMode(F_DEP_LUMP_MODE) && ( ! isIgnoreFLumping );

	if( isRunMode(CTMC_MODE) || isRunMode(CMRM_MODE))
	{
		/*The CTMC case*/
		if( type == TIME_UNBOUNDED_FORM ){
                        double * result;
                        if ( isLumping )
                                result = unbounded_until_lumping(phi, psi);
                        else
                                result = unbounded_until(phi, psi);
                        if ( NULL == result ) {
                                err_msg_8(err_CALLBY, "until(%d,%p[%d],%p[%d],"
                                        "%g,%g,%d)", type, (const void *) phi,
                                        bitset_size(phi), (const void *) psi,
                                        bitset_size(psi), subi, supi,
                                        isIgnoreFLumping, NULL);
                        }
                        return result;
		}else{
			if( type == TIME_INTERVAL_FORM ){
				if( subi == 0.0 ){
					return ( isLumping ? bounded_until_lumping(phi, psi, supi) : bounded_until(phi, psi, supi) );
				}else{
					if( subi > 0.0 && supi >= subi ){
						return ( isLumping ? interval_until_lumping(phi, psi, subi, supi) : interval_until(phi, psi, subi, supi) );
					}else{
						printf("ERROR: The formula time bounds [%f,%f] are inappropriate.\n", subi, supi);
                                                return (double*) calloc((size_t)
                                                        bitset_size(phi),
                                                        sizeof(double));
					}
				}
			}else{
				printf("ERROR: An unknown Until formula type %d.\n", type);
                                exit(EXIT_FAILURE);
			}
		}
	}else if( isRunMode(DTMC_MODE) || isRunMode(DMRM_MODE)){
		/*The DTMC case*/
		if( type == TIME_UNBOUNDED_FORM ){
			return ( isLumping ? dtmc_unbounded_until_lumping(phi, psi) : dtmc_unbounded_until(phi, psi) );
		}else{
			if( type == TIME_INTERVAL_FORM ){
				if( subi == 0.0 ){
					return ( isLumping ? dtmc_bounded_until_lumping(phi, psi, supi) : dtmc_bounded_until(phi, psi, supi) );
				}else{
					printf("ERROR: This formula with general time bounds is not supported for PCTL.\n");
                                        return (double*) calloc(
                                                (size_t) bitset_size(phi),
                                                sizeof(double));
				}
			}else{
				printf("ERROR: An unknown Until formula type %d.\n", type);
                                exit(EXIT_FAILURE);
			}
		}
	}else if( isRunMode(CTMDPI_MODE) ){
		/* We only support bounded until formulae of the form U[0, t] */
		if( type == TIME_INTERVAL_FORM  && subi == 0.0 ){
		  const int method = get_method_ctmdpi_transient();
		  if (CTMDPI_HD_AUTO_METHOD == method) {
                        const
			NDSparseMatrix *ctmdp = get_mdpi_state_space();
			if (ctmdp->uniform) {
				return ( transient_ctmdpi_hd_uni(supi, psi) );
			} else {
				return ( transient_ctmdpi_hd_non_uni(supi, psi) );
			}
		  } else if (CTMDPI_HD_UNI_METHOD == method) {
			return ( transient_ctmdpi_hd_uni(supi, psi) );
		  } else if (CTMDPI_HD_NON_UNI_METHOD == method) {
			return ( transient_ctmdpi_hd_non_uni(supi, psi) );
		  } else {
			/* should not happen - just to be sure that
			 * we did not forget something */
			printf("ERROR: Invalid method for transient analysis.\n");
                        exit(EXIT_FAILURE);
		  }
		}else{
			printf("ERROR: This Until formula type is not supported for CTMDPI.\n");
                        return (double *) calloc((size_t) bitset_size(phi),
                                        sizeof(double));
		}
	}else{
			printf("ERROR: The logic does not support the Until formula.\n");
                        exit(EXIT_FAILURE);
	}
	/* NOTE: This should not be happening */
	return NULL;
}

/**
* Solve next formula.
* @param type one of: TIME_UNBOUNDED_FORM, TIME_INTERVAL_FORM
* @param phi satisfaction relation for phi formula.
* @param subi sub I (time).
* @param supi sup I (time).
* @return result of the next formula for all states.
*/
double * next(const int type, const bitset * phi, double subi, double supi) {
	if( isRunMode(CTMC_MODE) || isRunMode(CMRM_MODE) ){
		/*The CTMC case*/
		if( type == TIME_UNBOUNDED_FORM ){
			return unbounded_next(phi);
		}else{
			if( type == TIME_INTERVAL_FORM ){
				if( subi == 0.0){
					return bounded_next(phi, supi);
				}else{
					if( subi > 0.0 && supi >= subi ){
						return interval_next(phi, subi, supi);
					}else{
						printf("ERROR: The formula time bounds [%f,%f] are inappropriate.\n", subi, supi);
                                                return (double*) calloc((size_t)
                                                        bitset_size(phi),
                                                        sizeof(double));
					}
				}
			}else{
				printf("ERROR: An unknown Next formula type %d.\n", type);
                                exit(EXIT_FAILURE);
			}
		}
	}else{
		if( isRunMode(DTMC_MODE) || isRunMode(DMRM_MODE)  ){
			/*The DTMC case*/
			if( type == TIME_UNBOUNDED_FORM ){
				return dtmc_unbounded_next(phi);
			}else{
				printf("ERROR: Formula is not supported by PCTL.\n");
                                return (double *) calloc((size_t)
                                        bitset_size(phi), sizeof(double));
			}
		}else{
			printf("ERROR: The logic does not support the Next formula.\n");
                        exit(EXIT_FAILURE);
		}
	}
	/* NOTE: This should not be happening */
	return NULL;
}

/**
* Solve until formula with rewards.
* @param phi satisfaction relation for phi formula.
* @param psi satisfaction relation for psi formula.
* @param subi sub I (time).
* @param supi sup I (time).
* @param subj sub J (reward).
* @param supj sup J (reward).
* @param isIgnoreFLumping if we want the formula dependent lumping to be ignored
* @param ppResultError  the pointer to the array of doubles. This array will
*                       store the errors
*			for the resulting probabilities. (This is the return variable)
* @return result of the until formula for all states.
*/
double * until_rewards(const bitset * phi, const bitset * psi, double subi,
                double supi, double subj, double supj, BOOL isIgnoreFLumping,
                double ** ppResultError)
{
	const BOOL isLumping = isRunMode(F_DEP_LUMP_MODE) && ( ! isIgnoreFLumping );

	if( isRunMode(CMRM_MODE) ){
		if( subi == 0.0 && subj == 0.0 && supi != 0.0 && supj != 0.0 ){
			return ( isLumping ? ctmrm_bounded_until_lumping(phi, psi, supi, supj, ppResultError) : ctmrm_bounded_until(phi, psi, supi, supj, ppResultError) );
		}else{
			printf("ERROR: Given parameters combination is not supported by CSRL.\n");
                        return (double *) calloc((size_t) bitset_size(phi),
                                        sizeof(double));
		}
	}else{
		if( isRunMode(DMRM_MODE) ){
			if( supi != 0.0 && supj != 0.0 ){
				return ( isLumping ? dtmrm_bounded_until_lumping(phi, psi, subi, supi, subj, supj) : dtmrm_bounded_until(phi, psi, subi, supi, subj, supj) );
			}else{
				printf("ERROR: Given parameters combination is not supported by PRCTL.\n");
                                return (double *) calloc((size_t)
                                        bitset_size(phi), sizeof(double));
			}
		}else{
			printf("ERROR: The logic does not support the Until formula.\n");
                        exit(EXIT_FAILURE);
		}
	}
	/* NOTE: This should not be happening */
	return NULL;
}

/**
* Solve next_rewards formula.
* @param phi satisfaction relation for phi formula.
* @param subi sub I (time).
* @param supi sup I (time).
* @param subj sub J (reward).
* @param supj sup J (reward).
* @return result of the next formula for all states.
* NOTE: ********without impulse rewards at present***********
*		  ********to be tested***********
*/
double * next_rewards(const bitset * phi, double subi, double supi, double subj,
                double supj)
{
	if( isRunMode(CMRM_MODE) ){
		/*The CMRM case*/
		if( supi != 0.0 ){
			return ctmrm_interval_next(phi, subi, supi, subj, supj);
		}else{
			printf("ERROR: Given parameters combination is not supported by CSRL.\n");
                        return (double *) calloc((size_t) bitset_size(phi),
                                        sizeof(double));
		}
	}else{
		if( isRunMode(DMRM_MODE) ){
			/*Only unbounded next is allowed (see parser for
			restriction of other cases) - The DMRM case*/
			if(subi==0 && supi==0 && subj==0 && supj==0){
				return dtmc_unbounded_next(phi);
			}else{
				printf("ERROR: Formula is not supported by PRCTL.\n");
                                return (double *) calloc((size_t)
                                        bitset_size(phi), sizeof(double));
			}
		}else{
			printf("ERROR: The logic does not support the Next formula.\n");
                        exit(EXIT_FAILURE);
		}
	}
}
