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
*	Authors: Ivan Zapreev, Christina Jansen
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
*	Source description: This file contains functions for generating
*		random numbers.
*/
#include "rand_num_generator.h"

#include "rng_app_crypt.h"
#include "rng_prism.h"
#include "rng_ciardo.h"
#include "rng_ymer.h"
#include "rng_gsl.h"

#include "macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/************************************************************************/
/******************************Common part*******************************/
/************************************************************************/

/**
* The typedef for the type of a function for generating non-unif distr random numbers.
*/
typedef double (*PTFGenRandNum)(void *);

/**
* The typedef for the type of a function for generating a seed
* for non-unif distr random numbers.
*/
typedef void (*PTFGenSeed)(void *, unsigned long);

/**
* The typedef for the type of a function for generating an exp rand number
*/
typedef double (*PTFGenExpRandNum)(void *, double lambda);

/* A predeclaration, needed for the initRNG method */
static double generateExpRandNumberNonGSL( void * pMethodGSLExp, double lambda);

/**
* The method initializes the main variables for any given RNG method.
* @param _method one from {RNG_APP_CRYPT_METHOD, RNG_PRISM_METHOD,
* RNG_CIARDO_METHOD, RNG_YMER_METHOD, RNG_GSL_RANLUX_METHOD, RNG_GSL_LFG_METHOD,
* RNG_GSL_TAUS_METHOD}
* @param ppFGenRandNum the pointer to pointer to random number generation function
* @param ppFGenExpRandNum the pointer to pointer to exp. rand number generation function
*				Should be set only when this function is used to define the method
*				for exp dist rand number. Otherwise should be NULL.
* @param ppFGenSeed the pointer to pointer to seed generation function
* @param pMethodGSL the pointer to pointer to GSL method structure
*/
static void initRNG(const int _method, PTFGenRandNum * ppFGenRandNum, PTFGenSeed * ppFGenSeed,
			PTFGenExpRandNum * ppFGenExpRandNum, void ** ppMethodGSL){
	switch( _method ){
		case RNG_APP_CRYPT_METHOD:
			*ppFGenRandNum = generateRandUnifAC;
			*ppFGenSeed = generateSeedAC;
			if ( ppFGenExpRandNum != NULL ){
				*ppFGenExpRandNum = generateExpRandNumberNonGSL;
			}
			break;
		case RNG_PRISM_METHOD:
			*ppFGenRandNum = generateRandUnifPR;
			*ppFGenSeed = generateSeedPR;
			if ( ppFGenExpRandNum != NULL ){
				*ppFGenExpRandNum = generateExpRandNumberNonGSL;
			}
			break;
		case RNG_CIARDO_METHOD:
			*ppFGenRandNum = generateRandUnifC;
			*ppFGenSeed = generateSeedC;
			if ( ppFGenExpRandNum != NULL ){
				*ppFGenExpRandNum = generateExpRandNumberNonGSL;
			}
			break;
		case RNG_YMER_METHOD:
			*ppFGenRandNum = generateRandUnifYM;
			*ppFGenSeed = generateSeedYM;
			if ( ppFGenExpRandNum != NULL ){
				*ppFGenExpRandNum = generateExpRandNumberNonGSL;
			}
			break;
		case RNG_GSL_RANLUX_METHOD:
		case RNG_GSL_LFG_METHOD:
		case RNG_GSL_TAUS_METHOD:
			*ppFGenRandNum = generateRandUnifGSL;
			*ppFGenSeed = generateSeedGSL;
			/* WARNING: Has to be freed before to avoid memory leaks */
			*ppMethodGSL = setGSLMethod( _method );
			if ( ppFGenExpRandNum != NULL ){
				*ppFGenExpRandNum = generateExpRandNumberGSL;
			}
			break;
		default:
			printf("ERROR: Invalid method for computing non-uniformly distributed random numbers!\n");
                        exit(EXIT_FAILURE);
	}
}

/************************************************************************/
/*************************Discrete Distribution**************************/
/************************************************************************/

/**
* The method for computing non-uniformly distributed random numbers
* for the discrete distribution.
* If no method was chosen, the default is set to RNG_APP_CRYPT_METHOD.
*/
static int method_discrete = RNG_APP_CRYPT_METHOD;

/**
* Contains the pointer to the GSL RNG method used for discrete distribution.
*/
static void * pMethodGSLDiscrete = NULL;

/**
* Defines the function pointer that will be assigned
* to the proper random number generation procedure later.
*/
static PTFGenRandNum pFGenRandNumDiscrete = NULL;

/**
* Defines the function pointer that will be assigned to the proper seeding
* procedure for the chosen random number generator later.
*/
static PTFGenSeed pFGenRandSeedDiscrete = NULL;

/**
* Sets the method for computing non-uniformly distributed random numbers for
* the discrete distribution.
* @param _method one from {RNG_APP_CRYPT_METHOD, RNG_PRISM_METHOD,
* RNG_CIARDO_METHOD, RNG_YMER_METHOD, RNG_GSL_RANLUX_METHOD, RNG_GSL_LFG_METHOD,
* RNG_GSL_TAUS_METHOD}
*/
void setRNGMethodDiscrete(const int _method){
	/* Free the previously allocated memory etc. */
	freeRNGDiscrete();

	/* Set the method */
	method_discrete = _method;

	/* Initialize the main variables */
	initRNG(method_discrete, &pFGenRandNumDiscrete, &pFGenRandSeedDiscrete, NULL, &pMethodGSLDiscrete);

	/* Genereate a new seed */
	generateNewSeedDiscrete();
}

/**
* Returns the selected method for computing non-uniformly distributed random
* numbers for the discrete distribution.
* NOTE: method should be one from {RNG_APP_CRYPT_METHOD, RNG_PRISM_METHOD,
* RNG_CIARDO_METHOD, RNG_YMER_METHOD, RNG_GSL_RANLUX_METHOD, RNG_GSL_LFG_METHOD,
* RNG_GSL_TAUS_METHOD}
*/
int getRNGMethodDiscrete(void) {
	return method_discrete;
}

/**
* This functions generates a new seed for the selected method for
* non-uniformly distributed random numbers for discrete distribution.
* NOTE: The stream for generating random numbers is seeded by the system clock!
* WARNING: Unless GSL is used this regenerates the seed for the exp. distr. rand.
* variables as well!
*/
void generateNewSeedDiscrete(void) {
	IF_SAFETY( pFGenRandSeedDiscrete != NULL )
		pFGenRandSeedDiscrete( pMethodGSLDiscrete, (unsigned)time( NULL) );
	ELSE_SAFETY
		printf("ERROR: The seed generator for the discrete distribution is not set.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* Helper function for choosing non-uniformly distributed random state.
*/
static int returnDistribStateIndex(const double * distribution, int size) {
	/* The cumulative probabilities form already considered states */
	double cumulative_prob = 0;
	/* Uniform distributed random number */
	double unif_rand = pFGenRandNumDiscrete( pMethodGSLDiscrete );
	int i = 0;

	do {
		cumulative_prob += distribution[i];
		i++;
	} while( i < size && cumulative_prob < unif_rand );

	return i - 1;
}

/**
* This function chooses a random value from the given array according
* to the related probability distribution.
* @param values the range of values one would like to choose from
* @param prob the related probability distribution
* @param size the size of the values and prob arrays
*/
int generateRandNumberDiscrete(const int * values, const double * distribution,
                int size)
{
	IF_SAFETY( ( values != NULL ) && ( distribution != NULL )  && ( size != 0 ) )
		IF_SAFETY( pFGenRandNumDiscrete != NULL )
                        return values[returnDistribStateIndex(distribution,
                                                size)];
		ELSE_SAFETY
			printf("ERROR: The random-number generator for the discrete distribution is not set.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: Attempting to simulate a random variable with input NULL parameters.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* This function frees all memory allocated for the random number generator.
*/
void freeRNGDiscrete(void) {
	if( pMethodGSLDiscrete != NULL ){
		freeRNGGSL(pMethodGSLDiscrete);
		pMethodGSLDiscrete = NULL;
	}
}

/************************************************************************/
/***********************Exponential Distribution*************************/
/************************************************************************/

/**
* Contains the pointer to the GSL RNG method used for exponential distribution.
*/
static void * pMethodGSLExp = NULL;

/**
* The method for computing non-uniformly distributed random numbers
* for the exponential distribution.
* If no method was chosen, the default is set to RNG_GSL_TAUS_METHOD.
*/
static int method_exp = RNG_GSL_TAUS_METHOD;

/**
* Defines the function pointer that will be assigned
* to the proper random number generation procedure later.
*/
static PTFGenRandNum pFGenRandNumExp = NULL;

/**
* Defines the function pointer that will be assigned to the proper seeding
* procedure for the chosen random number generator later.
*/
static PTFGenSeed pFGenRandSeedExp = NULL;

/**
* Defines the function pointer that will be assigned to the proper
* exp. rand. number generator.
*/
static PTFGenExpRandNum pFGenExpRandNum = NULL;

/**
* This method is a complement for the generateExpRandNumberGSL method
* that is based on a non-GSL function.
* @param pMethodGSLExp is not used, can be NULL
* @param lambda the lambda parameter of the exponential distribution
*/
static double generateExpRandNumberNonGSL(void * pMethodGSLExp_local,
                double lambda)
{
	IF_SAFETY( pFGenRandNumExp != NULL )
                return -1/lambda * log(pFGenRandNumExp(pMethodGSLExp_local));
	ELSE_SAFETY
		printf("ERROR: The random-number generator for the discrete distribution is not set.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* Sets the method for computing non-uniformly distributed random numbers for
* the exponential distribution.
* @param _method one from {RNG_APP_CRYPT_METHOD, RNG_PRISM_METHOD,
* RNG_CIARDO_METHOD, RNG_YMER_METHOD, RNG_GSL_RANLUX_METHOD, RNG_GSL_LFG_METHOD,
* RNG_GSL_TAUS_METHOD}
*/
void setRNGMethodExp(int _method){
	/* Free the previously allocated memory etc. */
	freeRNGExp();

	/* Set the method */
	method_exp = _method;

	/* Initialize the main variables */
	initRNG(method_exp, &pFGenRandNumExp, &pFGenRandSeedExp, &pFGenExpRandNum, &pMethodGSLExp);

	/* Generate a new seed */
	generateNewSeedExp();
}

/**
* Returns the selected method for computing non-uniformly distributed random
* numbers for the exponential distribution.
* NOTE: method should be one from {RNG_APP_CRYPT_METHOD, RNG_PRISM_METHOD,
* RNG_CIARDO_METHOD, RNG_YMER_METHOD, RNG_GSL_RANLUX_METHOD, RNG_GSL_LFG_METHOD,
* RNG_GSL_TAUS_METHOD}
*/
int getRNGMethodExp(void) {
	return method_exp;
}
/**
* This functions generates a new seed for the selected method for
* non-uniformly distributed random numbers for exponential distribution.
* NOTE: The stream for generating random numbers is seeded by the system clock!
* WARNING: Unless GSL is used this regenerates the seed for the discrete
* rand. variables as well!
*/
void generateNewSeedExp(void) {
	IF_SAFETY( pFGenRandSeedExp != NULL )
		pFGenRandSeedExp( pMethodGSLExp, (unsigned)time( NULL) );
	ELSE_SAFETY
		printf("ERROR: The seed generator for the exponential distribution is not set.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* Get a random exponentially distributed number
* @param lambda: inverse scale of exponential distribution
* @return the random number (exponentially distributed)
*/
inline double generateRandNumberExp(double lambda){
	IF_SAFETY( lambda > 0 )
		IF_SAFETY( pFGenExpRandNum != NULL )
			return pFGenExpRandNum( pMethodGSLExp, lambda );
		ELSE_SAFETY
			printf("ERROR: The random-number generator for the exponential distribution is not set.\n");
                        exit(EXIT_FAILURE);
		ENDIF_SAFETY
	ELSE_SAFETY
		printf("ERROR: Trying to simulate an exponential distribution with lambda == 0.\n");
                exit(EXIT_FAILURE);
	ENDIF_SAFETY
}

/**
* This function frees all memory allocated for the random number generator
* used for the exponential distribution.
*/
void freeRNGExp(void) {
	if( pMethodGSLExp != NULL ){
		freeRNGGSL(pMethodGSLExp);
		pMethodGSLExp = NULL;
	}
}

/****************************************************************************/
/*******************PRINT THE SIMULATION RUNTIME PARAMETERS******************/
/****************************************************************************/

static const char * getTheRNGMethodName(int rng_method_type) {
        const
	char * result_name = NULL;
	switch( rng_method_type ){
		case RNG_APP_CRYPT_METHOD:
			result_name = RNG_METHOD_APP_CRYPT_STR;
			break;
		case RNG_PRISM_METHOD:
			result_name = RNG_METHOD_PRISM_STR;
			break;
		case RNG_CIARDO_METHOD:
			result_name = RNG_METHOD_CIARDO_STR;
			break;
		case RNG_YMER_METHOD:
			result_name = RNG_METHOD_YMER_STR;
			break;
		case RNG_GSL_RANLUX_METHOD:
			result_name = RNG_METHOD_GSL_RANLUX_STR;
			break;
		case RNG_GSL_LFG_METHOD:
			result_name = RNG_METHOD_GSL_LFG_STR;
			break;
		case RNG_GSL_TAUS_METHOD:
			result_name = RNG_METHOD_GSL_TAUS_STR;
			break;
		default :
			result_name = RNG_METHOD_UNKNOWN_STR;
	}
	return result_name;
}

/**
* Print the random-number generator runtime parameters for discrete r.v.
*/
void printRuntimeRNGInfoDiscrete(void) {
	printf(" RNG discrete dist.\t = %s\n", getTheRNGMethodName( method_discrete ) );
}

/**
* Print the random-number generator runtime parameters for exponential r.v.
*/
void printRuntimeRNGInfoExp(void) {
	printf(" RNG exponential dist.\t = %s\n", getTheRNGMethodName( method_exp ) );
}
