/**
*	WARNING: Do Not Remove This Section
*
*       $LastChangedRevision: 434 $
*       $LastChangedDate: 2011-01-12 09:13:25 +0100 (Mi, 12. Jan 2011) $
*       $LastChangedBy: nguyen $
*
*	MRMC is a model checker for discrete-time and continuous-time Markov
*	reward models. It supports reward extensions of PCTL and CSL (PRCTL
*	and CSRL), and allows for the automated verification of properties
*	concerning long-run and instantaneous rewards as well as cumulative
*	rewards.
*
*	Copyright (C) The University of Twente, 2004-2008.
*	Copyright (C) RWTH Aachen, 2008-2009.
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
*	Source description: Contains the main function, starts the execution,
*	does initial validation of input parameters, loads the required input
*	files.
*/

# include "runtime.h"
# include "read_rewards.h"
# include "read_tra_file.h"
# include "read_lab_file.h"
# include "read_mdpi_file.h"
# include "read_impulse_rewards.h"
# include "execute_cmd_script.h"
# include "write_res_file.h"
# include "lump.h"
# include "parser_to_core.h"
# include "steady.h"
#include "rand_num_generator.h"

#include <string.h>
#include <time.h>
#ifndef __APPLE__
#       include <malloc.h>
#endif

extern
int yyparse (void);

/* This is the list of possible models */
#define CTMC_MODE_STR "ctmc"
#define DTMC_MODE_STR "dtmc"
#define DMRM_MODE_STR "dmrm"
#define CMRM_MODE_STR "cmrm"
#define CTMDPI_MODE_STR "ctmdpi"

/* This is the list of possible options */
#define F_IND_LUMP_MODE_STR "-ilump"
#define F_DEP_LUMP_MODE_STR "-flump"

/* This "logic" is used for testing vector */
/* matrix and matrix vector multiplications */
#define TEST_VMV_MODE_STR "vmv_test"

/* This is the list of possible file extensions */
#define TRA_FILE_EXT ".tra"
#define LAB_FILE_EXT ".lab"
#define REW_FILE_EXT ".rew"
#define CMD_FILE_EXT ".cmd"
#define RES_FILE_EXT ".res"
#define REWI_FILE_EXT ".rewi"
#define CTMDPI_FILE_EXT ".ctmdpi"

/**
* An extension can be one of:
*	.rew, .rewi, .tra, .lab, .ctmdpi
* plus at least one symbol of the name
*/
#define MIN_FILE_NAME_LENGTH 5
/**
* Here once again the max length of the extensions:
*	.rew, .rewi, .tra, .lab, .ctmdpi
* is 5 symbols, this we need to copy the extension out
* of the possible file name.
*/
#define MAX_FILE_EXT_LENGTH 7
#define MIN_FILE_EXT_LENGTH 4

/* These are the globals used only in this file. */
static sparse *space = NULL;		/* The statespace */
static NDSparseMatrix *mdpi = NULL;	/* The MDP (int. non-det) statespace */
static labelling *labels = NULL;	/* The labelling */

/**
* This is a list of boolean variables which indicate
* whether certain files were loaded or not
*/
static BOOL is_tra_present  = FALSE;
static BOOL is_lab_present  = FALSE;
static BOOL is_rew_present  = FALSE;
static BOOL is_cmd_present  = FALSE;
static BOOL is_res_present  = FALSE;
static BOOL is_rewi_present = FALSE;
static BOOL is_ctmdpi_present = FALSE;

/**
* Here we will store pointers to the input files
*/
static const char * tra_file  = NULL;
static const char * lab_file  = NULL;
static const char * rew_file  = NULL;
static const char * cmd_file  = NULL;
extern const char * res_file;
static const char * rewi_file = NULL;
static const char * ctmdpi_file = NULL;

/**
* This part simply prints the program usage info.
*/
static void usage(void)
{
	printf("Usage: mrmc <model> <options> <.tra file> <.ctmdpi file> <.lab file> <.rew file> <.rewi file> <.cmd file> <.res file>\n");
	printf("\t<model>\t\t- could be one of {%s, %s, %s, %s, %s}.\n",CTMC_MODE_STR, DTMC_MODE_STR, DMRM_MODE_STR, CMRM_MODE_STR, CTMDPI_MODE_STR);
	printf("\t<options>\t- could be one of {%s, %s}, optional.\n", F_IND_LUMP_MODE_STR, F_DEP_LUMP_MODE_STR);
	printf("\t<.tra file>\t- is the file with the matrix of transitions (for DMRM/CMRM, DTMC/CTMC).\n");
	printf("\t<.ctmdpi file>\t- is the file with the transition matrix and transition labels (for CTMDPI).\n");
	printf("\t<.lab file>\t- contains labeling.\n");
	printf("\t<.rew file>\t- contains state rewards (for DMRM/CMRM).\n");
	printf("\t<.rewi file>\t- contains impulse rewards (for CMRM, optional).\n");
	printf("\t<.cmd file>\t- contains script to execute (optional).\n");
	printf("\t<.res file>\t- filename where write_res_file writes the results to (optional).\n");
	printf("\nNote: In the '.tra' and '.ctmdpi' file transitions should be ordered by rows and columns!\n\n");
}

/**
* This part simply prints the programm Intro message
*/
static void print_intro(void)
{
	printf(" ------------------------------------------------------------------- \n");
	printf("|                    Markov Reward Model Checker                    |\n");
	printf("|                         MRMC Version 1.5                          |\n");
	printf("|               Copyright (C) RWTH Aachen, 2006-2011.               |\n");
	printf("|         Copyright (C) The University of Twente, 2004-2008.        |\n");
	printf("|                            Authors:                               |\n");
	printf("|      Ivan S. Zapreev (2004-2011), David N. Jansen (since 2007),   |\n");
	printf("|         E. Moritz Hahn (2007-2010), Falko Dulat (2009-2010),      |\n");
	printf("|         Christina Jansen (2007-2008), Sven Johr (2006-2007),      |\n");
	printf("|          Tim Kemna (2005-2006), Maneesh Khattri (2004-2005)       |\n");
	printf("|           MRMC is distributed under the GPL conditions            |\n");
	printf("|           (GPL stands for GNU General Public License)             |\n");
	printf("|          The product comes with ABSOLUTELY NO WARRANTY.           |\n");
	printf("|  This is a free software, and you are welcome to redistribute it. |\n");
	printf(" ------------------------------------------------------------------- \n");
}

/**
* This method simply tests matrix vector and vector matrix multiplication
* @param space this is an implicit parameter, which is a matrix to test
*/
static void testVMV(void)
{
	int i;
        const int N = mtx_rows(space);
	double *pV,*pR;
        unsigned long elapsed_time;
	/* Allocate arrays */
        pV = (double *) calloc((size_t) N, sizeof(double));
        pR = (double *) calloc((size_t) N, sizeof(double));
	/* Initialize the vector */
	for(i=0; i < N; i++)
	{
		pV[i] = 2.0;
	}

	/* Do computations M*v */
        startTimer();
        if ( err_state_iserror(multiply_mtx_MV(space, pV, pR)) ) {
                exit(err_macro_0(err_CALLBY, "testVMV()", (free(pR), free(pV),
                                free_sparse_ncolse(space), EXIT_FAILURE)));
        }

        elapsed_time = stopTimer();
        printf("Time for M*v: %lu000 micro sec(s)\n", elapsed_time);

	/* Do computations v*M */
        startTimer();
        if ( err_state_iserror(multiply_mtx_TMV(space, pV, pR)) ) {
                exit(err_macro_0(err_CALLBY, "testVMV()", (free(pR), free(pV),
                                free_sparse_ncolse(space), EXIT_FAILURE)));
        }

        elapsed_time = stopTimer();
        printf("Time for v*M: %lu000 micro sec(s)\n", elapsed_time);
	free(pV);
	free(pR);
        if ( err_state_iserror(free_sparse_ncolse(space)) ) {
                exit(err_macro_0(err_CALLBY, "testVMV()", EXIT_FAILURE));
        }
}

/**
* This method determines the running mode depending on the command line parameters
* It checks whether we will work with CTMC, DTMC, CMRM or DMRM.
* @return FALSE if it failed to recognize the input parameter as am MRMC
*	   mode parameter otherwise TRUE
*/
static BOOL setRunningMode(const char * mode)
{
	BOOL result = FALSE;

	if( strcmp(mode,CTMC_MODE_STR) == 0 ){
		result = TRUE;
		if( !isRunModeSet() )
		{
			addRunMode(CTMC_MODE);
			print_run_mode( TRUE, FALSE );
		}else{
			printf("WARNING: The model type has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}else if( strcmp(mode,DTMC_MODE_STR) == 0  && !isRunModeSet() ){
		result = TRUE;
		if( !isRunModeSet() )
		{
			addRunMode(DTMC_MODE);
			print_run_mode( TRUE, FALSE );
		}else{
			printf("WARNING: The model type has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}else if( strcmp(mode,DMRM_MODE_STR) == 0  && !isRunModeSet() ){
		result = TRUE;
		if( !isRunModeSet() )
		{
			addRunMode(DMRM_MODE);
			print_run_mode( TRUE, FALSE );
		}else{
			printf("WARNING: The model type has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}else if( strcmp(mode,CMRM_MODE_STR) == 0  && !isRunModeSet() ){
		result = TRUE;
		if( !isRunModeSet() )
		{
			addRunMode(CMRM_MODE);
			print_run_mode( TRUE, FALSE );
		}else{
			printf("WARNING: The model type has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}else if( strcmp(mode,CTMDPI_MODE_STR) == 0 && !isRunModeSet() ){
		result = TRUE;
		if( !isRunModeSet() )
		{
			addRunMode(CTMDPI_MODE);
			print_run_mode( TRUE, FALSE);
		}else{
			printf("WARNING: The model type has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}else if( strcmp(mode,TEST_VMV_MODE_STR) == 0 && !isRunModeSet() ){
		result = TRUE;
		if( !isRunModeSet() )
		{
			addRunMode(TEST_VMV_MODE);
			print_run_mode( TRUE, FALSE );
		}else{
			printf("WARNING: The model type has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}else if( strcmp(mode,F_IND_LUMP_MODE_STR) == 0 ){
		result = TRUE;
		if( !isRunMode(F_DEP_LUMP_MODE) ){
			addRunMode(F_IND_LUMP_MODE);
			print_run_mode( FALSE, TRUE );
		}else{
			printf("WARNING: The lumping mode has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}else if( strcmp(mode,F_DEP_LUMP_MODE_STR) == 0 ){
		result = TRUE;
		if( !isRunMode(F_IND_LUMP_MODE) ){
			addRunMode(F_DEP_LUMP_MODE);
			print_run_mode( FALSE, TRUE );
		}else{
			printf("WARNING: The lumping mode has already been set, the parameter '%s' is skipped.\n", mode);
		}
	}
	return result;
}

/**
* This function should find the extension and then return it and its length.
* @param p_pot_filename the input parameter which we suspect to be a file name
* @param p_extension the pointer to an extension string (a return parameter)
* @param p_ext_length the pointer to the extension string length (a return parameter)
* @return TRUE if we can subtract a file extension which we think is one of required
*/
static BOOL isValidExtension(const char * p_pot_filename, char * p_extension,
                int * p_ext_length)
{
	BOOL result = FALSE;
        const
	char *p;
	int length;

	length = strlen(p_pot_filename);
	p = strrchr(p_pot_filename,'.');     /* The last occurance of '.' in the file name */
	*p_ext_length = length - (p - p_pot_filename); /* including '.', excluding '\0' */

	/** @todo change this! */
	if( length >= MIN_FILE_NAME_LENGTH &&
	    ( *p_ext_length >= MIN_FILE_EXT_LENGTH) &&
	    ( *p_ext_length <= MAX_FILE_EXT_LENGTH)) {
                /* Get the extension */
		if (*p_ext_length == ( MAX_FILE_EXT_LENGTH - 2)){
			strncpy( p_extension, p , MAX_FILE_EXT_LENGTH - 1);
		}else if( *p_ext_length == ( MAX_FILE_EXT_LENGTH - 1 ) ){
			strncpy( p_extension, p , MAX_FILE_EXT_LENGTH );
		}else{
			strncpy( p_extension, p , MAX_FILE_EXT_LENGTH + 1 );
		}

		result = TRUE;
		strncpy( p_extension, p, MAX_FILE_EXT_LENGTH);
	}

	return result;
}

/**
* Load the .rewi impulse rewards file
*/
static void loadImpulseRewards(const char * file_name) {
	if( is_rewi_present ){
		if ( isRunMode(CMRM_MODE) )
		{
			sparse * rewi;

			printf("Loading the '%s' file, please wait.\n", file_name);
                        rewi = read_impulse_rewards(file_name, mtx_rows(space));
			if ( rewi == NULL )
			{
				printf("ERROR: The '%s' file '%s' was not found.\n", REWI_FILE_EXT, file_name);
                                exit(EXIT_FAILURE);
			}
			setImpulseRewards(rewi);
		}else{
			printf("WARNING: Impulse rewards are not supported in the current mode, skipping the '%s' file loading.\n", file_name);
		}
	}
}

/**
* Load the .rew state rewards file
*/
static void loadStateRewards(const char * file_name) {
	if( is_rew_present ){
		if ( isRunMode(DMRM_MODE) || isRunMode(CMRM_MODE) )
		{
			double * rew;

			printf("Loading the '%s' file, please wait.\n", file_name);
                        rew = read_rew_file(mtx_rows(space), file_name);
			if ( rew == NULL )
			{
				printf("ERROR: The '%s' file '%s' was not found.\n", REW_FILE_EXT, file_name);
                                exit(EXIT_FAILURE);
			}
			setStateRewards(rew);
		}else{
			printf("WARNING: State rewards are not supported in the current mode, skipping the '%s' file loading.\n", file_name);
		}
	}
}

/**
* Load the .ctmdpi transitions file
*/
static void loadCTMDPI(const char * file_name) {
	if(is_ctmdpi_present) {
		printf("Loading the '%s' file, please wait.\n", file_name);
		mdpi = read_NDSparseMatrix_file(file_name);
		if(mdpi == NULL ){
			printf("ERROR: The '%s' file '%s' was not found or is incorrect\n", CTMDPI_FILE_EXT, file_name);
                        exit(EXIT_FAILURE);
		}
		set_mdpi_state_space(mdpi);
	}
}

/**
* Load the .tra transitions file
*/
static void loadTransitions(const char * file_name) {
	if( is_tra_present ){
		printf("Loading the '%s' file, please wait.\n", file_name);
		space = read_tra_file(file_name);
		if( space == NULL )
		{
			printf("ERROR: The '%s' file '%s' was not found.\n", TRA_FILE_EXT, file_name);
                        exit(EXIT_FAILURE);
		}
		set_state_space(space);
	}
}

/**
* Load the .lab labelling file
* WARNING: Presumes that transitions have been loaded and
* sparse variable is properly initialized
*/
static void loadLabels(const char * file_name) {
	if( is_lab_present ){
		printf("Loading the '%s' file, please wait.\n", file_name);

		if( isRunMode(CTMDPI_MODE) ){
			labels = read_lab_file(mdpi->n, file_name);
		}else{
                        labels = read_lab_file(mtx_rows(space), file_name);
		}
		if(labels == NULL)
		{
			printf("ERROR: The '%s' file '%s' was not found.\n", LAB_FILE_EXT, file_name);
                        exit(EXIT_FAILURE);
		}
		set_labeller(labels);
	}
}

/**
* This method is used to check that all required files
* are present in the command line parameters and also that other
* command line options are used in valid combinations.
*/
static void checkConsistency(void) {
        const
	char *missing_file;

	/* Check for the model type present */
	if( !isRunModeSet() ){
		printf("ERROR: The <model> parameter is undefined.\n");
		usage();
                exit(EXIT_FAILURE);
	}

	/* Check for the option combinations. */
	/* Formula independent lumping for CSRL is not supported, the option will be ignored. */
	if( isRunMode(CMRM_MODE) && isRunMode(F_IND_LUMP_MODE) && is_rewi_present){
		printf("WARNING: Formula independent lumping for CMRM with impulse rewards is not supported, the '%s' option will be ignored.\n",F_IND_LUMP_MODE_STR);
		clearRunMode(F_IND_LUMP_MODE);
		printf("The formula independent lumping is OFF.\n");
	}
	/* Formula dependent lumping in CSRL is not supported for the case with impulse rewards */
	if( isRunMode(CMRM_MODE) && isRunMode(F_DEP_LUMP_MODE) && is_rewi_present ){
		printf("WARNING: Formula dependent lumping for CMRM with impulse rewards is not supported, the '%s' option will be ignored.\n",F_DEP_LUMP_MODE_STR);
		clearRunMode(F_DEP_LUMP_MODE);
		printf("The formula dependent lumping is OFF.\n");
	}
	/* Formula independent lumping in PRCTL is not supported for the case with impulse rewards */
	if( isRunMode(DMRM_MODE) && isRunMode(F_IND_LUMP_MODE) && is_rewi_present ){
		printf("WARNING: Formula independent lumping for DMRM with impulse rewards is not supported, the '%s' option will be ignored.\n",F_IND_LUMP_MODE_STR);
		clearRunMode(F_IND_LUMP_MODE);
		printf("The formula independent lumping is OFF.\n");
	}
	/* Formula dependent lumping in PRCTL is not supported for the case with impulse rewards */
	if( isRunMode(DMRM_MODE) && isRunMode(F_DEP_LUMP_MODE) && is_rewi_present ){
		printf("WARNING: Formula dependent lumping for DMRM with impulse rewards is not supported, the '%s' option will be ignored.\n",F_DEP_LUMP_MODE_STR);
		clearRunMode(F_DEP_LUMP_MODE);
		printf("The formula dependent lumping is OFF.\n");
	}
	/* Formula dependent lumping for CTMDPI is not supported, the option will be ignored. */
	if( isRunMode(CTMDPI_MODE) && isRunMode(F_DEP_LUMP_MODE) ){
		printf("WARNING: Formula dependent lumping for CTMDPI is not supported, the '%s' option will be ignored.\n",F_DEP_LUMP_MODE_STR);
		clearRunMode(F_DEP_LUMP_MODE);
		printf("The formula dependent lumping is OFF.\n");
	}
	/* Formula independent lumping for CTMDPI is not supported, the option will be ignored. */
	if( isRunMode(CTMDPI_MODE) && isRunMode(F_IND_LUMP_MODE) ){
		printf("WARNING: Formula independent lumping for CTMDPI is not supported, the '%s' option will be ignored.\n",F_IND_LUMP_MODE_STR);
		clearRunMode(F_IND_LUMP_MODE);
		printf("The formula independent lumping is OFF.\n");
	}

	/* Check for the presence of all required files. */
	missing_file = NULL;
	if( !isRunMode(CTMDPI_MODE) && !is_tra_present ){
		missing_file = TRA_FILE_EXT;
	}else if( isRunMode(CTMDPI_MODE) && !is_ctmdpi_present ){
		missing_file = CTMDPI_FILE_EXT;
	}else if( !is_lab_present ){
		missing_file = LAB_FILE_EXT;
	}else if( ( isRunMode(DMRM_MODE) || isRunMode(CMRM_MODE) ) && !is_rew_present ){
		missing_file = REW_FILE_EXT;
	}
	if( missing_file != NULL ){
		printf("ERROR: The '%s' file is required to be present.\n", missing_file);
		usage();
                exit(EXIT_FAILURE);
	}
}

/**
* This method is used for parsing the input paramters of the MRMC tool
* and also for loading the data files.
* NOTE: This method terminates the programm in case of incorrect input parameters
*/
static void parseInParamsAndInitialize(int argc, char *argv[])
{
	int i = 0, ext_length;
	char expension[MAX_FILE_EXT_LENGTH+1]; /* +1 because we need also to store the \0 symbol */

	/* NOTE: Note it is important to load files in the following order: */
	/*	  .tra, .lab, .rew, .rewi */
	/* So we first parse the input parameters, sort them out and only */
	/* then read the files. */
	for(i = 1; i < argc; i++)
	{
		if( !setRunningMode( argv[i] ) ){
			if( isValidExtension( argv[i], expension, &ext_length ) ){
				if( strcmp(expension, TRA_FILE_EXT) == 0 ){
					if( !is_tra_present ){
							is_tra_present = TRUE;
							tra_file = argv[i];
					}else{
						printf("WARNING: The '%s' file has been noticed before, skipping the '%s' file.\n", tra_file, argv[i]);
					}
				}else if ( strcmp(expension, LAB_FILE_EXT) == 0 ){
					if( !is_lab_present ){
							is_lab_present = TRUE;
							lab_file = argv[i];
					}else{
						printf("WARNING: The '%s' file has been noticed before, skipping the '%s' file.\n", lab_file, argv[i]);
					}
				}else if ( strcmp(expension, REW_FILE_EXT) == 0 ){
					if( !is_rew_present ){
							is_rew_present = TRUE;
							rew_file = argv[i];
					}else{
						printf("WARNING: The '%s' file has been noticed before, skipping the '%s' file.\n", rew_file, argv[i]);
					}
				}else if ( strcmp(expension, REWI_FILE_EXT) == 0 ){
					if( !is_rewi_present ){
							is_rewi_present = TRUE;
							rewi_file = argv[i];
					}else{
						printf("WARNING: The '%s' file has been noticed before, skipping the '%s' file.\n", rewi_file, argv[i]);
					}
				}else if ( strcmp(expension, CTMDPI_FILE_EXT) == 0) {
					if( !is_ctmdpi_present){
						is_ctmdpi_present = TRUE;
						ctmdpi_file = argv[i];
					}else{
						printf("WARNING: The '%s' file has been noticed before, skipping the '%s' file.\n", ctmdpi_file, argv[i]);
					}
				}else if ( strcmp(expension, CMD_FILE_EXT) == 0 ){
					if( !is_cmd_present ){
							is_cmd_present = TRUE;
							cmd_file = argv[i];
					}else{
						printf("WARNING: The '%s' file has been noticed before, skipping the '%s' file.\n", cmd_file, argv[i]);
					}
				}else if ( strcmp(expension, RES_FILE_EXT) == 0 ){
					if( !is_res_present ){
							is_res_present = TRUE;
							res_file = argv[i];
					}else{
						printf("WARNING: The '%s' file has been noticed before, skipping the '%s' file.\n", res_file, argv[i]);
					}
				}else {
				    printf("WARNING: An unknown file type '%s' for input file '%s', skipping.\n", expension, argv[i]);
				}
			} else {
				printf("WARNING: An unknown model, option or file '%s', skipping.\n", argv[i]);
			}
		}
	}

        /* Check that the parameters are consistent, */
	/* may cause the program termination. */
        checkConsistency();

	/* Load the files if present */
	loadTransitions( tra_file );
	loadCTMDPI( ctmdpi_file );
	loadLabels( lab_file );
	loadStateRewards( rew_file );
	loadImpulseRewards( rewi_file);
	
	/* Initialize the array to store the requested states to write */
	write_res_file_initialize();

	/* If we are in the test mode then just test and exit. */
	if( isRunMode(TEST_VMV_MODE) )
	{
		/* The tra file has to be preloaded already */
		testVMV();
                exit(EXIT_SUCCESS);
	}
}

/**
* Does formula independent lumping if the isRunMode(F_IND_LUMP_MODE_STR) == TRUE
*/
static void doFormulaIndependentLumping(void) {
	if( isRunMode( F_IND_LUMP_MODE ) ){
		partition * P;
		sparse *Q1;
		double * p_old_rew, * p_new_rew;
                unsigned
		long elapsed_time;

		printf("Formula Independent Lumping...\n");
		/* Starting the timer, if it has not been started yet */
		startTimer();

		/* Create initial partitioning */
		P = init_partition(labels);

		/* Lump the state space */
		Q1 = lump(P, space);

		/* Change labelling */
		change_labelling(labels, P);

		/* Change state rewards */
		p_old_rew = getStateRewards();
		p_new_rew = change_state_rewards(p_old_rew, P);
		freeStateRewards();
		setStateRewards(p_new_rew);

		/* free the original state space */
                /* if ( ! isRunMode(CTMDPI_MODE) ) { */
                        if ( err_state_iserror(free_sparse_ncolse(space)) ) {
                                exit(err_macro_0(err_CALLBY,
                                        "doFormulaIndependentLumping()",
                                        EXIT_FAILURE));
                        }
                /* } */
		free_row_sums();
		/* Reset the state space */
		space = Q1;
		/* Store into the runtime.c */
		set_state_space(space);
		setPartition(P);

		/* Stop timer and print the elapsed time */
		elapsed_time = stopTimer();
		printf("The Total Elapsed Lumping Time is %ld milli sec(s).\n", elapsed_time);
	}
}

int main(int argc, char *argv[])
{
#ifndef __APPLE__
	unsigned int sp1=0,sp2=0;
#endif
	int exitAfterCmdScript = 0;

	/* Pring the introduction message */
	print_intro();

#ifndef __APPLE__
	sp1 = mallinfo().uordblks;
#endif

	/* Parse and validate the input parameters, */
	/* load required files */
	parseInParamsAndInitialize(argc, argv);

#       ifndef __APPLE__
		sp2 = mallinfo().uordblks;

		printf("The Occupied Space is %d Bytes.\n", (sp2-sp1));
#       else
		printf("The Occupied Space is ??? Bytes.\n");
#       endif

	/* Do formula independent lumping, if required */
	doFormulaIndependentLumping();

	printf("Type 'help' to get help.\n>>");
	
	if (is_cmd_present)
	{
		exitAfterCmdScript = execute_cmd_script(cmd_file);
	}
	
	if (exitAfterCmdScript==0)
	{
		while( yyparse() )
		{
			printf(">>");
		}
	}

	if ( ! isRunMode(CTMDPI_MODE) ) {
			if ( err_state_iserror(free_sparse_ncolse(space)) ) {
					err_msg_2(err_CALLBY, "main(%d,%p)", argc, (void*) argv,
							EXIT_FAILURE);
			}
	}
	space = NULL; set_state_space(NULL);
	/* NOTE: At the moment there is no need to call free_row_sums */
	/* method because it is done when we do set_state_space(NULL)*/
	free_labelling(labels); labels = NULL; set_labeller(NULL);

	/* If something was allocated for impulse rewards */
	freeImpulseRewards();
	/* If something was allocated for state rewards */
	freeStateRewards();
	/* If something was allocated for lumping */
	freePartition();

	/* If something was allocated for model checkig the steady-state operator */
	freeSteady();

	/*This is done to deallocate the memory used by RNG methods.*/
	/* Free the random-number generator data, especially needed by GSL functions */
	freeRNGDiscrete();
	freeRNGExp();

	freeNDSparseMatrix();

        return EXIT_SUCCESS;
}
