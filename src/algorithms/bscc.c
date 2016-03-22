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
*	Source description: This file contains functions for the BSCCs
*		detection. The algorithm is based on the Tarjan's algorithm
*		for searching SCCs.
*/

#include "bscc.h"

#include "stack.h"

#include "runtime.h"

/*
The stack array structure that stores all the needed stacks,
The bscc_stack is one of its elements
*/
static stack *stackarray = NULL;

/*Stack needed for the bscc search */
static int *bscc_stack = NULL;

/******************THE BSCCs SEARCH FUNCTIONS*********************************/

/*This is a holder for the rate matrix of the CTMC*/
static const sparse * pStateSpace = NULL;

/*This bitset stores the visited states;*/
static bitset *pVisitedStates = NULL;

/*This bitset stores states that belong to some component
or from which we can reach a component;*/
static bitset *pInComponentStates = NULL;

/*This array stores mapping from node ids to the BSCC ids
The index of the array is the index of node and value in the corresponding
Cell is the id of the BSCC to which the given node belongs
If it is 0 then the BSCC was not found yet or it doesn't belong to any BSCC*/
static int * pBSCCs = NULL;

/*The skip variable is used to exit from the recursion
in the visit(int) function when the search should be stopped.*/
static BOOL bGlobalSkip = FALSE;

/*This variable stores the DFS order*/
static int dfs_order = 1;

/*This variable stores the id of the lately found bscc*/
static int bscc_counter = 0;

/**
* @param pBSCCsHolder the structure to retrieve bscc_counter from
* @return the number of BSCCs found so far.
*/
int getBSCCCounter(const TBSCCs * pBSCCsHolder)
{
	return pBSCCsHolder->bscc_counter;
}

/**
* This method returns a mapping array from states into the BSCCs' ids.
* @param pBSCCsHolder the structure to retrieve BSCCs mapping from
* @return The array which stores mapping from node ids into the BSCC ids
* The index of the array is the index of node and value in the corresponding
* Cell is the id of the BSCC to which the given node belongs
* If it is 0 then the BSCC was not found yet or it doesn't belong to any BSCC
*/
const int * getStatesBSCCsMapping(const TBSCCs * pBSCCsHolder)
{
	return pBSCCsHolder->pBSCCs;
}

/**
* This function simply instantiates and returns the pointer to TBSCCs struct
* @param pStateSpaceTmp a state space to work with
* @return an empty TBSCCs structure
*/
TBSCCs * allocateTBSCCs(const sparse * pStateSpaceTmp)
{
	/*Allocate structure*/
        TBSCCs * pBSCCsHolder = (TBSCCs *) calloc((size_t) 1, sizeof(TBSCCs));

        pBSCCsHolder->pStateSpace = pStateSpaceTmp;
        pBSCCsHolder->size = mtx_rows(pStateSpaceTmp);

	/*Create a visited states set if this is the first search*/
	pBSCCsHolder->pVisitedStates = get_new_bitset(pBSCCsHolder->size);
	pBSCCsHolder->pInComponentStates = get_new_bitset(pBSCCsHolder->size);
        pBSCCsHolder->pBSCCs = (int *) calloc((size_t) pBSCCsHolder->size,
                        sizeof(int));
	pBSCCsHolder->dfs_order = 1;
	pBSCCsHolder->bscc_counter = 0;

	return pBSCCsHolder;
}

/**
* This method is used to deallocate TBSCCs and all associated data.
* pBSCCsHolder the pointer to the structure to be cleaned
*/
void freeTBSCC(TBSCCs * pBSCCsHolder)
{
	if( pBSCCsHolder ) {
		if( pBSCCsHolder->pInComponentStates ) {
			free_bitset(pBSCCsHolder->pInComponentStates);
		}
		if( pBSCCsHolder->pVisitedStates ) {
			free_bitset(pBSCCsHolder->pVisitedStates);
		}
		if( pBSCCsHolder->pBSCCs ) {
			free(pBSCCsHolder->pBSCCs);
		}
		free(pBSCCsHolder);
	}
}

/**
* This functions checks whether the given state was already visited or not.
* @param i the state id
* @return true if the state was not yet visited.
*/
static inline BOOL isNotVisited(int i)
{
	return ( get_bit_val( pVisitedStates, i ) ? FALSE : TRUE );
}

/**
* Marks the state i as visitedfreeStack()
* @param ippBSCCs the state that was visited
*/
static inline void setVisited(int i)
{
	set_bit_val(pVisitedStates, i, BIT_ON);
}

/***************THE ROOT ARRAY MANAGEMENT PROCEDURES***************************/

/*This array is used to store the root values of the nodes*/
static int * pRoot = NULL;

/**
* This method initializes the root array of length size
*/
static inline void initRoot(int size)
{
        pRoot = (int *) calloc((size_t) size, sizeof(int));
}

/**
* This method is used to free the memory used by the root array
*/
static void freeRoot(void)
{
	free(pRoot);
	pRoot = NULL;
}

/**
* Retrieves the root value for the given node
* @param v the node id
* @return the root value
*/
static inline int getRoot(int v)
{
	return pRoot[v];
}

/**
* Stores the root value for the node v
* @param v the node id
* @param val the root value
*/
static inline void setRoot(int v, int val)
{
	pRoot[v] = val;
}

/****************THE BSCCS LIST MANAGEMENT*************************************/
/*The temporary storage for the lately found BSCCs' ids*/
static int ** ppNewBSCCs = NULL;

/**
* This function is used to return the number of BSCCs known in _ppBSCCs.
* This number is stored in the 0 element of the ppNewBSCCs array, so it
* is a pointer (int *) which is used to store an integer value only.
* @param _ppBSCCs the structure containing BSCCs along
* with the number of them, which is stored in _ppBSCCs[0]
* @return the number of BSCCs newly found (_ppBSCCs[0])
**/
int getBSCCCount(const int * const * _ppBSCCs) {
        const int * num_tmp = (const int *) (const void *) _ppBSCCs[0];
	return (int) ( (unsigned long int) num_tmp );
}

/**
* This method is used to increase the BSCCs counter in _ppBSCCs[0].
* It increases _ppBSCCs[0] by one.
* @param _ppBSCCs the structure containing BSCCs along
* with the number of them, which is stored in _ppBSCCs[0]
*/
static inline void increaseBSCCCount(int ** _ppBSCCs){
	int* num_tmp = (int *) ( (void *) _ppBSCCs[0] );
	int num_tmp_tmp = (int) ( (unsigned long int) num_tmp );
	_ppBSCCs[0] = (int *) ( (unsigned long int) ( num_tmp_tmp + 1) );
}


/**
* This metod initializes the ppBSCCs array. It allocates the one
* value and stores the zero value there as it is the initial number
* of BSCCs.
*/
static void initBSCCsList(void)
{
        ppNewBSCCs = (int **) calloc((size_t) 1, sizeof(int *));
	/*Store the initial amount of BSCCs (it is zero)*/
	ppNewBSCCs[0] = 0;
}

/**
* This method is used to free memory used for storing BSCCs
* @param ppLastBSCCs the set of BSCCs returned by the
*                   int ** getNewBSCCs(bitset *)
*                   function.
*/
void freeBSCCs(int ** ppLastBSCCs)
{
	if(ppLastBSCCs) {
                int i, num = getBSCCCount((const int **) ppLastBSCCs);
		for( i = 1; i <= num ; i++ ) {
			free(ppLastBSCCs[i]);
		}
		free(ppLastBSCCs);
	}
}

/**
* Returns the found BSCCs list
* @return the list of BSCCS
*/
static int ** getNewBSCCsList(void)
{
	return ppNewBSCCs;
}

/**
* This method adds a new BSCC component.
* Increases the bscc_counter value by 1
*/
static void addBSCCToTheList(void)
{
	int lid;

	increaseBSCCCount(ppNewBSCCs);
        lid = getBSCCCount((const int **) ppNewBSCCs);

	/*Extend the ppBSCCs array length, we say ppNewBSCCs[0]+2 because we
	so far have ppBSCCs[0]+1 elements in it and we need one extra*/
	ppNewBSCCs = (int **) realloc( ppNewBSCCs, ( lid + 1 ) * sizeof( int *) );
	/*Allocate memory to store the id and the number of nodes of the BSCC*/
        ppNewBSCCs[lid] = (int *) calloc((size_t) 2, sizeof(int));
	/*Increase the number of found BSCCs*/
	ppNewBSCCs[ lid ][0] = ++bscc_counter;
	ppNewBSCCs[ lid ][1] = 0;
}

/**
* increases the counter of nodes of the current BSCC
*/
static void addBSCCNode(void)
{
        ppNewBSCCs[getBSCCCount((const int **) ppNewBSCCs)][1]++;
}

/**
* This function checks if the BSCC consist of 1 node and if yes then
* it stores it's id in the ppNewBSCCs[ppNewBSCCs[0]][2] element
*/
static inline void checkForSingleNode(int w)
{
        int lid = getBSCCCount((const int **) ppNewBSCCs);
	if( ppNewBSCCs[ lid ][1] == 1 ) {
		/*Add new element*/
		ppNewBSCCs[ lid ]=(int *)realloc(ppNewBSCCs[lid],3*sizeof(int));
		/*Store value*/
		ppNewBSCCs[lid][2] = w;
	}
}

/**************THE BSCC SEARCH PROCEDURES**************************************/

/**
* This method indicates if the w belongs to the component or not and if
* there exists a path from w to some component or not.
* @param w the node to be tested
* @return TRUE if w is in a component or there is a path from w to some component.
*/
static inline BOOL isInComponent(int w)
{
	return ( get_bit_val(pInComponentStates, w) ?  TRUE : FALSE );
}

/**
* This procedure sets the corresponding value to the v element of the
	*pInComponentStates
* bit set. This is used to mark the component as belonging to some state.
*/
static inline void setInComponent(int v, const BITSET_BLOCK_TYPE bit)
{
	set_bit_val(pInComponentStates, v, bit);
}

/**
* This method retrieves a BSCC from the stack.
* @param *cur_stack the pointer to the stack the BSCC should be read
* @param v the root node of the BSCC
* NOTE: this function uses the global int pointer bscc_stack
*/
static void getBSCCFromStack(int v)
{
	int w;
	addBSCCToTheList();
	do {
		w = popStack(bscc_stack);
		setInComponent(w, BIT_ON);
                if ( pBSCCs[w] == bscc_counter )
                        continue;
		pBSCCs[w] = bscc_counter;
		addBSCCNode();
	}while(w != v);

	checkForSingleNode(w);
}

/**
* This recursive procedure is used to pass through the successive nodes
* in order to detect all the Bottom Strongly Connected Components
* containing the given node. If there is a path from i to some other
* component (SCC) then the procedure exits.
* @param v the initial node
* @param *bscc_stack the stack an eventually found bscc is stored on
* NOTE: this function uses the global int pointer bscc_stack
*/
static void visit_rec(int v)
{
        int initial_root = dfs_order++;

	setRoot(v,initial_root);
	setInComponent(v, BIT_OFF);
	setVisited(v);
	bscc_stack = pushStack(bscc_stack, v);

	/*Iterate through all the successive nodes*/
        mtx_walk_row(pStateSpace, (const int) v, w, dummy) {
		if( isNotVisited(w) ) {
			/*Start recursion*/
			visit_rec(w);
			/*Check if we need to exit search*/
			if( bGlobalSkip ) break;
		}
		if( ! isInComponent(w)) {
			setRoot(v, MIN(getRoot(v), getRoot(w)));
		} else {
			/*There is a way from v to some component to which w belongs*/
			/*So v can not be a part of BSCC, thus skip.*/
			bGlobalSkip = TRUE;
			/*setInComponent(v, BIT_ON);*/

			/*The bscc_stack pointer might have changed due to the stack
			reallocation, therefore we have to reflect this change in the
			stack array here, since it is not going to be changed while
			this recursion any more*/
			stackarray->stackp[0] = bscc_stack;
			break;
		}
	}
        end_mtx_walk_row;
	/*If we did not meet any component yet then it means that there
	can be a BSCC in the stack*/
	if( ! bGlobalSkip ) {
		if ( getRoot(v) == initial_root ) {
			/*Found a BSCC let's get it from the stack.*/
			getBSCCFromStack(v);
		}
	}
}

/**
* This procedure processes the calculations for the nodes in path_stack needed
* for recognizing BSCCs and in case of a BSCC found also gets it from the
* bscc_stack
* @param *path_stack the stack of nodes to do calculations for
*/
static inline void node_calculations(int *path_stack)
{
        int v, v_initial_root;
	v = popStack(path_stack);
	v_initial_root = getRoot(v);
	/*Iterate through all the successive nodes*/
        mtx_walk_row(pStateSpace, (const int) v, w, dummy) {
		if( ! isInComponent(w)) {
			setRoot(v, MIN(getRoot(v), getRoot(w)));
		} else {
			/*There is a way from v to some component to which w belongs*/
			/*So v can not be a part of BSCC, thus skip.*/
			bGlobalSkip = TRUE;
			/*setInComponent(v, BIT_ON);*/
			break;
		}
	}
        end_mtx_walk_row;
	/*If we did not meet any component yet then it means that there
	can be a BSCC in the stack*/
	if( ! bGlobalSkip ) {
		if ( getRoot(v) == v_initial_root ) {
			/*Found a BSCC let's get it from the stack.*/
			getBSCCFromStack(v);
		}
	}
}

/**
* This procedure passes through the successive nodes (depth-first-search)
* in order to detect all the Bottom Strongly Connected Components
* containing the given node. If there is a path from v to some other
* component (SCC) then the procedure exits.
* @param v the initial node
* NOTE: this function uses the global int pointer bscc_stack
*/
static void visit_non_rec(int v)
{
	int w = 0, succ_counter=0, num=0, initial_root = dfs_order;
	int* dfs_stack = stackarray->stackp[1];
	/* Stores tuples of the form (node, number of node's successors) */
	int* path_stack = stackarray->stackp[2];

	/* Push starting node on dfs stack */
	dfs_stack = pushStack(dfs_stack, v);

	/* While dfs stack not empty, process with next element on stack */
	while( ( ( v = popStack(dfs_stack) ) != EMPTY_STACK ) && (!bGlobalSkip)) {
		if( ( w = popStack(path_stack) ) != EMPTY_STACK ) {
			/* Decrease successor counter of predecessor */
			path_stack = pushStack(path_stack, --w);
		}

		initial_root++;
		setRoot(v,initial_root);
		setInComponent(v, BIT_OFF);
		setVisited(v);
		/* Push v on bscc stack */
		bscc_stack = pushStack(bscc_stack, v);
		succ_counter = 0;

		/*Iterate through all the successive nodes*/
                mtx_walk_row(pStateSpace, (const int) v, ww, dummy) {
                        if( isNotVisited(ww) ) {
				/*Push successive node on dfs-stack*/
                                dfs_stack = pushStack(dfs_stack, ww);
				/* A non-visited successor found */
				succ_counter++;
			}
		}
                end_mtx_walk_row;
		/* Push a tuple of node and number of node's successors on the path stack*/
		path_stack = pushStackTuple(path_stack, v, succ_counter);
		num = 0;

		/* As long as all successive nodes have been visited, retrieve elements from
		path stack and do bscc calculations for this nodes */
		while( ((succ_counter = popStack(path_stack)) != EMPTY_STACK ) && (!bGlobalSkip)) {
			/* All successors have been visited */
			if (succ_counter == 0) {
				node_calculations(path_stack);
			} else {
				/*At least one non-visited successor is found,
				process dfs search with next element */
				path_stack = pushStack(path_stack, succ_counter);
				break;
			}
		}
	}
	/*Clean stacks */
	cleanStack(path_stack);
	cleanStack(dfs_stack);

	/*The stack pointers might have changed, therefore we
	have to update the stack array*/
	stackarray->stackp[2] = path_stack;
	stackarray->stackp[1] = dfs_stack;
}


/**
* This function returns the list of ids for a newly found BSCCs
* @param pBSCCsHolder the work structure to store all required data for
*						BSCCs search for a particular state space
* @param pStates the set of states
* @return the two dimensional array (result[][])
* result[1][0] the id of a newly found bsccs,
* result[1][1] the number of elements in a bscc
* result[1][2] the id of the node if result[1][1] = 1
*/
int ** getNewBSCCs(TBSCCs * pBSCCsHolder, const bitset *pStates)
{
	const int size = pBSCCsHolder->size;
	int i, elem;

	/* Get the BSCC method from the runtime.c */
	const int method = get_method_bscc();

	/*Define the function pointer that will be assigned
	to the proper VISIT procedure later*/
	void (*visit)(int) = NULL;

	/*Init the BSCCs list*/
	initBSCCsList();

	/*Allocate the root array;*/
	initRoot(size);

	/*Copy data from TBSCC struct to internal variables*/
	pStateSpace = pBSCCsHolder->pStateSpace;
	pVisitedStates = pBSCCsHolder->pVisitedStates;
	pInComponentStates = pBSCCsHolder->pInComponentStates;
	pBSCCs = pBSCCsHolder->pBSCCs;
	bscc_counter = pBSCCsHolder->bscc_counter;
	dfs_order = pBSCCsHolder->dfs_order;

	/* Choice for visit_rec or visit_non_rec */
	if( method==REC ) {
		printf( "WARNING: Running BSCC search in recursive mode! ");
		printf( "Segmentation fault\n may occur because of insufficient stack size. " );
		printf( "If it does, switch\n to non-recursive mode instead.\n");
		/*Allocate the stackarray*/
		stackarray = getNewStackArray(1);
		visit = visit_rec;
	/* Should be non-recursive otherwise */
	} else {
		/*Allocate the stackarray*/
		stackarray = getNewStackArray(3);
		visit = visit_non_rec;
	}
	bscc_stack = stackarray->stackp[0];

        /*Pass through all the states from the *pStates.*/
        /* get_idx_next_non_zero() is more efficient than checking every bit in
           pStates individually. David N. Jansen. */
        i = state_index_NONE;
        while ( (i = get_idx_next_non_zero(pStates, i)) != state_index_NONE ) {
			if( isNotVisited(i) ) {
				bGlobalSkip = FALSE;

				visit(i);

                                /*Mark all the components in stack as belonging
				to some component*/
				while( ( elem = popStack(bscc_stack) ) != EMPTY_STACK ){
					setInComponent(elem, BIT_ON);
				}
			}
	}

	/*The bscc stack pointer might have changed, therefore we
	have to update the stack array*/
	stackarray->stackp[0] = bscc_stack;

	/*Copy stack variable values back to TBSCCs*/
	pBSCCsHolder->bscc_counter = bscc_counter;
	pBSCCsHolder->dfs_order = dfs_order;

	/*Free the stack memory*/
	freeStackArray(stackarray);
	stackarray = NULL;
	/*The bscc_stack gets freed when the stackarray is freed*/
	bscc_stack = NULL;
	freeRoot();

	return getNewBSCCsList();
}

/**
* This function searches for BSCCs that contain good states (pGoodStates)
* @param pStateSpace the sparse matrix of the CTMC
* @param pGoodStates the bitset containing the good states
* @param pppNonTrivBSCCBitSets (a return parameter) the pointer to an array of bitsets
*				that will contain non-trivial BSCCs with good states (pGoodStates)
* @param pNumberOfNonTrivBSCCs (a return parameter) the pointer to the number of non-trivial
*				BSCCs containing good states (pGoodStates), i.e. the size of
*				(* pppNonTrivBSCCBitSets) array
* @return the bitset containing the good states (pGoodStates) each of which is a trivial BSCCs
*/
bitset * getGoodStateBSCCs(const sparse * pStateSpace_local,
                                                        const bitset *
                                                        pGoodStates,
							bitset *** pppNonTrivBSCCBitSets, int * pNumberOfNonTrivBSCCs ){
	int i, index;
	bitset * pTmpBitSet = NULL;
	/*Create the initial structure for storing BSCC search data*/
        TBSCCs * pBSCCsHolder = allocateTBSCCs(pStateSpace_local);
        int ** ppNewBSCCs_local;
        const
	int * bscc_mapping;

        /* The structure that stores the set of states belonging to BSCC i */
	bitset ** ppBSCCBitSets;
	/* The bitset indicating the non-trivial BSCCs containing Psi states */
	bitset * pAllowedNonTrivBSCCsBitSet;
	/* The bitset indicating the trivial BSCCs containing Psi states */
	bitset * pAllowedTrivialBSCCsBitSet;
	/* The bitset indicating the Psi states, that form a trivial BSCC */
	bitset * pTrivialBSCCBitSet = NULL;
	/* The structure, which stores the set of states belonging to non-trivial BSCC i */
	bitset ** ppNonTrivBSCCBitSets = NULL;
	/* The state space dimension */
        const int n_states = mtx_rows(pStateSpace_local);
	/* The number of BSCCs found */
	int numberOfBSCCs;
	/* The number of non-trivial BSCCs containing good states */
	int numberOfNonTrivBSCCs = 0;

	/* Find BSCCs starting in Psi states */
	/* NOTE: We might find BSCCs that do not contain Psi states !!!*/
	printf("Find new BSCCs....\n");
        ppNewBSCCs_local = getNewBSCCs(pBSCCsHolder, pGoodStates);

	numberOfBSCCs = getBSCCCounter( pBSCCsHolder );
	pAllowedTrivialBSCCsBitSet = get_new_bitset(numberOfBSCCs);
	pAllowedNonTrivBSCCsBitSet = get_new_bitset(numberOfBSCCs);

	/* Initialize the structure, that stores the set of states belonging to BSCC with array index i */
        ppBSCCBitSets = (bitset **) calloc((size_t) numberOfBSCCs,
                        sizeof(bitset *));
	for( i = 0; i < numberOfBSCCs; i++ ){
		ppBSCCBitSets[i] = get_new_bitset(n_states);
	}

	bscc_mapping = getStatesBSCCsMapping(pBSCCsHolder);

	/* Initialize ppBSCCBitSets by sorting out the states belonging to the same BSCC */
	for( i = 0; i < n_states; i++ ){
		/* Does state i belong to a BSCC? */
		if(  bscc_mapping[i] != 0 ){
			set_bit_val(ppBSCCBitSets[bscc_mapping[i] -1], i, BIT_ON);
		}
	}

	/* Compute the BSCCs containing Psi states */
	for(i = 0; i < numberOfBSCCs; i++){
		/* Check if the BSCC contains at least one Psi state */
		pTmpBitSet = and( ppBSCCBitSets[i], pGoodStates );
		/* If the bitset contains elements, i.e. there are Psi states in the given BSCC */
		if( !is_bitset_zero( pTmpBitSet ) ){
			/* Is the found BSCC containing Psi states trivial? */
			if( count_non_zero( ppBSCCBitSets[i] ) > 1 ){
				numberOfNonTrivBSCCs++;
				set_bit_val( pAllowedNonTrivBSCCsBitSet, i, BIT_ON );
			}else{
				set_bit_val( pAllowedTrivialBSCCsBitSet, i, BIT_ON);
			}
		}else{
			/* If there are no Psi states than there is no reason of keeping this BSCC bitset */
			free_bitset( ppBSCCBitSets[i] );
			/* WARNING: We should set it to NULL in order to prevent double deallocation later */
			ppBSCCBitSets[i] = NULL;
		}
		free_bitset( pTmpBitSet ); pTmpBitSet = NULL;
	}

	/* Allocate memory for the structures containing the return values in the end */
	pTrivialBSCCBitSet = get_new_bitset( n_states );
        ppNonTrivBSCCBitSets = (bitset **) calloc((size_t) numberOfNonTrivBSCCs,
                        sizeof(bitset *));

	/* Part all BSCCs that contain good states into trivial and non-trivial */
	index = -1; i = 0;
	while( ( index = get_idx_next_non_zero( pAllowedNonTrivBSCCsBitSet, index ) ) != -1 ){
		/* For every non-trivial BSCC containing good states copy the bitset pointer */
		ppNonTrivBSCCBitSets[i] = ppBSCCBitSets[index];
		/* WARNING: Remove the pointer to this BSCC because it will be used in the caller function */
		ppBSCCBitSets[index] = NULL;
		i++;
	}

	/* Manage the trivial BSCCs that contain Good states */
	index = -1;
	while( ( index = get_idx_next_non_zero( pAllowedTrivialBSCCsBitSet, index ) ) != -1 ){
		/* Remember the good states forming trivial BSCCs */
                i = get_idx_next_non_zero(ppBSCCBitSets[index],
				state_index_NONE);
		set_bit_val(pTrivialBSCCBitSet, i, BIT_ON);
		/* WARNING: Free the BSCC right away, to avoid doing it later */
		free_bitset( ppBSCCBitSets[index] ); ppBSCCBitSets[index] = NULL;
	}

	/* Free memory */
        freeBSCCs(ppNewBSCCs_local); ppNewBSCCs_local = NULL;
	freeTBSCC( pBSCCsHolder ); pBSCCsHolder = NULL;

	free_bitset( pAllowedNonTrivBSCCsBitSet ); pAllowedNonTrivBSCCsBitSet = NULL;
	free_bitset( pAllowedTrivialBSCCsBitSet ); pAllowedTrivialBSCCsBitSet = NULL;

	/* WARNING: We do not have to clean the BSCC sets of this array because either */
	/* it has been done before or we use them to be returned in pppBSCCs */
	free( ppBSCCBitSets ); ppBSCCBitSets = NULL;

	/* Assign the return values */
	( * pNumberOfNonTrivBSCCs ) = numberOfNonTrivBSCCs;
	( * pppNonTrivBSCCBitSets ) =  ppNonTrivBSCCBitSets;

	return pTrivialBSCCBitSet;
}
