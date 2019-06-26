/* This function sets structure for individual choice data. The main 
 * output is the set of choice-agent index pairs. Irrelevant pairs are 
 * excluded (those, in which the choice is not in the agent's choice set). 
 * Pairs are sorted lexicographically by agent and inverse preference 
 * order (that is, most preferred option comes last). */

#include <stdint.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	/*
	* RHS[0] = choice set matrix, num_choices X num_agents
	* RHS[1] = lists of ranked preferences, list_length X num_agents
	*/
	size_t num_agents, i, j, nnz, last;
    uint64_t *jc, *agent_idx;
	uint16_t num_choices, preflist_length, nnz_j, *preflist, 
		*num_skipped, *num_listed, *choice_idx, *tail;
	mxLogical *choice_set, *choice_set_first, *choice_set_j;
	
    /* The matrix of eligibility dummies. Column -- agent, row -- choice */
	choice_set = (mxLogical *)mxGetData(prhs[0]);
	choice_set_first = choice_set;
	num_choices = mxGetM(prhs[0]);
	num_agents = mxGetN(prhs[0]);

	preflist = (uint16_t *)mxGetData(prhs[1]);
	preflist_length = mxGetM(prhs[1]);
	
	/* Output */
	plhs[2] = mxCreateNumericMatrix(num_agents, 1, mxUINT16_CLASS, mxREAL);
	num_skipped = (uint16_t *)mxGetData(plhs[1]);
	
	plhs[3] = mxCreateNumericMatrix(num_agents, 1, mxUINT16_CLASS, mxREAL);
	num_listed = (uint16_t *)mxGetData(plhs[2]);
	
	/* Find how much memory we have to allocate for non-zero elements 
	   and compute the choice set size for all agents */
	nnz = 0;
	
	for(j=0; j<num_agents; j++) {
		
		nnz_j = 0;
		for(i=0; i<num_choices; i++){
			nnz_j += *choice_set++;
		}
		nnz += nnz_j;
		num_skipped[j] = nnz_j;
		
	}
	choice_set = choice_set_first;
	
    /* Choice index in the sparse format. Similar to IR in the standard 
     * MATLAB sparse array structure. Using MATLAB indexing. */
	plhs[0] = mxCreateNumericMatrix(nnz, 1, mxUINT16_CLASS, mxREAL);
	choice_idx = (uint16_t *)mxGetData(plhs[0]);
    
    plhs[1] = mxCreateNumericMatrix(num_agents+1, 1, mxUINT64_CLASS, mxREAL);
    jc = (uint64_t *)mxGetData(plhs[1]);
    /* Use MATLAB indexing (start from 1) */
    jc[0] = 1;
	
    /* Workspace for tracking choice sets of individual agents */
    choice_set_j = (mxLogical *)malloc(num_choices*sizeof(mxLogical));
    
	for(j=0; j<num_agents; j++) {
		
        /* We will need to modify j's choice set. Copy it to avoid
         * modifying the input argument */
        memcpy(choice_set_j, choice_set, sizeof(mxLogical)*num_choices);
        
		/* First, cut out preference lists */
		tail = choice_idx + num_skipped[j] - 1;
		for (i=0; i<preflist_length; i++){
			/* Zeros code skipped/invalid items */
			if(*preflist){
				
				/* Is agent j eligible for the next entry? If some choice 
                 * is double-listed, this drops all lower-ranked entries 
                 * with the same choice */
				if(choice_set_j[*preflist - 1]){
					/* If so, tick off this entry's dummy */
					choice_set_j[*preflist - 1] = 0;
					/* Add the entry to the tail, move the tail pointer up */
					*tail = *preflist;
					tail--;
				}
			}
			/* Next list item */
			preflist++;
		}
		/* Number of listed items = how far up we moved the tail pointer */
		num_listed[j] = choice_idx + num_skipped[j] - 1 - tail;
		
		/* Second, process the remaining choice set items. The order 
         * doesn't matter here */
		for (i=0; i<num_choices; i++){
			if(*choice_set_j++){
				*choice_idx = i+1;
				choice_idx++;
			}
		}
		/* Jump over the tail */
		choice_idx += num_listed[j];
        /* Save index information */
		num_skipped[j] -= num_listed[j];
        jc[j+1] = jc[j] + num_skipped[j] + num_listed[j];  
        
        /* Go to the next agent's choice set */
        choice_set += num_choices;
	}
    
    
    /* If needed, produce a full-size agent index */
    if(nlhs>4) {
        plhs[4] = mxCreateNumericMatrix(nnz, 1, mxUINT64_CLASS, mxREAL);
        agent_idx = (uint64_t *)mxGetData(plhs[4]);
        
        for(j=0; j<num_agents; j++) {
            last = jc[j+1] - 1;
            for(i=0; i<last; i++){
                *agent_idx = j+1;
                agent_idx++;
            }
        }
    }
    
    free(choice_set_j);
    
	return;
	
}
