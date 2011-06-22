/*-----------------------------------------------------------------
/
/   pga.c - a simple binary GA package
/
c  Author:
c
c  Richard Judson
c  rsjuds@ca.sandia.gov  
c  (510)294-1438           
c  FAX (510)294-2234
c
c  Center for Computational Engineering
c  Sandia National Laboratories
c  P.O. Box 969
c  Livermore, CA 94551-0969
c
c-----------------------------------------------------------------------
c
c                       Disclaimer
c
c  We don't claim this code is good for anything - if you think it is,
c  great, but it's up to you to decide. If the code doesn't work: tough.
c  If you lose a million because the code messes up, it's you that's out
c  a million, not us. If you don't like this disclaimer: tough. We
c  reserve the right to do the absolute minimum provided by law, up to
c  and including nothing.
c
c-----------------------------------------------------------------------
c
c
c  Copyright 1995, Sandia Corporation. The U.S. Government retains
c  a limited license in this software.
c
c  The U.S. Government retains, in this software, a paid-up,
c  nonexclusive, irrevocable worldwide license to reproduce, prepare
c  derivative works, perform publicly and display publicly by or for the
c  Government, including the right to distribute to other Government
c  contractors. 
c
c  Neither the United States, the U.S. Dept. of Energy, nor any of their
c  employees, makes any warranty, express or implied, or assumes any
c  legal liability or responsibility for the accuracy, completeness, or
c  usefulness of any information, apparatus, product, or product
c  disclosed, or represents that its use would not infringe privately
c  owned rights. 
c
c------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "pga.h"	
#include "eval.h"	


/*------------------------------------------------------------------*/
main()
{
  GA_main();
}

/*------------------------------------------------------------------
/
/   GA_main
/
/------------------------------------------------------------------*/
void GA_main()
{
  int         igen, iniche, jniche;
  population  *pop, newpop;
  bestpop     *best, newbest;
  parm_struct parm;
  int         i;
  FILE *fpout;


  /* Read in the input parameters.
  ----------------------------- */
  read_GAparms(&parm);
  fpout = fopen("GAout.dat","w");
  write_GAparms(&parm,fpout);
  write_GAparms(&parm,stdout);

  /* Allocate space for the populations 
  ---------------------------------------*/
  pop  = (population *)malloc(parm.niches*sizeof(population));
  best = (bestpop *)malloc(parm.niches*sizeof(bestpop));

  /* Perform the initialization
  ------------------------------------------------------------ */
  GA_setup(&parm,pop,&newpop,best,&newbest,fpout);

  /* Perform the requested number of generations of GA analysis.
  ------------------------------------------------------------ */
  for(igen=1;igen<=parm.ngen;igen++) {
    GA_loop(&parm,pop,&newpop,best,&newbest,fpout,igen);
  }
  /* End of loop over generations
     -------------------------------------*/
  fclose(fpout);
  free(pop);
  free(best);
}

/*------------------------------------------------------------------
/
/   GA_loop
/
/------------------------------------------------------------------*/
void GA_loop(parm_struct *parm, population *pop, population* newpop,
             bestpop *best, bestpop* newbest, FILE *fpout, int igen)
{
  int iniche,jniche,i;
  double param_list[MAXWORDS];

  /* Evaluate the fitnesses.
  ------------------------ */
  for(iniche=0;iniche<parm->niches;iniche++)  {
    for (i=0; i<pop->npop; i++) {
      filter(parm, pop[iniche].indiv[i].chromosome, param_list);
      pop[iniche].indiv[i].fitness = 
	eval_user(iniche,parm->eval,parm->words,param_list);
    }
  }
  
  /* If parm->ni_gen[0] is flagged then a sharing of best values should
     occur before any crossovers. Note this should only be done on
     restarted runs.
  ------------------------------------------------------------- */
  if (parm->restart==1 && parm->niches>1 && parm->ni_gen[0]==0  && igen==1) {
    for(iniche=0;iniche<parm->niches;iniche++) 
      for(jniche=0;jniche<parm->niches;jniche++) {
	pop[iniche].indiv[parm->npop-1-jniche]=pop[jniche].indiv[0];
      }      
  }
  
  /* Perform the selection.
  ----------------------- */
  for(iniche=0;iniche<parm->niches;iniche++) {
    if(parm->select_method==ROULETTE) 
      select_pop_roulette(parm, &pop[iniche]);
    else if(parm->select_method==STEP_FTN) 
      select_pop_stepftn(parm, &pop[iniche]);
    else 
      fprintf(stdout,"Illigal selection flag: %d\n",parm->select_method);
  }
  
  /* Gather statistics and update the list of best individuals.
  ----------------------------------------------------------- */
  for(iniche=0;iniche<parm->niches;iniche++)  {
    best_pop(&pop[iniche], &best[iniche], newbest);
    best[iniche] = (*newbest);
    if(igen%parm->save_freq ==0) {
      statistics_pop(igen, iniche, parm, &pop[iniche], &best[iniche],fpout);
      statistics_pop(igen, iniche, parm, &pop[iniche], &best[iniche],stdout);
      print_pop(iniche,parm, &pop[iniche]);
      print_best(iniche,parm, &best[iniche]);
    }
  }  
  
  /* Pass good individuals across niches if the generation is selected.
  ------------------------------------------------------------------- */
  if (parm->niches > 1) 
    for (i=0; i<parm->niflag; i++) 
      if (igen==parm->ni_gen[i])  {
	fprintf(stdout,"\nAt generation %d share best values\n\n",igen);
	fprintf(fpout,"\nAt generation %d share best values\n\n",igen);
	for(iniche=0;iniche<parm->niches;iniche++) 
	  for(jniche=0;jniche<parm->niches;jniche++) {
	    pop[iniche].indiv[parm->npop-1-jniche]=pop[jniche].indiv[0];
	    if(parm->select_method==ROULETTE) 
	      select_pop_roulette(parm, &pop[iniche]) ;
	    else                   
	      select_pop_stepftn(parm, &pop[iniche]) ;
	  }      
      }  
  
  /* Perform crossover.
  -------------------- */
  for(iniche=0;iniche<parm->niches;iniche++)  {
    crossover_pop(parm, &pop[iniche], newpop,igen);
    pop[iniche] = (*newpop);
  }
  
  /* Perform mutation.
  ------------------- */
  for(iniche=0;iniche<parm->niches;iniche++) 
    mutate_pop(parm, &pop[iniche]);
}
/*------------------------------------------------------------------
/
/   GA_setup
/
/------------------------------------------------------------------*/
void GA_setup(parm_struct *parm, population *pop, population* newpop,
	      bestpop *best, bestpop* newbest, FILE *fpout)
{
  int iniche;


   /* If indicated share the best values among the niches at the
      beginning of a restarted run.
   ------------------------------------------------------------ */
  if (parm->restart==1 && 
      parm->niches>1 && 
      parm->niflag!=0 && 
      parm->ni_gen[0]==0) {
    fprintf(stdout,"At beginning of restart share best values\n\n");
    fprintf(fpout,"\nAt beginning of restart share best values\n\n");
  }

  for(iniche=0;iniche<parm->niches;iniche++) {
    best[iniche].nbest = parm->nbest;
  }

  GAran(parm->MSEED);

   /* Create the initial population.
   ------------------------------ */
  if (!parm->restart)
    for(iniche=0;iniche<parm->niches;iniche++) {
      pop[iniche].npop=parm->npop;
      initialize_pop(parm,&pop[iniche],iniche);
    }
  else {

    /* If this is a restarted run read from GApop.dat files.
    ------------------------------------------------------ */
    for(iniche=0;iniche<parm->niches;iniche++) {
      initialize_pop(parm,&pop[iniche],iniche);
    }
  }
      
  /* Create the first "best" population.
  ------------------------------------- */   
  if (!parm->restart) 
    for(iniche=0;iniche<parm->niches;iniche++)  
      initialize_best(parm,&best[iniche],iniche);

  else {

    /* If this is a restarted run read from GAbest.dat files.
    -------------------------------------------------------- */
    for(iniche=0;iniche<parm->niches;iniche++) {
      initialize_best(parm,&best[iniche],iniche);
    }
  }

  /* Define newpop and newbest data structures.
  ------------------------------------------- */  
  initialize_pop(parm,newpop,-1);
  initialize_best(parm,newbest,-1);
}
