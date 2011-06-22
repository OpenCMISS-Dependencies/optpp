//------------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_STREAM
#define WANT_MATH

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <malloc.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "opt.h"
#include "cblas.h"
#include "precisio.h"

#include "pga.h"	
#include "eval.h"	
#include "optga.h"

void GA_loop(NLP0 *nlp, parm_struct *parm, population*,population*,bestpop*,
	     bestpop*,FILE*,int);

void OptGA::InitOpt()
{
  nlp->initFcn();

  if (debug_)
    nlp->setDebug;

  ret_code = 0;

  if(nlp->hasConstraints()){
    CompoundConstraint* constraints = nlp->getConstraints();
    ColumnVector xstart = nlp->getXc();
    double feas_tol = tol.getCTol();
    bool feasible = constraints->amIFeasible(xstart, feas_tol);
    if (!feasible) {
      *optout << "ga WARNING:  Initial guess not feasible.\n"
	      << "GA may be unable to make progress." << endl;
    }
  }

  fprintf(outfile,"\n\t\t%s Method\n", method);
  fprintf(outfile,"==============================================================================\n");
  fprintf(outfile,"  Generation  Niche     Best           Avg.        avg.best        conv   lost\n");
  fprintf(outfile,"------------------------------------------------------------------------------\n");
  fflush(outfile);
}
//
// This is the main routine for the Genetic Algorithm method
//
void OptGA::Optimize()
{

// Allocate local vectors 

  int ndim = dim;
  ColumnVector x(ndim);

  int         igen;
  population  *pop, newpop;
  bestpop     *best, newbest;

  FILE *fpout = outfile;

  parm_struct *parm;

// Initialize the GA

  InitOpt();

  if (ret_code == 0) {

    /* Allocate space for the populations 
       ---------------------------------------*/

    parm = ga_parms;
    pop  = (population *)malloc(parm->niches*sizeof(population));
    best = (bestpop *)malloc(parm->niches*sizeof(bestpop));

    /* Perform the initialization
       ------------------------------------------------------------ */
    GA_setup(parm,pop,&newpop,best,&newbest,fpout);

    /* Perform the requested number of generations of GA analysis.
       ------------------------------------------------------------ */
    for(igen=1;igen<=parm->ngen;igen++) {
      GA_loop(nlp,parm,pop,&newpop,best,&newbest,fpout,igen);
    }
    /* End of loop over generations
       -------------------------------------*/
    iter_taken = parm->ngen;
    fcn_evals  = (parm->ngen)*(parm->niches)*(parm->npop);

    free(pop);
    free(best);
  }
}

void OptGA::PrintStatus(char *s) // Set Message
{
  fprintf(outfile,"\nOptGA: %s\n",s);
  fprintf(outfile,"Optimization method            = %s\n",method);
  fprintf(outfile,"Dimension of the problem       = %d\n",dim);
  fprintf(outfile,"Return code                    = %d\n",ret_code);
  fprintf(outfile,"%s\n",mesg);
  fprintf(outfile,"Number of iterations taken     = %d\n",iter_taken);
  fprintf(outfile,"Number of function evaluations = %d\n",fcn_evals);

  nlp->PrintState(s);
  tol.PrintTol();

  write_GAparms(ga_parms,outfile);

}

int OptGA::CheckConvg() // Check convergence
{
  int    n;
  double stol, ftol, rftol;
  double xnorm, snorm;
  ColumnVector xc;

  double step_tol, fvalue;

  n  = nlp->GetDim();
  xc = nlp->GetXc();
  fvalue = nlp->GetF();

  xnorm =  Norm2(xc);
  
// Test 1. step tolerance 

  ColumnVector step(n);
  step = xc - xprev;
  step_tol = tol.GetStepTol();
  snorm = Norm2(step);
  stol  = step_tol*max(1.0,xnorm);
  if (snorm  <= stol) {
    strcpy(mesg,"Algorithm converged - Norm of last step is less than step tolerance");
    fprintf(outfile,"snorm = %12.4e, stol = %12.4e\n", snorm, stol);
    return 1;
  }
  
// Test 2. function tolerance
  Real deltaf = fprev - fvalue;
  ftol = tol.GetFTol();
  rftol = ftol*max(1.0,fabs(fvalue));
  if (deltaf <= rftol) {
    strcpy(mesg,"Algorithm converged - Difference in successive fcn values is less than fcn tolerance");
    fprintf(outfile,"deltaf = %12.4e, ftol = %12.4e\n", deltaf, ftol);
    return 2;
  }
  

  // Nothing to report 

  return 0;

}

/*------------------------------------------------------------------
/
/   GA_loop
/
/------------------------------------------------------------------*/
void GA_loop(NLP0 *nlp, parm_struct *parm, population *pop, population* newpop,
             bestpop *best, bestpop* newbest, FILE *fpout, int igen)
{
  int iniche,jniche,i,k;
  double param_list[MAXWORDS];

  int ndim = nlp->GetDim();
  ColumnVector xparm(ndim);
  int best_niche;
  double best_fit;

  /* Evaluate the fitnesses.
  ------------------------ */
  for(iniche=0;iniche<parm->niches;iniche++)  {
    for (i=0; i<pop->npop; i++) {
      filter(parm, pop[iniche].indiv[i].chromosome, param_list);
      for (k = 0; k < ndim; k++) xparm(k + 1) = param_list[k];
      pop[iniche].indiv[i].fitness = nlp->EvalF(xparm);
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
      fprintf(stdout,"Illegal selection flag: %d\n",parm->select_method);
  }
  
  /* Gather statistics and update the list of best individuals.
  ----------------------------------------------------------- */
  best_fit = 1.e6;
  for(iniche=0;iniche<parm->niches;iniche++)  {
    best_pop(&pop[iniche], &best[iniche], newbest);
    best[iniche] = (*newbest);
    if(igen%parm->save_freq ==0) {
      statistics_pop(igen, iniche, parm, &pop[iniche], &best[iniche],fpout);
//      statistics_pop(igen, iniche, parm, &pop[iniche], &best[iniche],stdout);
      print_pop(iniche,parm, &pop[iniche]);
      print_best(iniche,parm, &best[iniche]);
//
//    Find and store best individual
//
      if(best[iniche].indiv[0].fitness < best_fit) {
	best_fit = best[iniche].indiv[0].fitness;
	best_niche = iniche;
	filter(parm, best[iniche].indiv[0].chromosome, param_list);
	for (k = 0; k < ndim; k++) xparm(k + 1) = param_list[k];
	nlp->SetX(xparm);
	nlp->SetF(best_fit);
      }
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

//-------------------------------------------------------------------------
// GANLF0 method routines.
//-------------------------------------------------------------------------
void GANLF0::InitFcn() // Initialize Function
{
  if (init_flag == NO)
  {
      init_fcn(dim, xc);
      init_flag = YES;
  }
  else
  {
    fprintf(outfile,"GANLF0:InitFcn: Warning - initialization called twice\n");
    init_fcn(dim, xc);
  }
}
double GANLF0::EvalF() // Evaluate Function
{
  fcn(dim, xc, fvalue);
  nfevals++;

  return fvalue;
}
double GANLF0::EvalF(ColumnVector& x) // Evaluate Function at x
{
  double fx;
  fcn(dim, x, fx);
  nfevals++;
  
  return fx;
}
