/*-----------------------------------------------------------------
/
/   pga_util.c - utility routines for pga
/
/------------------------------------------------------------------*/
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

/*------------------------------------------------------------------
/
/   read_GAparms
/
/------------------------------------------------------------------*/
void read_GAparms(parm_struct* parm, FILE *fp)
{
  char line[120];
  int istop, i, start, iniche;
  
  fgets(line, 120, fp);
  sscanf(line,"%s",parm->eval);
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->restart);
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->niches);
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->npop);
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->ngen);
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->select_method);
  
  /* Determine at which, if any, generations should best arrays 
     should be communicated among the various niches.
    ----------------------------------------------------------- */
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->niflag);
  if (parm->niflag>MAX_NI)  {
    fprintf(stdout,"ERROR: niflag=%d is greater than MAX_NI=%d\n",
	    parm->niflag,MAX_NI);
    exit(1) ;
  }
  if (parm->niflag>0) {
    for (i=0; i<parm->niflag; i++) {
      fgets(line,80,fp);
      sscanf(line,"%d",&parm->ni_gen[i]) ;
      if (parm->ni_gen[i]>parm->ngen) {
	fprintf(stdout,"warning: niche interaction assigned for igen>ngen\n");
	fprintf(stdout,"niche interaction set to ngen\n");
	parm->ni_gen[i]=parm->ngen;
      }
    }
  }
  if (parm->niches==1 && parm->niflag>0)  {
    fprintf(stdout,"Warning: Only 1 niche, niche interaction turned off\n");
    parm->niflag =0;
  }
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->print_flag);
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->nbest);
  if (parm->npop <parm->nbest) {
    fprintf(stdout,"Warning: the population size=%d is less than Nbest=%d\n",
	    parm->npop,parm->nbest) ;
    fprintf(stdout,"Nbest set to the population size\n");
    parm->nbest=parm->npop;
  }
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->MSEED);
  
  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->save_freq);
  
  fgets(line, 120, fp);
  sscanf(line,"%lf",&parm->select_rate);
  if (parm->select_rate*(double) (parm->npop) < 2.0) {
    fprintf(stdout,"ERROR: population size=%d is too small for the select_rate\n",
	    parm->npop) ;
    exit(1);
  }
    
  fgets(line, 120, fp);
  sscanf(line,"%lf",&parm->mutate_rate);
  
  fgets(line, 120, fp);
  sscanf(line,"%lf",&parm->crossover_rate);

  fgets(line, 120, fp);
  sscanf(line,"%lf",&parm->threshold);

  fgets(line, 120, fp);
  sscanf(line,"%d",&parm->words);
  
  start = 0;
  for(i=0;i<parm->words;i++) {
    fgets(line, 120, fp);
    sscanf(line,"%d%d%lf%lf",&parm->length[i],
	   &parm->type[i],
	   &parm->min[i],
	   &parm->max[i]);
    parm->start[i] = start;
    start += parm->length[i];
  }
  
  fprintf(stdout,"It is suggested that mutate_rate be approx.= %12.5f\n",
	 (1.0 / (double) (parm->words)) ) ;
  
  parm->chr_length = start;
  fclose(fp);
  
  istop =0;
  if( parm->chr_length > MAXLOCI) {
    fprintf(stdout," chr_length > MAXLOCI\n");
    istop=1;
  }
  
  if( parm->npop > MAXPOP) {
    fprintf(stdout," npop > MAXLPOP\n");
    istop=1;
  }
  
  if( parm->nbest > MAXBEST) {
    fprintf(stdout," nbest > MAXBEST\n");
    istop=1;
  }
  
  if( parm->words > MAXWORDS) {
    fprintf(stdout," words > MAXWORDS\n");
    istop=1;
  }
  
  if(istop==1) exit(0);
}
/*------------------------------------------------------------------
/
/   write_GAparms
/
/------------------------------------------------------------------*/
void write_GAparms(parm_struct* parm,FILE* fp)
{
  int i;

  fprintf(fp,"Evaluation function:    %s\n",parm->eval);
  fprintf(fp,"Restart flag    :       %5d\n",parm->restart);
  fprintf(fp,"Number of niches:       %5d\n",parm->niches);
  fprintf(fp,"Niche interaction flag  %5d\n",parm->niflag);
  
  for (i=0; i<parm->niflag; i++) {
    fprintf(fp,"\nShare bests at generation = %d\n",parm->ni_gen[i]) ;
  }
  
  fprintf(fp,"Print statistics every %d generation\n",parm->print_flag);
  
  fprintf(fp,"Population size:        %5d\n",parm->npop);
  fprintf(fp,"Number of generations:  %5d\n",parm->ngen);
  fprintf(fp,"Nbest:                  %5d\n",parm->nbest);
  fprintf(fp,"Length of chromosome:   %5d\n",parm->chr_length);
  fprintf(fp,"MSEED:                  %12d\n",parm->MSEED);
  fprintf(fp,"Save frequency:         %5d\n",parm->save_freq);
  fprintf(fp,"Selection rate:         %10.6f\n",parm->select_rate);
  fprintf(fp,"Mutation rate:          %10.6f\n",parm->mutate_rate);
  fprintf(fp,"Convergence threshold:  %10.6f\n",parm->threshold);
  fprintf(fp,"Crossover rate:         %10.6f\n",parm->crossover_rate);
  fprintf(fp,"Number of words:        %5d\n",parm->words);
  fprintf(fp,"Select Method flag:     %5d\n",parm->select_method);
  fprintf(fp,"====================================================\n");
  fprintf(fp," word  start  length    type      min       max\n");
  fprintf(fp,"----------------------------------------------------\n");
  for(i=0;i<parm->words;i++) {
    fprintf(fp,"  %3d  %4d     %3d     %2d    %8.2f  %8.2f\n",
	    i, 
	    parm->start[i],
	    parm->length[i],
	    parm->type[i],
	    parm->min[i],
	    parm->max[i]);
  }
  fprintf(fp,"==============================================================================\n");
  fprintf(fp,"  Generation  Niche     Best           Avg.        avg.best        conv   lost\n");
  fprintf(fp,"------------------------------------------------------------------------------\n");
  fflush(fp);
}
/*------------------------------------------------------------------
/
/   print_best
/
/------------------------------------------------------------------*/
void print_best(int iniche,parm_struct* parm,bestpop*  best)
{
  FILE *fp;
  int i, j, k;
  double   array[MAXWORDS];
  char name[80];

  make_name("GAbest",iniche,".dat",name);
  
  if((fp = fopen(name,"w"))==NULL) {
    fprintf(stdout,"Cannot open %s in print_best()\n",name);
    return;
  }
  
  for(i=0;i<best->nbest;i++) {
    filter(parm, best->indiv[i].chromosome, array);
    
    fprintf(fp,"-----------------------------------------------\n");
    fprintf(fp,"Individual: %5d    Fitness: %12.6f \n",
	    i, best->indiv[i].fitness);
    for(j=0;j<parm->words;j++) {
      for(k=parm->start[j];k<parm->start[j]+parm->length[j];k++) 
	fprintf(fp,"%1d",best->indiv[i].chromosome[k]);
      
      for(k=0;k<25-parm->length[j];k++) 
	fprintf(fp," ");
      
      fprintf(fp,"%10.4f\n",array[j]);
    }
  }
  fclose(fp);

  make_name("GAarray",iniche,".dat",name);

  if((fp = fopen(name,"w"))==NULL) {
    fprintf(stdout,"Cannot open %s in print_best()\n",name);
    return;
  }
  
  for(i=0;i<best->nbest;i++) {
    filter(parm, best->indiv[i].chromosome, array);
    fprintf(fp,"%10.4e\t",best->indiv[i].fitness);
    for(j=0;j<parm->words;j++) 
      fprintf(fp,"%15.8e ",array[j]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}
/*------------------------------------------------------------------
/
/   print_pop
/
/------------------------------------------------------------------*/
void print_pop(int iniche,parm_struct* parm,population* pop)
{
  FILE *fp;
  int i, j;
  char name[80];

  make_name("GApop",iniche,".dat",name);
  if((fp = fopen(name,"w"))==NULL) {
    fprintf(stdout,"Cannot open %s in print_pop()\n",name);
    return;
  }
  
  for(i=0;i<pop->npop;i++) {
    fprintf(fp,"%5d Fitness: %12.6f ",i, pop->indiv[i].fitness);
    for(j=0;j<parm->chr_length;j++) 
      fprintf(fp,"%1d",pop->indiv[i].chromosome[j]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}
/*------------------------------------------------------------------
/
/   statistics_pop
/
/------------------------------------------------------------------*/
void statistics_pop(int igen,int iniche,parm_struct* parm,
		    population* pop,bestpop* best,FILE* fp)
{
  double bestfit, avgfit, avg_best,conval;
  int i,j, counter, nconv, nlost;
  
  avgfit = 0.0;
  for(i=0;i<pop->npop;i++)
    avgfit += pop->indiv[i].fitness;
  
  /* Calculate the average of the best fitnesses.
  --------------------------------------------- */   
  avg_best = 0.0;
  for(i=0;i<parm->nbest;i++)  
    avg_best += best->indiv[i].fitness;
  avg_best /= parm->nbest;
  
  nconv=nlost=0;
  for(j=0;j<parm->chr_length;j++) {
    counter=0;
    for(i=0;i<pop->npop;i++) 
      counter += pop->indiv[i].chromosome[j];
    conval = (double)counter/(double)pop->npop;
    if(conval>parm->threshold || conval<(1.0-parm->threshold))
      nconv++;
    if(counter==0 || counter==pop->npop)
      nlost++;
  }
  
  avgfit /= pop->npop;
  bestfit = best->indiv[0].fitness;
  
  if ((igen%parm->print_flag)==0)
    fprintf(fp," %6d        %2d  %12.6f   %12.4e   %12.4e   %5d  %5d\n",
	    igen,iniche,bestfit,avgfit,avg_best,nconv,nlost);
  fflush(fp);
}
/*------------------------------------------------------------------
/
/   evaluate_pop
/
/------------------------------------------------------------------*/
void evaluate_pop(int iniche,parm_struct* parm,population* pop)
{
  int i,j;

  for (i=0; i<pop->npop; i++) 
    pop->indiv[i].fitness = 
      evaluate_individual(iniche,parm,pop->indiv[i].chromosome);
}
/*------------------------------------------------------------------
/
/   evaluate_individual
/
/------------------------------------------------------------------*/
double evaluate_individual(int iniche,parm_struct* parm,int* chromosome)
{
  double fitness, param_list[MAXWORDS];

  filter(parm, chromosome, param_list);

  fitness = eval_user(iniche,parm->eval,parm->words,param_list);
  return fitness;
}
/*------------------------------------------------------------------
/
/   make_name
/
/------------------------------------------------------------------*/
void make_name(char* prefix,int n,char* suffix,char* name)
{
  char number[80];

  strcpy(name,prefix);
  sprintf(number,"%d",n);
  strcat(name,number);
  strcat(name,suffix);
}
/*------------------------------------------------------------------
/
/   make_name_2
/
/------------------------------------------------------------------*/
void make_name_2(char* prefix,int n,int j,char* suffix,char* name)
{
  char number[80];

  strcpy(name,prefix);
  sprintf(number,"%d",n);
  strcat(name,number);
  strcat(name,"_");
  sprintf(number,"%d",j);
  strcat(name,number);
  strcat(name,suffix);
}
/*------------------------------------------------------------------
/
/   filter
/
/------------------------------------------------------------------*/
void filter(parm_struct* parm,int* chromosome,double* array)
{
  int i, start, n, chr_bin[MAXLOCI];
  double min, max;

  for(i=0;i<parm->words;i++) {
    start = parm->start[i];
    n     = parm->length[i];
    degray(chromosome, chr_bin, start, n);
  }

  for(i=0;i<parm->words;i++) {
    start = parm->start[i];
    n     = parm->length[i];
    min   = parm->min[i];
    max   = parm->max[i];

    if(parm->type[i] == 0) 
      array[i] = (double)str_to_int(chr_bin,start,n,min,max);
    else if(parm->type[i] == 1) 
      array[i] = str_to_double(chr_bin,start,n,min,max);
  }
}
/*------------------------------------------------------------------
/
/   chr_to_double
/
/------------------------------------------------------------------*/
double str_to_double(int* chromosome,int start,int n,double min,double max)
{
  int i, factor, value;
  double doublevalue;

  factor = 1;
  value  = 0;
  for(i=start+n-1;i>=start;i--) {
    value += factor * chromosome[i];
    factor *= 2;
  }
  doublevalue = (double)value*(max-min)/(double)(factor-1) + min;
  return doublevalue;
}
/*------------------------------------------------------------------
/
/   chr_to_int
/
/------------------------------------------------------------------*/
int str_to_int(int* chromosome,int start,int n,double min,double max)
{
  int i, factor, value;
  double doublevalue;

  factor = 1;
  value  = 0;
  for(i=start+n-1;i>=start;i--) {
    value += factor * chromosome[i];
    factor *= 2;
  }
  
  doublevalue = (double)value*(max-min)/(double)(factor-1) + min;
  value = (int)doublevalue;

  return value;
}
/*------------------------------------------------------------------
/
/   initialize_best
/
/------------------------------------------------------------------*/
void initialize_best(parm_struct* parm, bestpop* best,int iniche)
{
  int i, j;
  char name[80];
  char cjunk1[80], cjunk2[80],line[80] ;
  double rjunk;
  int ijunk,k;
  FILE   *fp_rstrt;

  best->nbest = parm->nbest;

  if (!parm->restart || iniche<0) {
    for(i=0;i<parm->nbest;i++) {
      for(j=0;j<parm->chr_length;j++) 
	best->indiv[i].chromosome[j] = 1;
      best->indiv[i].fitness       = 100000000.0;
    }
  }
  
  else {

    /* Open GAbest<iniche>.dat file.
    ----------------------------- */
    make_name("GAbest",iniche,".dat",name);
    if((fp_rstrt = fopen(name,"r"))==NULL) {
      fprintf(stdout,"Cannot open %s in main()\n",name);
      return;
    }
    
    /* Read the best fitnesses and their chromosomes.
    ---------------------------------------------- */
    for(i=0;i<parm->nbest;i++) {
      fgets(line, 120, fp_rstrt);
      fgets(line, 120, fp_rstrt);
      sscanf(line,"%s%d%s%lf",&cjunk1,&ijunk,&cjunk2,
	     &(best->indiv[i].fitness)) ;
      
      for(j=0;j<parm->words;j++) {
	for(k=parm->start[j];k<parm->start[j]+parm->length[j];k++) {
	  fscanf(fp_rstrt,"%1d",&(best->indiv[i].chromosome[k]) ) ;
	  /*fprintf(stdout,"%1d",best->indiv[i].chromosome[k]) ;*/
	}
	fscanf(fp_rstrt,"%lf\n",&rjunk);
      }
    }
    fclose(fp_rstrt);
  }
}
/*------------------------------------------------------------------
/
/   initialize_pop
/
/------------------------------------------------------------------*/
void initialize_pop(parm_struct* parm,population* pop,int iniche)
{
  int i, j;
  double ran_num;
  int ijunk;
  char cjunk[15] ;
  char name[80];
  double rjunk, dummy;
  FILE   *fp_rstrt;
  
  pop->npop = parm->npop;

  if (!parm->restart || iniche<0) {
    for(i=0;i<parm->npop;i++) {
      for(j=0;j<parm->chr_length;j++) {
	ran_num = GAran(0);
	if(ran_num < 0.5) 
	  pop->indiv[i].chromosome[j] = 0;
	else
	  pop->indiv[i].chromosome[j] = 1;
      }
    }
  }
  else {
    /* Open GApop<iniche>.dat file.
    ----------------------------- */
    make_name("GApop",iniche,".dat",name);
    if((fp_rstrt = fopen(name,"r"))==NULL) {
      fprintf(stdout,"Cannot open %s in main()\n",name);
      return;
    }
    
    /* Read the last population of chromosomes from the previous run.
    --------------------------------------------------------------- */
    for(i=0;i<parm->npop;i++) {
      fscanf(fp_rstrt,"%d%s%lf",&ijunk,&cjunk,&(pop->indiv[i].fitness)) ;
      for(j=0;j<parm->chr_length;j++) {      
	fscanf(fp_rstrt,"%1d",&(pop->indiv[i].chromosome[j])) ;
      }
      fscanf(fp_rstrt,"\n") ;
    }
  }
}
/*------------------------------------------------------------------
/
/   mutate_pop
/
/------------------------------------------------------------------*/
void mutate_pop(parm_struct *parm, population *pop)
{
  int i, npop, j, flag;
  double ran_num;

  npop = pop->npop;
  for(i=1;i<npop;i++) {
    flag = 0;
    for(j=0;j<parm->chr_length;j++) {
      ran_num = GAran(0);
      if(ran_num < parm->mutate_rate && flag==0) {
	flag = 0;
	if(pop->indiv[i].chromosome[j] == 1)
	  pop->indiv[i].chromosome[j] = 0;
	else
	  pop->indiv[i].chromosome[j] = 1;
      }
    }
  }
}
/*------------------------------------------------------------------
/
/   best_pop
/
/------------------------------------------------------------------*/
void best_pop(population *pop, bestpop *oldbest,bestpop *newbest)
{
  int i, number;
  best_sort_struct bsort[2*MAXBEST];
  
  for(i=0;i<oldbest->nbest;i++) {
    bsort[i].number  = i;
    bsort[i].fitness = oldbest->indiv[i].fitness;
    bsort[i+oldbest->nbest].number  = i+oldbest->nbest;
    bsort[i+oldbest->nbest].fitness = pop->indiv[i].fitness;
  }
  
  qsort(bsort, 2*oldbest->nbest, sizeof(best_sort_struct), compare_bsort);
  newbest->nbest = oldbest->nbest;

  for(i=0;i<newbest->nbest;i++) {
    number = bsort[i].number;
    if(number<newbest->nbest)
      newbest->indiv[i] = oldbest->indiv[number];
    else
      newbest->indiv[i] = pop->indiv[number-oldbest->nbest];
  }
}
/*------------------------------------------------------------------
/
/   compare_bsort
/
/------------------------------------------------------------------*/
int compare_bsort(a,b)
best_sort_struct *a,*b;
{
  int flag;
  if(a->fitness <  b->fitness) flag = -1;
  if(a->fitness >  b->fitness) flag = 1;
  if(a->fitness == b->fitness) flag = 0;
  return flag;
}
/*------------------------------------------------------------------
/
/   select_pop_stepftn
/
/------------------------------------------------------------------*/
void select_pop_stepftn(parm_struct *parm,population *pop)
{
  int i, j, imax, scount, idum;
  int ichange[MAXPOP] ;
  double ran_num;
  
  for(i=0;i<pop->npop;i++)
    pop->indiv[i].select_flag = 0;

  /* Sort the fitnesses of the current population in order (low to high).
  -------------------------------------------------------------------- */
  qsort(pop->indiv, pop->npop, sizeof(individual), compare_pop);

  /* Determine which entries are duplicate structures.
  -------------------------------------------------- */
  for(i=0;i<pop->npop;i++) {
    ichange[i] = 0 ;
    
    if (i>0 && pop->indiv[i].fitness<1.e8 && pop->indiv[i].fitness==
	pop->indiv[i-1].fitness) {
      ichange[i] = 1;
    }
  }
  
  /* Preferentially select low energy non-duplicates
  ---------------------------------------------------- */
  for(i=0;i<pop->npop;i++) {
    pop->indiv[i].duplicate=0;
    if (ichange[i]==1) pop->indiv[i].duplicate = 1;
  }
  
  imax = min(pop->npop,(int)(pop->npop*parm->select_rate));
  
  scount=0;
  for(i=0;i<pop->npop;i++) {
    if(scount<=imax) {
      pop->indiv[i].select_flag = 1;
      scount++;
    }
  }
  /* hit the duplicates with a mutation 
  -------------------------------------*/
  for(i=0;i<pop->npop;i++) {
    if(pop->indiv[i].duplicate==1) {
      for(j=0;j<parm->chr_length;j++) {
	ran_num = GAran(0);
	if(ran_num < parm->mutate_rate) {
	  ran_num = GAran(0);
	  if(pop->indiv[i].chromosome[j]==0) pop->indiv[i].chromosome[j]=1;
	  else pop->indiv[i].chromosome[j]=0;
	}
      }
      pop->indiv[i].duplicate=0;
      pop->indiv[i].fitness = 
	evaluate_individual(0,parm,pop->indiv[i].chromosome);
    }
  }

}
/*------------------------------------------------------------------
/
/   select_pop_stepftn_old
/
/------------------------------------------------------------------*/
void select_pop_stepftn_old(parm_struct *parm,population *pop)
{
  int i, imax;
  int ichange[MAXPOP] ;
  
  for(i=0;i<pop->npop;i++)
    pop->indiv[i].select_flag = 0;

  /* Sort the fitnesses of the current population in order (low to high).
  -------------------------------------------------------------------- */
  qsort(pop->indiv, pop->npop, sizeof(individual), compare_pop);

  /* Determine which entries are duplicate structures.
  -------------------------------------------------- */
  for(i=0;i<pop->npop;i++) {
    ichange[i] = 0 ;
    
    if (i>0 && pop->indiv[i].fitness<1.e8 && pop->indiv[i].fitness==
	pop->indiv[i-1].fitness) {
      ichange[i] = 1;
    }
  }
  
  /* Set duplicate structures to the maximum fitness.
  ---------------------------------------------------- */
  for(i=0;i<pop->npop;i++) {
    if (ichange[i]==1) pop->indiv[i].fitness = 1.e8;
  }
  
  /* Re-sort the fitnesses to push the new max values to the bottom.
  -------------------------------------------------------------------- */
  qsort(pop->indiv, pop->npop, sizeof(individual), compare_pop);

  imax = min(pop->npop,(int)((double)pop->npop*parm->select_rate));
  
  for(i=0;i<imax;i++) pop->indiv[i].select_flag = 1;
}
/*------------------------------------------------------------------
/
/   select_pop_roulette
/
/------------------------------------------------------------------*/
void select_pop_roulette(parm_struct *parm,population *pop)
{
  int i, j,k,ptr,sample[MAXPOP];
  int ichange[MAXPOP] ;
  population  newpop;
  double Worst, sum, factor, expected;

  
  for(i=0;i<pop->npop;i++)
    pop->indiv[i].select_flag = 0;

  /* Sort the fitnesses of the current population in order (low to high).
  -------------------------------------------------------------------- */
  qsort(pop->indiv, pop->npop, sizeof(individual), compare_pop);

#ifdef REMOVE_DUP
  /* Determine which entries are duplicate structures.
  -------------------------------------------------- */
  for(i=0;i<pop->npop;i++) {
    ichange[i] = 0 ;
    
    if (i>0 && pop->indiv[i].fitness<1.e8 && pop->indiv[i].fitness==
	pop->indiv[i-1].fitness) {
      ichange[i] = 1;
    }
  }
  
  /* Set fitness of duplicate structures to the maximum fitness.
  ---------------------------------------------------- */
  Worst= -1.0e8;
  for(i=0;i<pop->npop;i++) {
    if (pop->indiv[i].fitness>Worst) Worst=pop->indiv[i].fitness;
  }
  for(i=0;i<pop->npop;i++) {
    if (ichange[i]==1) pop->indiv[i].fitness = Worst*1.05;;
  }
  
  /* Re-sort the fitnesses to push the new max values to the bottom.
  -------------------------------------------------------------------- */
  qsort(pop->indiv, pop->npop, sizeof(individual), compare_pop);
#endif

  /* The following code taken from select.c from GENESIS1.2 
   ----------------------------------------------------------*/
  /* denominator for selection probabilities */
  Worst = pop->indiv[pop->npop-1].fitness;
  for (sum = i = j = 0; i < pop->npop; i++)
    if (pop->indiv[i].fitness <= Worst) {
      sum += pop->indiv[i].fitness;
      j++;
    }
  
  factor = pop->npop/(Worst*j - sum);

  k = 0;	  /* index of next Selected structure */
  ptr = GAran(0);   /* spin the wheel one time */

  for (i=0; i<pop->npop;i++) 
    sample[i] = i;

  for (sum = i = 0; i < pop->npop; i++) {
    if (pop->indiv[i].fitness <= Worst)
      expected = (Worst - pop->indiv[i].fitness) * factor;
    else expected = 0.0;
    
    for (sum += expected; sum > ptr; ptr++)
      sample[k++] = i;
    if(k==pop->npop) break;
  }
  if (k != pop->npop) {
/*    fprintf(stdout,"select_pop_roulette: internal scaling error: %d\n",k);*/
  }
  
  /* form the new population */
  for (i = 0; i < pop->npop; i++) {
    k = sample[i];
    newpop.indiv[i] = pop->indiv[k];
  }
  
  for (i = 0; i < pop->npop; i++) {
    pop->indiv[i]= newpop.indiv[i];
    pop->indiv[i].select_flag = 1;
  }

}
/*------------------------------------------------------------------
/
/   compare_pop
/
/------------------------------------------------------------------*/
int compare_pop(a,b)
individual *a,*b;
{
  int flag;
  if(a->fitness <  b->fitness) flag = -1;
  if(a->fitness >  b->fitness) flag = 1;
  if(a->fitness == b->fitness) flag = 0;
  return flag;
}
/*------------------------------------------------------------------
/
/   crossover_pop
/
/------------------------------------------------------------------*/
void crossover_pop(parm_struct* parm,population* oldpop,population* newpop, int igen)
{
  int i,j,jniche,npop, nsel, n1, n2, ngap,count,accept,cham,nelite;
  individual in1, in2;
  
  count=0;

  npop = oldpop->npop;
  newpop->npop = npop;

  /* Use elitest strategy, i.e. save a copy of the best individuals
     from the current genration for the next generation .
  -------------------------------------------------------------- */
  nelite=1;
  for(i=0;i<nelite;i++) {
    newpop->indiv[i] = oldpop->indiv[0];
    count++;
  }
  
  ngap=(int)(npop*(1.0-parm->crossover_rate));
  for(i=0;i<ngap;i++) {
    j=parm->npop*GAran(0);
    if(oldpop->indiv[j].select_flag) {
      newpop->indiv[i+count] = oldpop->indiv[j];
      count++;
    }
  }
  
  /* If best values were shared between the niches save those values.
  ---------------------------------------------------------------- */
  if (parm->niches > 1) 
    for (i=0; i<parm->niflag; i++) 
      if (igen==parm->ni_gen[i])  
	for(jniche=1;jniche<parm->niches;jniche++) {
	  newpop->indiv[count] = oldpop->indiv[jniche] ;
	  count++ ;
	}

  while(count<npop-1) {
    n1 = min(npop-1,(int)(GAran(0)*(double)npop));
    n2 = min(npop-1,(int)(GAran(0)*(double)npop));
    if((n1!=n2) &&
       oldpop->indiv[n1].select_flag == 1 &&
       oldpop->indiv[n2].select_flag == 1 
       ) {

      crossover(parm, &oldpop->indiv[n1], &oldpop->indiv[n2],&in1,&in2);
      newpop->indiv[count]=in1;
      count++;
      newpop->indiv[count]=in2;
      count++;
    }
  }
}
/*------------------------------------------------------------------
/
/   crossover_pop_old
/
/------------------------------------------------------------------*/
void crossover_pop_old(parm_struct* parm,population* oldpop,population* newpop, int igen)
{
  int i,jniche,npop, n1, n2, count,accept,cham,nelite;
  individual in1, in2;
  double ham_cut;
 
  ham_cut = 0.4;
  count=0;

  npop = oldpop->npop;
  newpop->npop = npop;

  /* Use elitest strategy, i.e. save a copy of the best individuals
     from the current genration for the next generation .
  -------------------------------------------------------------- */
  nelite=npop/10;
  for(i=0;i<nelite;i++) {
    newpop->indiv[i] = oldpop->indiv[0];
    count++;
  }
  
  for(i=0;i<nelite;i++) {
    newpop->indiv[i+count] = oldpop->indiv[i+1];
    count++;
  }
  
  /* If best values were shared between the niches save those values.
  ---------------------------------------------------------------- */
  if (parm->niches > 1) 
    for (i=0; i<parm->niflag; i++) 
      if (igen==parm->ni_gen[i])  
	for(jniche=1;jniche<parm->niches;jniche++) {
	  newpop->indiv[count] = oldpop->indiv[jniche] ;
	  count++ ;
	}

  cham=0;
  while(count<npop-1) {
    cham++;
    if(cham >= npop*10) {
      ham_cut*=0.8;
      fprintf(stdout,"Rescale ham_cut after %6d steps to: %8.4f\n",cham,ham_cut);
      cham=0;
    }
    if(ham_cut<1.0e-6) ham_cut= -1.0;
    n1 = min(npop-1,(int)(GAran(0)*(double)npop));
    n2 = min(npop-1,(int)(GAran(0)*(double)npop));
    
    if((n1!=n2) &&
       oldpop->indiv[n1].select_flag == 1 &&
       oldpop->indiv[n2].select_flag == 1 ) {
      crossover(parm, &oldpop->indiv[n1], &oldpop->indiv[n2],&in1,&in2);
      accept=1;
      for(i=0;i<count;i++) 
	if(hamming_distance(parm,&in1,&newpop->indiv[i])<ham_cut) accept=0;
      if(accept) {
	newpop->indiv[count]=in1;
	count++;
      }
      accept=1;
      for(i=0;i<count;i++) 
	if(hamming_distance(parm,&in2,&newpop->indiv[i])<ham_cut) accept=0;
      if(accept) {
	newpop->indiv[count]=in2;
	count++;
      }
    }
  }
  fprintf(stdout,"Number of crossovers: %6d vs. population: %6d\n",cham,npop);
}
/*------------------------------------------------------------------
/
/   crossover
/
/------------------------------------------------------------------*/
void crossover(parm_struct* parm,individual* parent1,individual* parent2,
	       individual* child1,individual* child2)
{
  int ncut, i;
  
  for(i=0;i<parm->chr_length;i++)
    child1->chromosome[i] = child2->chromosome[i] = 0;

  ncut = (int)(GAran(0) * (double)parm->chr_length);

  if(ncut < parm->chr_length && ncut > 0) {
    for(i=0;i<ncut;i++) {
      child1->chromosome[i] = parent1->chromosome[i];
      child2->chromosome[i] = parent2->chromosome[i];
    }

    for(i=ncut;i<parm->chr_length;i++) {
      child1->chromosome[i] = parent2->chromosome[i];
      child2->chromosome[i] = parent1->chromosome[i];
    }
  }
 
  else {
    for(i=0;i<parm->chr_length;i++) {
      child1->chromosome[i] = parent1->chromosome[i];
      child2->chromosome[i] = parent2->chromosome[i];
    }
  }
}
/*------------------------------------------------------------------
/
/   hamming_distance
/
/------------------------------------------------------------------*/
double hamming_distance(parm_struct* parm, individual* i1,individual* i2)
{
  double r;
  int i, ch;

  ch=0;
  for(i=0;i<parm->chr_length;i++) {
    if(i1->chromosome[i] != i2->chromosome[i])
      ch++;
  }
  r=(double)ch/(double)parm->chr_length;

  return r;
}
/*------------------------------------------------------------------
/
/   GAran
/
/------------------------------------------------------------------*/
double GAran(int iseed)
{
  double ran_num;
  int idum;

  if(iseed == 0 ) {
    idum    = 1;
    ran_num = ran3(&idum);
  }
  else {
    idum   = -1;
    ran_num = ran3(&idum);
  }
  return ran_num;
}
/*------------------------------------------------------------------
/
/   ran3
/
/------------------------------------------------------------------*/
#define MBIG 1000000000
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(int* idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  static int MSEED;

  if(*idum > 0)
    MSEED = (*idum);

  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}

#undef MBIG
#undef MZ
#undef FAC

/*------------------------------------------------------------------
/
/   degray
/
/------------------------------------------------------------------*/
void degray(int* chr_gray,int* chr_bin,int start,int n)
{
  int i;
  int last;

  last = 0;
  for (i = start; i < start+n; i++) {
    if (chr_gray[i])  chr_bin[i] = !last;
    else chr_bin[i] = last;
    last = chr_bin[i];
  }
}

int max(int i,int j)
{
  if(i<j) 
    return j;
  else
    return i;
}

int min(int i,int j)
{
  if(i<j) 
    return i;
  else
    return j;
}

double rmax(double i,double j)
{
  if(i<j) 
    return j;
  else
    return i;
}

double rmin(double i,double j)
{
  if(i<j) 
    return i;
  else
    return j;
}
