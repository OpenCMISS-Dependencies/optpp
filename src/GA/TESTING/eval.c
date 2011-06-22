/*-----------------------------------------------------------------
/
/   eval.c - sample eval file
/------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include <sys/types.h>
#include <sys/stat.h>
#include "eval.h"
#include "vector.h"
#include "polymer.h"

/*------------------------------------------------------------------
/
/   eval_user - To add your own eval function, called, for instance
/               my_eval, 
/
/   1. add the function my_eval(int,double*) to the bottom of this file
/
/   2. put the prototype in eval.h 
/
/   3. in eval_user add the lines as in the sample below, 
/
/   4. recompile (make pga)
/
/   5. change the first line in GAin.dat to my_eval
/
/   6. set the number and range of parameters at the bottom of GAin.dat
/
/   7. run pga
/
/------------------------------------------------------------------*/
double eval_user(int iniche,char *f,int n, double* params)
{
  double fitness;

  if(!strcmp(f,"eval0"))
    fitness=eval0(n,params);
  else if(!strcmp(f,"eval1"))
    fitness=eval1(n,params);

/* sample new lines ..._____________
  else if(!strcmp(f,"my_eval"))
    fitness=my_eval(n,params);
 ..._________________________________*/

  else {
    fprintf(stdout,"You have specified a bogus function name: %s\n",f);
    exit(0);
  }
  return fitness;
}
/*------------------------------------------------------------------
/
/   eval1 - The 10 compound QSAR function
/
/------------------------------------------------------------------*/
double eval1(int nx, double* x)
{
  double v[10]={250.0,250.0,245.0,245.0,265.0,265.0,260.0,260.0,100.0,120.0};
  int    nhb[10]={5,4,5,4,5,4,5,4,3,3};
  static double logki[10];
  double xtrue[4]={0.35,0.78,4.07,1.45};
  int i, n=10;
  static int ifirst=1;
  double logki_calc,fitness=0.0;

  /* calculate the "true" logki values
  --------------------------------------*/
  if(ifirst) {
    
    fprintf(stdout,"--------------------------------------------------\n");
    fprintf(stdout,"    True parameters for test compounds\n");
    fprintf(stdout,"--------------------------------------------------\n");
    for(i=0;i<10;i++) {
      logki[i]=10.0+xtrue[0]*pow(v[i],xtrue[1]) - xtrue[2]*pow(nhb[i],xtrue[3]);
      fprintf(stdout,"compound:  %3d Mol. vol. %8.4f  nhb: %3d logKi: %8.4f\n",
	      i,v[i],nhb[i],logki[i]);
    }
    ifirst=0;
  }

  for(i=0;i<n;i++) {
    logki_calc=10.0+x[0]*pow(v[i],x[1]) - x[2]*pow(nhb[i],x[3]);
    fitness+=(logki[i]-logki_calc)*(logki[i]-logki_calc);
  }
  fitness=sqrt(fitness/n);
  return fitness;
}
/*------------------------------------------------------------------
/
/   eval0 - this is the required format of the user defined
/           eval function. n is the number of double precision
/           paramters contained in params
/
/------------------------------------------------------------------*/
double eval0(int n, double* params)
{
  int i;
  double fitness;

  fitness=0.0;
  for(i=0;i<n;i++)
    fitness+=params[i]*params[i];

  return fitness;
}



