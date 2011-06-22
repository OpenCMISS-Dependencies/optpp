/*------------------------------------------------------------------------
/
/   polymer.c
/
/   This file contains functions that retunr the energy of a 2-D 
/   polymer after either a folding seqeunce (eval_fold) or a specification
/   of a set of internal angles (eval_cart). In both cases, gradient
/   minimization is done using CG routines from the book Numerical
/   Recipes in C by Press, et.al. Input for the function is given 
/   in their individual preambles. Other files needed are dnumrec.c,
/   vector.c, polymer.h and sequence.dat.
/
/   The standard test cases are 7, 19, 37 and 61 atoms for which the 
/   global energy minima are:
/
/    7:   -6.5340
/   19:  -27.3433  
/   37:  -62.4579
/   61: -111.8109
/
/
/------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
c  Copyright (C) 1995
c
c  Author:
c
c  Richard Judson       
c  rsjuds@ca.sandia.gov  
c  (510)294-1438           
c  FAX (510)294-2234
c
c  Center for Computational Engineering, M.S. 9214
c  Sandia National Laboratories
c  P.O. Box 969
c  Livermore, CA 94551-0969
c
c------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "pga.h"
#include "vector.h"
#include "polymer.h"

int         n_atoms;
atom_struct polymer[MAXATOMS];
FILE        *fpe;
int         spec1= -1, spec2= -1;
double      beste[MAXBEST];
double      bestfitness[MAXBEST];
extern void dfpmin();
double dist_matrix[MAXATOMS*MAXATOMS];
double dist_matrix_0[MAXATOMS*MAXATOMS];
int constrain;

/*-----------------------------------------------------------
/
/   eval_cart - Calculate the energy given a set of angles
/
/   input: iniche (int) - tells which of several 
/             populations the current individual is in.
/          angles (double *) - array holding the n_atoms-2
/             angles. Values go from [0->2pi].
/
/   input file: sequence.dat
/             first line is the number of atoms.
/             subsequent lines give forcefield parameters.
/             For the "standard problem", each line of 
/             sequence.dat should be:
/
/             0.50  100.0  0.50  1.0
/
/
/------------------------------------------------------------*/
double eval_cart(int n, double *angles)
{
  int        i, j, k, iter, itmax, sign;
  double     energy, theta, x[2*MAXATOMS], ftol, fret;
  char name[80], number[10];
  static int ifirst=1;
  static double PI;

  /* read the sequence on first pass */
  if(ifirst) {
    PI = acos(-1.0);
    read_seq();
    ifirst=0;
    for(i=0;i<MAXBEST;i++) {
      beste[i]=1000000.0;
      bestfitness[i]=1000000.0;
    }

  }


  /* make the polymer */
  make_polymer( angles );

  j = 0;
  for(k=0;k<n_atoms;k++) {
    x[++j] = polymer[k].R.x;
    x[++j] = polymer[k].R.y;
  }

  /* calculate the basic energy */
  energy=e_total(-1,-1);
  ftol = 0.001; 
  itmax = 100;  
  spec1=spec2= -1;
/*  dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func,dfunc); */

  j = 0;
  for(k=0;k<n_atoms;k++) {
    polymer[k].R.x = x[++j];
    polymer[k].R.y = x[++j];
  }

/*  if(PRINT_ENERGY_2) printf("FINAL ENERGY: %8.2f\n",energy); */

/*  if(energy < beste[iniche]) {
 *   beste[iniche] = energy;
 *   strcpy(name,"best");
 *   sprintf(number,"%d",iniche);
 *   strcat(name,number);
 *   strcat(name,".dat");
 *   write_mm(beste[iniche],name);
 * }
 */

  return energy;
}
/*-----------------------------------------------------------
/
/   eval_fold - Calculate the energy after a folding
/               sequence.
/
/   input: iniche (int) - tells which of several 
/             populations the current individual is in.
/          foldseq (int *) - array holding the folding
/             folding sequence. The first integer gives
/             the number of folds to perform. Subsequent
/             pairs give (1) the first of the two atoms
/             to be folded together and (2) the distance 
/             down the chain where teh second atom sits.
/             For instance [2,6] would fold atoms 2 and 8
/             together. For 19 atoms, the range of the first
/             value in each pair is [0->18] (atom numbers 
/             start at 0) and the range of the second value is
/             [2->17].
/
/   input file: sequence.dat
/             first line is the number of atoms.
/             subsequent lines give forcefield parameters.
/             For the "standard problem", each line of 
/             sequence.dat should be:
/
/             0.50  100.0  0.50  1.0
/
/   output:
/        energy: energy after the final gradient minimization
/
/------------------------------------------------------------*/
double eval_fold(int iniche, double* foldseq)
{
  int        i, j, k, iter, itmax, sign;
  double     PI, energy, theta, x[2*MAXATOMS], 
             ftol, fret, angles[MAXATOMS];
  char name[80], number[10];
  static int ifirst=1;

  /* read the sequence and open history.dat on first pass */
  if(ifirst) {
    for(i=0;i<MAXBEST;i++) {
      beste[i]=1000000.0;
      bestfitness[i]=1000000.0;
    }

    read_seq();
    if(HISTORY) {
      strcpy(name,"history");
      sprintf(number,"%d",iniche);
      strcat(name,number);
      strcat(name,".dat");
      fpe=fopen(name,"w");
      fprintf(fpe,"%4d\n",n_atoms);
    }
    ifirst=0;
  }

  PI = acos(-1.0);

  /* make the polymer */
  sign = 1;
  for(i=0;i<n_atoms-2;i++)  {
    sign *= -1;
    angles[i] = sign * 0.1;
  }

  /* print the folding sequence into history.dat */
  if(!ifirst && HISTORY ) {
    fprintf(fpe,"%3d",(int)foldseq[0]);
    
    for(i=1;i<=foldseq[0];i++) 
      fprintf(fpe,"%3d%3d",(int)foldseq[2*i-1],(int)foldseq[2*i]);
    fflush(fpe);
  }

  make_polymer( angles );

  if(PRINT_ENERGY_1) {
    printf("-------------------------------------------------\n");
    printf("Folds: %3d",(int)foldseq[0]);
    for(i=1;i<=(int)foldseq[0];i++)  {
      spec1 = (int)foldseq[2*i-1];
      spec2 = (spec1 + (int)foldseq[2*i])%n_atoms;
      printf(" [%2d %2d]",spec1,spec2);
    }
    printf("\n");
  }

  ftol = 0.001;
  itmax = 100;
  j = 0;
  for(k=0;k<n_atoms;k++) {
    x[++j] = polymer[k].R.x;
    x[++j] = polymer[k].R.y;
  }

  if(RELAX) {
    for(i=1;i<=(int)foldseq[0];i++) {
      spec1 = (int)foldseq[2*i-1];
      spec2 = (spec1 + (int)foldseq[2*i])%n_atoms;
      dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func,dfunc);

      if(PRINT_ENERGY_2) {
	spec1=spec2= -1;
	energy=e_total(spec1, spec2); 
	printf("Iterations: %5d Energy: %8.2f\n",iter,energy);
      }
    }

    j = 0;
    for(k=0;k<n_atoms;k++) {
      theta  = 2.0*PI*k / n_atoms;
      x[++j]+= 0.25*cos(theta);
      x[++j]+= 0.25*sin(theta);
    }
  }

  /* do a final minimization */
  ftol = 0.001; 
  itmax = 100;  
  spec1=spec2= -1;
  dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func,dfunc);

  if(PRINT_ENERGY_2) printf("FINAL ENERGY: %8.2f\n",energy);

  /* print the energy into history.dat */
  if(!ifirst && HISTORY) 
    fprintf(fpe," %12.4f\n",energy);

  if(energy < beste[iniche]) {
    beste[iniche] = energy;
    strcpy(name,"best");
    sprintf(number,"%d",iniche);
    strcat(name,number);
    strcat(name,".dat");
    write_mm(beste[iniche],name);
  }

  return energy;
}
/*-----------------------------------------------------------
/
/   write_pdb - write a pdb file
/
/------------------------------------------------------------*/
void write_pdb(double energy, char* fname)
{
  int i,j, ip;
  FILE *fp;
  char name[5],resname[5],resid[5];

  if((fp = fopen(fname,"w"))==NULL) {
    printf("Cannot open %s in write_pdb()\n",fname);
    return;
  }

  fprintf(fp,"REMARK  %4d      %s     %8.2f\n",n_atoms,fname,energy);

  for (i=0;i<n_atoms;i++) {
    strcpy(name,"   C");
    strcpy(resname,"GLY");
    strcpy(resid,"   1");
    fprintf(fp,"ATOM  %5d %4.4s %3.3s  %4.4s    ",
	    i+1,name,resname,resid);
    fprintf(fp," %7.3f %7.3f %7.3f %5.2f %5.2f\n",
	    polymer[i].R.x, 
	    polymer[i].R.y, 
	    polymer[i].R.z,
	    0.0,0.0);
  }
  i=0;
  fprintf(fp,"END\n");
  fclose(fp);
}
/*-----------------------------------------------------------
/
/   write_dat - write an xmgr file
/
/------------------------------------------------------------*/
void write_dat(double energy, char* fname)
{
  int i,j, ip;
  FILE *fp;
  char name[5],resname[5],resid[5];

  if((fp = fopen(fname,"w"))==NULL) {
    printf("Cannot open %s in write_pdb()\n",fname);
    return;
  }

  for (i=0;i<n_atoms;i++) 
    fprintf(fp," %7.3f %7.3f\n",polymer[i].R.x, polymer[i].R.y);
  fclose(fp);
}
/*-----------------------------------------------------------
/
/   write_mm - write a macromodel file
/
/------------------------------------------------------------*/
void write_mm(double energy, char* fname)
{
  int i,j, ip;
  FILE *fp;

  if((fp = fopen(fname,"w"))==NULL) {
    printf("Cannot open %s in write_mm()\n",fname);
    return;
  }

  fprintf(fp,"%4d      %s     %8.2f\n",n_atoms,fname,energy);

  i=0;
  fprintf(fp,"%4d  %4d%2d  %4d%2d  %4d%2d  %4d%2d  %4d%2d  %4d%2d%12.6f%12.6f%12.6f\n",
	  42,i+2,1,0,0,0,0,0,0,0,0,0,0,
	  polymer[i].R.x, 
	  polymer[i].R.y, 
	  polymer[i].R.z); 

  for(i=1;i<n_atoms-1;i++) 
    fprintf(fp,"%4d  %4d%2d  %4d%2d  %4d%2d  %4d%2d  %4d%2d  %4d%2d%12.6f%12.6f%12.6f\n",
	    42,i+2,1,i,1,0,0,0,0,0,0,0,0,
	    polymer[i].R.x, 
	    polymer[i].R.y, 
	    polymer[i].R.z); 

  i=n_atoms-1;
  fprintf(fp,"%4d  %4d%2d  %4d%2d  %4d%2d  %4d%2d  %4d%2d  %4d%2d%12.6f%12.6f%12.6f\n",
	  42,i,1,0,0,0,0,0,0,0,0,0,0,
	  polymer[i].R.x, 
	  polymer[i].R.y, 
	  polymer[i].R.z); 


  fclose(fp);
}
/*-----------------------------------------------------------
/
/   func - energy function for the numerical 
/          recipes minimization routines
/
/------------------------------------------------------------*/
double func(double *x)
{
  int i, j;
  double energy;

  j = 0;
  for(i=0;i<n_atoms;i++) {
    polymer[i].R.x = x[++j];
    polymer[i].R.y = x[++j];
  }
  energy=e_total(spec1, spec2);
  if(DEBUGEVAL) printf(" [func] energy: %8.2f\n",energy);
  return energy;
}
/*-----------------------------------------------------------
/
/   dfunc - derivative function for the numerical 
/           recipes minimization routines
/
/------------------------------------------------------------*/
void dfunc(double *x,double * dx)
{
  int i, j;
  double energy;

  j = 0;
  for(i=0;i<n_atoms;i++) {
    polymer[i].R.x = x[++j];
    polymer[i].R.y = x[++j];
  }

  energy=e_total(spec1, spec2);
  if(DEBUGEVAL) printf(" [dfunc] energy: %8.2f\n",energy);

  j = 0;
  for(i=0;i<n_atoms;i++) {
    dx[++j] = polymer[i].f.x;
    dx[++j] = polymer[i].f.y;
  }
}
/*-----------------------------------------------------------
/
/   func2 - energy function for the numerical 
/          recipes minimization routines
/
/------------------------------------------------------------*/
double func2(double *x)
{
  int i, j;
  double energy;

  j = 0;
  for(i=0;i<n_atoms;i++) {
    polymer[i].R.x = x[++j];
    polymer[i].R.y = x[++j];
  }
  energy=e_constrain(constrain);
  if(DEBUGEVAL) printf(" [func] energy: %8.2f\n",energy);
  return energy;
}
/*-----------------------------------------------------------
/
/   dfunc2 - derivative function for the numerical 
/           recipes minimization routines
/
/------------------------------------------------------------*/
void dfunc2(double *x,double * dx)
{
  int i, j;
  double energy;

  j = 0;
  for(i=0;i<n_atoms;i++) {
    polymer[i].R.x = x[++j];
    polymer[i].R.y = x[++j];
  }

  energy=e_constrain(constrain);
  if(DEBUGEVAL) printf(" [dfunc] energy: %8.2f\n",energy);

  j = 0;
  for(i=0;i<n_atoms;i++) {
    dx[++j] = polymer[i].f.x;
    dx[++j] = polymer[i].f.y;
    if(DEBUGEVAL) printf(" [dfunc] j %3d fx,fy: %8.2f  %8.2f\n",j,dx[j-1],dx[j-2]);
  }
}
/*-----------------------------------------------------------
/
/   func3 - energy function for the numerical 
/          recipes minimization routines
/
/------------------------------------------------------------*/
double func3(double *x)
{
  int i, j;
  double energy;

  j = 0;
  for(i=0;i<n_atoms;i++) {
    polymer[i].R.x = x[++j];
    polymer[i].R.y = x[++j];
  }
  energy=e_total_q(constrain);
  if(DEBUGEVAL) printf(" [func] energy: %8.2f\n",energy);
  return energy;
}
/*-----------------------------------------------------------
/
/   dfunc3 - derivative function for the numerical 
/           recipes minimization routines
/
/------------------------------------------------------------*/
void dfunc3(double *x,double * dx)
{
  int i, j;
  double energy;

  j = 0;
  for(i=0;i<n_atoms;i++) {
    polymer[i].R.x = x[++j];
    polymer[i].R.y = x[++j];
  }

  energy=e_total_q(constrain);
  if(DEBUGEVAL) printf(" [dfunc] energy: %8.2f\n",energy);

  j = 0;
  for(i=0;i<n_atoms;i++) {
    dx[++j] = polymer[i].f.x;
    dx[++j] = polymer[i].f.y;
    if(DEBUGEVAL) printf(" [dfunc] j %3d fx,fy: %8.2f  %8.2f\n",j,dx[j-1],dx[j-2]);
  }
}
/*-----------------------------------------------------------
/
/   make_polymer 
/
/------------------------------------------------------------*/
void make_polymer( double *Phi )
{
  int i;
  double R0, Theta[MAXATOMS];

  /* set the first two atoms */
  polymer[0].R.x = polymer[0].R.y = polymer[0].R.z = 0.0;
  polymer[1].R.y = polymer[1].R.z = 0.0;
  R0 = polymer[0].Rbond + polymer[1].Rbond;
  polymer[1].R.x = R0;

  Theta[0] = Phi[0];
  for(i=1;i<n_atoms-2;i++)
    Theta[i] = Phi[i] + Theta[i-1];

  /* loop over the rest of the atoms */
  for(i=2;i<n_atoms;i++) {

    R0 = polymer[i-1].Rbond + polymer[i].Rbond;
    polymer[i].R.x = polymer[i-1].R.x + R0*cos(Theta[i-2]);
    polymer[i].R.y = polymer[i-1].R.y + R0*sin(Theta[i-2]);
    polymer[i].R.z = 0.0;

  }
}
/*-----------------------------------------------------------
/
/   distance - give the distance between two atoms
/
/------------------------------------------------------------*/
double distance(atom_struct* atom_1,atom_struct* atom_2)
{
  double r;
  r = sqrt( (atom_1->R.x - atom_2->R.x) * (atom_1->R.x - atom_2->R.x) +
	    (atom_1->R.y - atom_2->R.y) * (atom_1->R.y - atom_2->R.y) +
	    (atom_1->R.z - atom_2->R.z) * (atom_1->R.z - atom_2->R.z) );
  return r;
}
/*-----------------------------------------------------------
/
/   cos_angle - returns the cosine of the angle 
/               defined by three atoms
/
/------------------------------------------------------------*/
double cos_angle(atom_struct* atom_1, atom_struct* atom_2,atom_struct*  atom_3)
{
  double r12, r23, r13, cosine;

  r12=distance(atom_1, atom_2);
  r23=distance(atom_2, atom_3);
  r13=distance(atom_1, atom_3);
  if( r12==0.0 || r23==0.0 || r12==0.0) {
    printf(" *** error: two atoms overlap ...\n");
    print_atom(atom_1);
    print_atom(atom_2);
    print_atom(atom_3);
    exit(0);
  }
  cosine = (r23*r23+r12*r12-r13*r13)/(2.0*r12*r23);
  return cosine;
}
/*-----------------------------------------------------------
/
/   e_vdw - gives the energy and derivatives for a 
/           Van der Waals interaction
/
/-----------------------------------------------------------*/
double e_vdw(double Rvdw,double Dvdw,vector* a1,vector* a2,vector* f1,vector* f2)
{
  double r,q,q6,e,dedr,f,rcut;
  vector rv;

  rv   = vec_diff(a2,a1);
  r    = vec_length(&rv);
  rcut = 0.2;
  if(r < rcut) r=rcut;
  if(r<5.0) {
    q     = Rvdw/r;
    q6    = q*q*q*q*q*q;
    e     = Dvdw * q6*(q6 - 2.0);
    dedr  = -12.0*Dvdw * q6*(q6 - 1.0)/r;

    f     = dedr/r;
    f2->x = f*(a2->x - a1->x);
    f2->y = f*(a2->y - a1->y);
    f2->z = 0.0;
  }
  else {
    e=0.0;
    f2->x = 0.0;
    f2->y = 0.0;
    f2->z = 0.0;
  }
  *f1 = vec_mult(f2,-1.0);

  return e;
}
/*-----------------------------------------------------------
/
/   e_q - gives the energy and derivatives for a 
/           Van der Waals interaction
/
/-----------------------------------------------------------*/
double e_q(double q1, double q2,vector* a1,vector* a2,vector* f1,vector* f2)
{
  double r,q,q6,e,dedr,f,rcut;
  vector rv;
  double efac=332.0636;

  rv   = vec_diff(a2,a1);
  r    = vec_length(&rv);
  rcut = 0.2;
  if(r < rcut) r=rcut;

  e     = efac*q1*q2/r;
  dedr  = -efac*q1*q2/r/r;

  f     = dedr/r;
  f2->x = f*(a2->x - a1->x);
  f2->y = f*(a2->y - a1->y);
  f2->z = 0.0;
  *f1 = vec_mult(f2,-1.0);

  return e;
}
/*-----------------------------------------------------------
/
/   e_bond - gives the energy and derivatives for a bond
/
/-----------------------------------------------------------*/
double e_bond(vector* a1,vector* a2,double r0,double k,
	      vector* f1,vector* f2)
{
  double r,e,d,f;
  vector rv;

  rv = vec_diff(a2,a1);
  r  = vec_length(&rv);

  d  = (r-r0);
  e  = 0.5*k*d*d;
  f  = k*d/r;
  f2->x = f*(a2->x - a1->x);
  f2->y = f*(a2->y - a1->y);
  f2->z = f*(a2->z - a1->z);
  *f1 = vec_mult(f2,-1.0);

  return e;
}
/*-----------------------------------------------------------
/
/   e_total - calculate the total energy and the forces
/             if desired, for the polymer
/
/------------------------------------------------------------*/
double e_total(int spec1, int spec2)
{
  double e, e1, e2, e10, e20, Rbond, Kbond, Rvdw, Dvdw;
  int    i, j;
  vector a1,a2,f1,f2;

  e=e1=e2=0.0;

  /* zero out the forces */
  for(i=0;i<n_atoms;i++) {
    polymer[i].f.x=0.0;
    polymer[i].f.y=0.0;
    polymer[i].f.z=0.0;
  }

  /* bonding contributions */
  for(i=0;i<n_atoms-1;i++) {
    a1    = polymer[i].R;
    a2    = polymer[i+1].R;
    Rbond = polymer[i].Rbond + polymer[i+1].Rbond;
    Kbond = (polymer[i].Kbond + polymer[i+1].Kbond)/2.0;
    e10=e_bond(&a1,&a2,Rbond,Kbond,&f1,&f2);
    e1+=e10;

    polymer[i].f.x += f1.x;    
    polymer[i].f.y += f1.y;    
    polymer[i].f.z += f1.z;    
    polymer[i+1].f.x += f2.x;
    polymer[i+1].f.y += f2.y;
    polymer[i+1].f.z += f2.z;
  }

  /* special folding bond */

  if(spec1 >= 0 && spec2 >= 0) {
    a1    = polymer[spec1].R;
    a2    = polymer[spec2].R;
    Rbond = polymer[spec1].Rbond + polymer[spec2].Rbond;
    Kbond = 0.25*(polymer[spec1].Kbond + polymer[spec2].Kbond)/2.0;
    e10=e_bond(&a1,&a2,Rbond,Kbond,&f1,&f2);
    e1+=e10;

    polymer[spec1].f.x += f1.x;    
    polymer[spec1].f.y += f1.y;    
    polymer[spec1].f.z += f1.z;    
    polymer[spec2].f.x += f2.x;
    polymer[spec2].f.y += f2.y;
    polymer[spec2].f.z += f2.z;
  }

  /* calculate the Van der Waals repulsion term */
  for(i=0;i<n_atoms-2;i++) {
    for(j=i+2;j<n_atoms;j++) {
      a1    = polymer[i].R;
      a2    = polymer[j].R;
      Rvdw  = polymer[i].Rvdw + polymer[j].Rvdw;
      Dvdw  = (polymer[i].Dvdw + polymer[j].Dvdw)/2.0;
      e20=e_vdw(Rvdw,Dvdw,&a1,&a2,&f1,&f2);
      e2+=e20;
      
      polymer[i].f.x += f1.x;    
      polymer[i].f.y += f1.y;    
      polymer[i].f.z += f1.z;    
      polymer[j].f.x += f2.x;
      polymer[j].f.y += f2.y;
      polymer[j].f.z += f2.z;
    }
  }
  e = e1+e2;
  return e;
}
/*--------------------------------------------------------------
/
/   read_seq - reads in the sequence file
/
--------------------------------------------------------------*/
void read_seq()
{
  int    i;
  FILE   *infile;
  char line[120];

  if((infile = fopen("sequence.dat","r"))==NULL) {
    printf("*** Fatal error: file sequence.dat cannot be opened\n");
    exit(0);
  }
  
  fgets(line,80,infile);
  sscanf(line,"%d",&n_atoms);
  if(n_atoms > MAXATOMS) {
    printf("n_atoms: %d > MAXATOMS: %d\n",n_atoms,MAXATOMS);
    exit(0);
  }

  for(i=0;i<n_atoms;i++) 
    fscanf(infile,"%lf %lf %lf %lf",
	   &polymer[i].Rbond,
	   &polymer[i].Kbond,
	   &polymer[i].Rvdw,
	   &polymer[i].Dvdw);
  fclose(infile);

}
/*-----------------------------------------------------------
/
/   print_atom - print the information on an atom
/
/------------------------------------------------------------*/
void print_atom(atom_struct *atom)
{
  int atom_id();

  printf("------------------------------ atom info ---------------------------\n");
  printf("| x,y,z:                   %8.3f %8.3f %8.3f\n",atom->R.x,atom->R.y,atom->R.z);
  printf("--------------------------------------------------------------------\n");
}
/*-----------------------------------------------------------
/
/   eval_state_table - Calculate the energy after a folding
/               sequence.
/
/   input: iniche (int) - tells which of several 
/             populations the current individual is in.
/   foldseq (int *) - 
/
/   input file: sequence.dat
/             first line is the number of atoms.
/             subsequent lines give forcefield parameters.
/             For the "standard problem", each line of 
/             sequence.dat should be:
/
/             0.50  100.0  0.50  1.0
/
/   output:
/        energy: energy after the final gradient minimization
/
/------------------------------------------------------------*/
double eval_state_table(int iniche, double* state_matrix)
{
  int     i, j, k, index,iter, itmax, sign, nsteps, set;
  double  PI, energy, x[2*MAXATOMS], ftol, fret, angles[MAXATOMS], 
          fitness0, fitness,gcount;
  char    name[80];
  vector  Rj,Rk;
  static int ifirst=1,ran_n;
  int         ranptr,ran_i;
  static double ranlist[1000000],ran_offset[MAXWORDS],ran_mag;

  nsteps=200; /* Evolution used 50 here */
  PI = acos(-1.0);
/*----------------------------------------------------------------
/
/   do initital setup first time only
/
/----------------------------------------------------------------*/
  if(ifirst) {
    for(i=0;i<MAXBEST;i++) {
      beste[i]=1000000.0;
      bestfitness[i]=1000000.0;
    }
    read_seq();
    ifirst=0;
    
    set=1;
    if(set==1) {
      angles[0] = 0.0;
      angles[1] = 60.0;
      angles[2] = 120.0;
      angles[3] = 0.0;
      angles[4] = 0.0;
      angles[5] = 300.0;
      angles[6] = 240.0;
      angles[7] = 0.0;
      angles[8] = 0.0;
      angles[9] = 0.0;
      angles[10] = 120.0;
      angles[11] = 60.0;
      angles[12] = 0.0;
      angles[13] = 0.0;
      angles[14] = 240.0;
      angles[15] = 300.0;
      angles[16] = 0.0;
    }
    else if(set==2) {
      angles[0] = 120.0;
      angles[1] = 60.0;
      angles[2] = 60.0;
      angles[3] = 60.0;
      angles[4] = 60.0;
      angles[5] = 0.0;
      angles[6] = 60.0;
      angles[7] = 60.0;
      angles[8] = 0.0;
      angles[9] = 60.0;
      angles[10] = 0.0;
      angles[11] = 60.0;
      angles[12] = 0.0;
      angles[13] = 60.0;
      angles[14] = 0.0;
      angles[15] = 60.0;
      angles[16] = 0.0;
    }
    for(i=0;i<17;i++)
      angles[i]*=PI/180.0;
    make_polymer( angles );
    write_pdb(energy,"p0.pdb");
    j = 0;
    for(k=0;k<n_atoms;k++) {
      x[++j] = polymer[k].R.x;
      x[++j] = polymer[k].R.y;
    }
    constrain=0;
    ftol = 0.001;
    itmax = 100;
    dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func2,dfunc2);
  
    fprintf(stdout,"[steps: %4d] PERFECT ENERGY: %8.2f\n",iter,energy);
    fflush(stdout);
    write_pdb(energy,"p1.pdb");

    index=0;
    for(j=0;j<n_atoms-2;j++) {
      Rj=polymer[j].R;
      for(k=j+2;k<n_atoms;k++) {
	Rk=polymer[k].R;
	dist_matrix_0[index]=vec_dist(&Rj,&Rk);
	index++;
      }
    }
    ran_n=21; /*22*/
    ran_mag=0.4; /*0.02*/
    for(i=0;i<17*ran_n*10;i++)
      ranlist[i]=GAran(0);
  }
/*--------------------------------------------------------------------
/
/   main loop over all initial conformations
/
/--------------------------------------------------------------------*/
  fitness=0.0;
  gcount=0.0;
  ranptr=0;
  for(ran_i=0;ran_i<ran_n;ran_i++) {
    if(ran_i==0) {
      for(k=0;k<17;k++) {
	ran_offset[k]=0.0;
      }
    }
    else {
      for(k=0;k<17;k++) {
	ran_offset[k]=(ranlist[ranptr]-0.5)*ran_mag;
	ranptr++;
      }
    }

    /* make the polymer */
    sign = 1;
    for(i=0;i<n_atoms-2;i++)  {
      sign *= -1;
      angles[i] = sign * 0.1 + ran_offset[i];
    }
    make_polymer( angles );

    j = 0;
    for(k=0;k<n_atoms;k++) {
      x[++j] = polymer[k].R.x;
      x[++j] = polymer[k].R.y;
    }
    
    for(i=0;i<nsteps;i++) {
      index= 0;
      for(j=0;j<n_atoms-2;j++) {
	Rj=polymer[j].R;
	for(k=j+2;k<n_atoms;k++) {
	  Rk=polymer[k].R;
	  dist_matrix[index]=vec_dist(&Rj,&Rk)+state_matrix[index];
	  index++;
	}
      }
      
      ftol = 0.01;
      itmax = 100;
      constrain=1;
      dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func2,dfunc2);
      if(PRINT_ENERGY_1) {
	constrain=0;
	energy=e_constrain(constrain); 
	printf("Step: %5d Iterations: %5d Energy: %8.2f\n",i,iter,energy);
      }
      j=0;
      for(k=0;k<n_atoms;k++) {
	polymer[k].R.x=x[++j];
	polymer[k].R.y=x[++j];
      }
      
      fitness0 = dist_fitness();
      if(PRINT_TRAJ) {
	fprintf(stdout,"step: %4d  distance: %8.2f\n",i,fitness0);
	fflush(stdout);
	
      }
    }
    
    /* do a final minimization */
#define FINAL
#ifdef FINAL
    constrain=0;
    ftol = 0.001; 
    itmax = 100;  
    dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func2,dfunc2);
    
    index=0;
    for(j=0;j<n_atoms-2;j++) {
      Rj=polymer[j].R;
      for(k=j+2;k<n_atoms;k++) {
	Rk=polymer[k].R;
	dist_matrix[index]=vec_dist(&Rj,&Rk);
	index++;
      }
    }
#endif
    fitness0 = dist_fitness();
    fitness += fitness0;
    if(fitness0<0.1) gcount++;
    if(PRINT_ENERGY_2) {
      printf("  initial state: %4d fitness: %8.4f gcount: %6.1f\n",
	     ran_i,fitness0,gcount);
      fflush(stdout);
    }
  }
  fitness/=ran_n;
  gcount/=ran_n;

#define USEGC
#ifdef USEGC
  fitness=1.0-gcount;
#endif

  if(PRINT_ENERGY_2) {
    printf("[steps: %4d] FINAL ENERGY: %8.2f dist fitness: %8.2f\n",
	   iter,energy,fitness);
    fflush(stdout);
  }
  
  if(energy < beste[iniche]) {
    beste[iniche] = energy;
    make_name("best",iniche,".pdb",name);
    write_pdb(beste[iniche],name);
  }
  
  if(fitness < bestfitness[iniche]) {
    bestfitness[iniche] = fitness;
    make_name("bfit",iniche,".pdb",name);
    write_pdb(bestfitness[iniche],name);
  }
  
  return fitness;
}
/*-----------------------------------------------------------
/
/   check_state_table - Calculate the energy after a folding
/               sequence.
/
/   input: iniche (int) - tells which of several 
/             populations the current individual is in.
/   foldseq (int *) - 
/
/   input file: sequence.dat
/             first line is the number of atoms.
/             subsequent lines give forcefield parameters.
/             For the "standard problem", each line of 
/             sequence.dat should be:
/
/             0.50  100.0  0.50  1.0
/
/   output:
/        energy: energy after the final gradient minimization
/
/------------------------------------------------------------*/
double check_state_table(int iniche, double* state_matrix, 
			double* ran_offset,int rancnt)
{
  int     i, j, k, index,iter, itmax, sign, nsteps, set;
  double  PI, energy, x[2*MAXATOMS], ftol, fret, angles[MAXATOMS], fitness;
  char    name[80];
  vector  Rj,Rk;
  static int ifirst=1;
  FILE *fp1;
  double dm[20][20];

  nsteps=200; /* Evolution used 50 here */
  PI = acos(-1.0);
  ftol = 0.01;
  itmax = 100;

  if(ifirst) {
    for(i=0;i<MAXBEST;i++) {
      beste[i]=1000000.0;
      bestfitness[i]=1000000.0;
    }
    read_seq();
    ifirst=0;
    
    set=1;
    if(set==1) {
      angles[0] = 0.0;
      angles[1] = 60.0;
      angles[2] = 120.0;
      angles[3] = 0.0;
      angles[4] = 0.0;
      angles[5] = 300.0;
      angles[6] = 240.0;
      angles[7] = 0.0;
      angles[8] = 0.0;
      angles[9] = 0.0;
      angles[10] = 120.0;
      angles[11] = 60.0;
      angles[12] = 0.0;
      angles[13] = 0.0;
      angles[14] = 240.0;
      angles[15] = 300.0;
      angles[16] = 0.0;
    }
    else if(set==2) {
      angles[0] = 120.0;
      angles[1] = 60.0;
      angles[2] = 60.0;
      angles[3] = 60.0;
      angles[4] = 60.0;
      angles[5] = 0.0;
      angles[6] = 60.0;
      angles[7] = 60.0;
      angles[8] = 0.0;
      angles[9] = 60.0;
      angles[10] = 0.0;
      angles[11] = 60.0;
      angles[12] = 0.0;
      angles[13] = 60.0;
      angles[14] = 0.0;
      angles[15] = 60.0;
      angles[16] = 0.0;
    }
    for(i=0;i<17;i++)
      angles[i]*=PI/180.0;
    make_polymer( angles );
    write_pdb(energy,"p0.pdb");
    write_dat(energy,"p0.dat");
    ftol = 0.001; itmax = 100;  
    constrain=0;
    j = 0;
    for(k=0;k<n_atoms;k++) {
      x[++j] = polymer[k].R.x;
      x[++j] = polymer[k].R.y;
    }
    dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func2,dfunc2);
  
    fprintf(stdout,"[steps: %4d] PERFECT ENERGY: %8.2f\n",iter,energy);
    fflush(stdout);
    write_pdb(energy,"p1.pdb");
    write_dat(energy,"p1.dat");

    index=0;
    for(j=0;j<n_atoms-2;j++) {
      Rj=polymer[j].R;
      for(k=j+2;k<n_atoms;k++) {
	Rk=polymer[k].R;
	dist_matrix_0[index]=vec_dist(&Rj,&Rk);
	index++;
      }
    }
  }

#define PRINT_DMATRIX
#ifdef PRINT_DMATRIX
  for(j=0;j<n_atoms;j++) {
    for(k=0;k<n_atoms;k++) {
      dm[j][k]=0.0;
    }
  }
  index= 0;
  for(j=0;j<n_atoms-2;j++) {
    for(k=j+2;k<n_atoms;k++) {
      dm[j][k]=state_matrix[index];
      index++;
    }
  }
  make_name("dmatrix",iniche,".dat",name);
  fp1=fopen(name,"w");
  
  fprintf(fp1,"     ");
  for(k=0;k<n_atoms;k++) {
    fprintf(fp1,"  %3d",k);
  }
  fprintf(fp1,"\n");
  for(k=0;k<n_atoms+1;k++) {
    fprintf(fp1,"_____",k);
  }
  fprintf(fp1,"\n");
  for(j=0;j<n_atoms;j++) {
    fprintf(fp1,"%4d ",j);
    for(k=0;k<n_atoms;k++) {
      if(fabs(dm[j][k])>0.001) 
	fprintf(fp1,"%5.1f",dm[j][k]);
      else
	fprintf(fp1,"     ");
      
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);
#endif

  /* make the polymer */
  sign = 1;
  for(i=0;i<n_atoms-2;i++)  {
    sign *= -1;
    angles[i] = sign * 0.1 + ran_offset[i];
/*    fprintf(stdout,"angle: %8.2f  offset: %8.2f\n",sign*0.1, ran_offset[i]);*/
  }
  make_polymer( angles );
  make_name_2("x0",iniche,rancnt,".dat",name);
  write_dat(fitness,name);
 
  j = 0;
  for(k=0;k<n_atoms;k++) {
    x[++j] = polymer[k].R.x;
    x[++j] = polymer[k].R.y;
  }


  for(i=0;i<nsteps;i++) {
    index= 0;
    for(j=0;j<n_atoms-2;j++) {
      Rj=polymer[j].R;
      for(k=j+2;k<n_atoms;k++) {
	Rk=polymer[k].R;
	dist_matrix[index]=vec_dist(&Rj,&Rk)+state_matrix[index];
	index++;
      }
    }


    constrain=1;
    dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func2,dfunc2);
    if(PRINT_ENERGY_1) {
      constrain=0;
      energy=e_constrain(constrain); 
      printf("Step: %5d Iterations: %5d Energy: %8.2f\n",i,iter,energy);
    }
    j=0;
    for(k=0;k<n_atoms;k++) {
      polymer[k].R.x=x[++j];
      polymer[k].R.y=x[++j];
    }

    fitness = dist_fitness();
    if(PRINT_TRAJ) {
      fprintf(stdout,"step: %4d  distance: %8.2f\n",i,fitness);
      fflush(stdout);
    }
  }

  /* do a final minimization */
#define FINAL
#ifdef FINAL
  ftol = 0.001; itmax = 100;  
  constrain=0;
  dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func2,dfunc2);
  
  index=0;
  for(j=0;j<n_atoms-2;j++) {
    Rj=polymer[j].R;
    for(k=j+2;k<n_atoms;k++) {
      Rk=polymer[k].R;
      dist_matrix[index]=vec_dist(&Rj,&Rk);
      index++;
    }
  }
#endif
  fitness = dist_fitness();

  if(PRINT_ENERGY_2) {
    printf("[steps: %4d] FINAL ENERGY: %8.2f dist fitness: %8.2f\n",
	   iter,energy,fitness);
    fflush(stdout);
  }
  
  make_name_2("f",iniche,rancnt,".pdb",name);
  write_pdb(fitness,name);
  make_name_2("x",iniche,rancnt,".dat",name);
  write_dat(fitness,name);

  return fitness;
}
/*-----------------------------------------------------------
/
/   eval_q - Calculate the energy after a fold with charges
/
/   input: iniche (int) - tells which of several 
/             populations the current individual is in.
/   foldseq (int *) - 
/
/   input file: sequence.dat
/             first line is the number of atoms.
/             subsequent lines give forcefield parameters.
/             For the "standard problem", each line of 
/             sequence.dat should be:
/
/             0.50  100.0  0.50  1.0
/
/   output:
/        energy: energy after the final gradient minimization
/
/------------------------------------------------------------*/
double eval_q(int iniche, double* charges)
{
  int     i, j, k, index,iter, itmax, sign, nsteps;
  double  PI, energy, x[2*MAXATOMS], ftol, fret, angles[MAXATOMS], fitness;
  char    name[80];
  vector  Rj,Rk;
  static int ifirst=1;

  nsteps=200;
  PI = acos(-1.0);
  ftol = 0.01;
  itmax = 100;

  if(ifirst) {
    for(i=0;i<MAXBEST;i++) {
      beste[i]=1000000.0;
      bestfitness[i]=1000000.0;
    }
    read_seq();
    ifirst=0;
    
    angles[0] = 0.0;
    angles[1] = 60.0;
    angles[2] = 120.0;
    angles[3] = 0.0;
    angles[4] = 0.0;
    angles[5] = 300.0;
    angles[6] = 240.0;
    angles[7] = 0.0;
    angles[8] = 0.0;
    angles[9] = 0.0;
    angles[10] = 120.0;
    angles[11] = 60.0;
    angles[12] = 0.0;
    angles[13] = 0.0;
    angles[14] = 240.0;
    angles[15] = 300.0;
    angles[16] = 0.0;
    for(i=0;i<17;i++)
      angles[i]*=PI/180.0;
    make_polymer( angles );
    write_pdb(energy,"p0.pdb");
    ftol = 0.001; itmax = 100;  
    constrain=0;
    j = 0;
    for(k=0;k<n_atoms;k++) {
      x[++j] = polymer[k].R.x;
      x[++j] = polymer[k].R.y;
    }
    dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func2,dfunc2);
  
    fprintf(stdout,"[steps: %4d] PERFECT ENERGY: %8.2f\n",iter,energy);
    fflush(stdout);
    write_pdb(energy,"p1.pdb");

    index=0;
    for(j=0;j<n_atoms-2;j++) {
      Rj=polymer[j].R;
      for(k=j+2;k<n_atoms;k++) {
	Rk=polymer[k].R;
	dist_matrix_0[index]=vec_dist(&Rj,&Rk);
	index++;
      }
    }
  }

  /* make the polymer */
  sign = 1;
  for(i=0;i<n_atoms-2;i++)  {
    sign *= -1;
    angles[i] = sign * 0.1;
  }
  make_polymer( angles );

  j = 0;
  for(k=0;k<n_atoms;k++) {
    x[++j] = polymer[k].R.x;
    x[++j] = polymer[k].R.y;
    polymer[k].q=charges[k];
  }

  /* do the minimization */
  ftol = 0.001; itmax = 100;  
  constrain=1;
  dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func3,dfunc3);
  
  index=0;
  for(j=0;j<n_atoms-2;j++) {
    Rj=polymer[j].R;
    for(k=j+2;k<n_atoms;k++) {
      Rk=polymer[k].R;
      dist_matrix[index]=vec_dist(&Rj,&Rk);
      index++;
    }
  }

  fitness = dist_fitness();

  if(PRINT_ENERGY_2) {
    printf("[steps: %4d] INIT. ENERGY: %8.2f dist fitness: %8.2f\n",
	   iter,energy,fitness);
    fflush(stdout);
  }

  ftol = 0.001; itmax = 100;  
  constrain=0;
  dfpmin(x,2*n_atoms,ftol,&iter,itmax,&energy,func3,dfunc3);
  
  index=0;
  for(j=0;j<n_atoms-2;j++) {
    Rj=polymer[j].R;
    for(k=j+2;k<n_atoms;k++) {
      Rk=polymer[k].R;
      dist_matrix[index]=vec_dist(&Rj,&Rk);
      index++;
    }
  }

  fitness = dist_fitness();

  if(PRINT_ENERGY_2) {
    printf("[steps: %4d] FINAL ENERGY: %8.2f dist fitness: %8.2f\n",
	   iter,energy,fitness);
    fflush(stdout);
  }

  
  if(energy < beste[iniche]) {
    beste[iniche] = energy;
    make_name("best",iniche,".pdb",name);
    write_pdb(beste[iniche],name);
  }
  
  if(fitness < bestfitness[iniche]) {
    bestfitness[iniche] = fitness;
    make_name("bfit",iniche,".pdb",name);
    write_pdb(bestfitness[iniche],name);
  }
  
  return fitness;
}
/*-----------------------------------------------------------
/
/   dist_fitness
/
/------------------------------------------------------------*/
double dist_fitness()
{
  int index,j,k,count;
  double fitness,cutoff;
  cutoff=2.0;
  index=0;
  fitness=0;
  count=0;
  for(j=0;j<n_atoms-2;j++) {
    for(k=j+2;k<n_atoms;k++) {
      if(dist_matrix_0[index]<=cutoff) {
	fitness+=(dist_matrix[index]-dist_matrix_0[index])*
	  (dist_matrix[index]-dist_matrix_0[index]);
	count++;
      }
      index++;
    }
  }
    
  fitness = sqrt(2.0*fitness/count);
  return fitness;
}
/*-----------------------------------------------------------
/
/   e_constrain - calculate the total energy and the forces
/                 if desired, for the polymer
/
/------------------------------------------------------------*/
double e_constrain(int constrain)
{
  double e, e1, e2, e10, e20, Rbond, Kbond, Rvdw, Dvdw, rcut;
  int    i, j, k, index;
  vector a1,a2,f1,f2,Rj,Rk;

  e=e1=e2=0.0;
  rcut=4.0;

  /* zero out the forces */
  for(i=0;i<n_atoms;i++) {
    polymer[i].f.x=0.0;
    polymer[i].f.y=0.0;
    polymer[i].f.z=0.0;
  }

  /* bonding contributions */
  for(i=0;i<n_atoms-1;i++) {
    a1    = polymer[i].R;
    a2    = polymer[i+1].R;
    Rbond = polymer[i].Rbond + polymer[i+1].Rbond;
    Kbond = (polymer[i].Kbond + polymer[i+1].Kbond)/2.0;
    e10=e_bond(&a1,&a2,Rbond,Kbond,&f1,&f2);
    e1+=e10;

    polymer[i].f.x += f1.x;    
    polymer[i].f.y += f1.y;    
    polymer[i].f.z += f1.z;    
    polymer[i+1].f.x += f2.x;
    polymer[i+1].f.y += f2.y;
    polymer[i+1].f.z += f2.z;
  }

  /* special folding bonds */

  if(constrain) {
    index= 0;
    for(j=0;j<n_atoms-2;j++) {
      Rj=polymer[j].R;
      for(k=j+2;k<n_atoms;k++) {
	Rk=polymer[k].R;
	Rbond=dist_matrix[index];
	if(Rbond<rcut) {
	  Rvdw  = polymer[i].Rvdw + polymer[j].Rvdw;
	  if(Rbond<Rvdw) Rbond=Rvdw;
	  Kbond = 0.25*(polymer[j].Kbond + polymer[k].Kbond)/2.0;
	  e10=e_bond(&Rj,&Rk,Rbond,Kbond,&f1,&f2);
	  e1+=e10;
	  
	  polymer[j].f.x += f1.x;    
	  polymer[j].f.y += f1.y;    
	  polymer[j].f.z += f1.z;    
	  polymer[k].f.x += f2.x;
	  polymer[k].f.y += f2.y;
	  polymer[k].f.z += f2.z;
	}
	index++;
      }
    }
  }

  /* calculate the Van der Waals repulsion term */
  for(i=0;i<n_atoms-2;i++) {
    for(j=i+2;j<n_atoms;j++) {
      a1    = polymer[i].R;
      a2    = polymer[j].R;
      Rvdw  = polymer[i].Rvdw + polymer[j].Rvdw;
      Dvdw  = (polymer[i].Dvdw + polymer[j].Dvdw)/2.0;
      e20=e_vdw(Rvdw,Dvdw,&a1,&a2,&f1,&f2);
      e2+=e20;
      
      polymer[i].f.x += f1.x;    
      polymer[i].f.y += f1.y;    
      polymer[i].f.z += f1.z;    
      polymer[j].f.x += f2.x;
      polymer[j].f.y += f2.y;
      polymer[j].f.z += f2.z;
    }
  }
  e = e1+e2;

  return e;

}
/*-----------------------------------------------------------
/
/   e_total_q - calculate the total energy and the forces
/                 if desired, for the polymer
/
/------------------------------------------------------------*/
double e_total_q(int constrain)
{
  double e, e1, e2, e10, e20, Rbond, Kbond, Rvdw, Dvdw, rcut, q1,q2;
  int    i, j, k, index;
  vector a1,a2,f1,f2,Rj,Rk;

  e=e1=e2=0.0;
  rcut=4.0;

  /* zero out the forces */
  for(i=0;i<n_atoms;i++) {
    polymer[i].f.x=0.0;
    polymer[i].f.y=0.0;
    polymer[i].f.z=0.0;
  }

  /* bonding contributions */
  for(i=0;i<n_atoms-1;i++) {
    a1    = polymer[i].R;
    a2    = polymer[i+1].R;
    Rbond = polymer[i].Rbond + polymer[i+1].Rbond;
    Kbond = (polymer[i].Kbond + polymer[i+1].Kbond)/2.0;
    e10=e_bond(&a1,&a2,Rbond,Kbond,&f1,&f2);
    e1+=e10;

    polymer[i].f.x += f1.x;    
    polymer[i].f.y += f1.y;    
    polymer[i].f.z += f1.z;    
    polymer[i+1].f.x += f2.x;
    polymer[i+1].f.y += f2.y;
    polymer[i+1].f.z += f2.z;
  }

  /* charges */

  if(constrain) {
    for(i=0;i<n_atoms-2;i++) {
      for(j=i+2;j<n_atoms;j++) {
	a1    = polymer[i].R;
	a2    = polymer[j].R;
	q1    = polymer[i].q;
	q2    = polymer[j].q;
	e20=e_q(q1,q2,&a1,&a2,&f1,&f2);
	e2+=e20;
	
	polymer[i].f.x += f1.x;    
	polymer[i].f.y += f1.y;    
	polymer[i].f.z += f1.z;    
	polymer[j].f.x += f2.x;
	polymer[j].f.y += f2.y;
	polymer[j].f.z += f2.z;
      }
    }
  }

  /* calculate the Van der Waals repulsion term */
  for(i=0;i<n_atoms-2;i++) {
    for(j=i+2;j<n_atoms;j++) {
      a1    = polymer[i].R;
      a2    = polymer[j].R;
      Rvdw  = polymer[i].Rvdw + polymer[j].Rvdw;
      Dvdw  = (polymer[i].Dvdw + polymer[j].Dvdw)/2.0;
      e20=e_vdw(Rvdw,Dvdw,&a1,&a2,&f1,&f2);
      e2+=e20;
      
      polymer[i].f.x += f1.x;    
      polymer[i].f.y += f1.y;    
      polymer[i].f.z += f1.z;    
      polymer[j].f.x += f2.x;
      polymer[j].f.y += f2.y;
      polymer[j].f.z += f2.z;
    }
  }
  e = e1+e2;

  return e;

}




