/*-----------------------------------------------------------
/
/   polymer.h
/
/------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define RELAX              0
#define PRINTALL           0
#define PRINT_ENERGY_1     0
#define PRINT_ENERGY_2     1
#define PRINT_VEC          0
#define PRINT_FORCES       0
#define PRINT_TRAJ         0
#define MAXATOMS           62
#define HISTORY            0
#define DEBUGEVAL          0

typedef struct {
  vector R;
  vector f;
  double Rbond;
  double Kbond;
  double Rvdw;
  double Dvdw;
  double q;
} atom_struct;

double e_total(int,int);
double e_constrain(int);
double e_total_q(int);
double eval_fold(int,double*);
double eval_cart(int,double*);
double eval_state_table(int,double*);
double check_state_table(int,double*,double*,int);
double eval_q(int,double*);
void write_mm(double,char*);
void write_pdb(double,char*);
void write_dat(double,char*);
double func(double*);
void dfunc(double*,double*);
double func2(double*);
void dfunc2(double*,double*);
double func3(double*);
void dfunc3(double*,double*);
void make_polymer(double*);
double distance(atom_struct*,atom_struct*);
double cos_angle(atom_struct*,atom_struct*,atom_struct*);
double e_vdw(double,double,vector*,vector*,vector*,vector*);
double e_q(double,double,vector*,vector*,vector*,vector*);
void read_seq();
void print_atom(atom_struct*);
double dist_fitness();





