#ifndef pga_h
#define pga_h
#define DEBUG     0
#define MAXLOCI   2000
#define MAXPOP    200
#define MAXBEST   100
#define MAXWORDS  200
#define EVAL      1
#define PRINTALL  0
#define MAX_NI    10

#define STEP_FTN  0
#define ROULETTE  1

typedef struct {
  int        chromosome[MAXLOCI];
  double     fitness;
  int        select_flag;
  int        niche;
  int        duplicate;
} individual;

typedef struct {
  int        npop;
  individual indiv[MAXPOP];
} population;

typedef struct {
  int nbest;
  individual indiv[MAXBEST];
} bestpop;

typedef struct {
  int    number;
  double   fitness;
} best_sort_struct;

typedef struct {
  char     eval[80];
  int      restart;
  int      niches;
  int      npop;
  int      ngen;
  int      chr_length;
  int      nbest;
  int      save_freq;
  int      print_flag;
  int      niflag;
  int      MSEED;
  int      ni_gen[MAX_NI];
  double   select_rate;
  int      select_method;
  double   mutate_rate;
  double   threshold;
  double   crossover_rate;
  int      words;
  int      start[MAXWORDS];
  int      length[MAXWORDS];
  int      type[MAXWORDS];
  double   min[MAXWORDS];
  double   max[MAXWORDS];
} parm_struct;

#ifdef __cplusplus
extern "C" {
#endif
void read_GAparms(parm_struct*,FILE *);
void write_GAparms(parm_struct*,FILE *);
void print_best(int,parm_struct*,bestpop*);
void print_pop(int,parm_struct*,population*);
void statistics_pop(int,int,parm_struct*,population*,bestpop*,FILE*);
void evaluate_pop(int,parm_struct*,population*);
double evaluate_individual(int,parm_struct*,int*);
double eval_user(int,char*,int,double*);
void make_name(char*,int,char*,char*);
void make_name_2(char*,int,int,char*,char*);
void filter(parm_struct*,int*,double*);
double str_to_double(int*,int,int,double,double);
int str_to_int(int*,int,int,double,double);
void initialize_pop(parm_struct*,population*,int);
void initialize_best(parm_struct*, bestpop*,int);
void mutate_pop(parm_struct*,population*);
void best_pop(population*,bestpop*,bestpop*);
void select_pop_stepftn(parm_struct*,population*);
void select_pop_roulette(parm_struct*,population*);
void crossover_pop(parm_struct*,population*,population*,int);
void crossover(parm_struct*,individual*,individual*,individual*,individual*);
double hamming_distance(parm_struct*,individual*,individual*);
double GAran(int);
double ran3(int*);
void degray(int*,int*,int,int);
/*int max(int,int);*/
/*int min(int,int);*/
double rmax(double,double);
double rmin(double,double);
int compare_bsort();
int compare_pop();
void GA_setup(parm_struct*,population*,population*,bestpop*,bestpop*,FILE*);
void GA_loop(parm_struct*,population*,population*,bestpop*,bestpop*,FILE*,int);
void GA_main();
#ifdef __cplusplus
}
#endif
#endif
