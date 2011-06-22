#ifndef optga_h
#define optga_h

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "opt.h"
#include "pga.h"	

//----------------------------------------------------------------------
// GA Method objects
//----------------------------------------------------------------------

class OptGA: public OptDirect {

protected:
  NLP0 *nlp;
  parm_struct *ga_parms;
  FILE *infile;

public:
  OptGA(){}
  OptGA(NLP0* p): OptDirect(p->GetDim()), nlp(p) 
    {strcpy(method,"GA");}
  OptGA(NLP0* p, TOLS t): OptDirect(p->GetDim(), t), nlp(p) 
    {strcpy(method,"GA");}

  virtual ~OptGA(){}
  virtual void AcceptStep(int k, int step_type)
    {DefaultAcceptStep(k, step_type);}

  void Set_Parms(parm_struct *p)       {ga_parms = p;}
  parm_struct *Get_Parms()         {return ga_parms;}

  virtual void UpdateModel(int k)
    {DefaultUpdateModel(k);}

  void SetInput(FILE *fp) { infile = fp;}
  FILE *GetInput() { return infile;}

//
// These are defined elsewhere

  void InitOpt();
  void Optimize();
  int CheckConvg();
  void PrintStatus(char *);
};

//
//  Derived from NLP's
//
class GANLF0: public NLP0 {
protected:
  USERFCN0 fcn;
  INITFCN init_fcn;
  int init_flag;
  parm_struct *ga_parms;
  int iniche;
  
public:
  GANLF0() {}  // Constructor
  virtual ~GANLF0() {}               // Destructor
  GANLF0(int ndim):NLP0(ndim){}
  GANLF0(int ndim, USERFCN0 f):NLP0(ndim){fcn = f;}
  GANLF0(int ndim, USERFCN0 f, INITFCN i):
      NLP0(ndim){fcn = f; init_fcn = i; init_flag = NO;}

  void Set_Parms(parm_struct *p)       {ga_parms = p;}
  parm_struct *Get_Parms()         {return ga_parms;}

  void Set_Iniche(int i)       {iniche = i;}
  int  Get_Iniche()            {return iniche;}


  void InitFcn();               // Initialize selected function
  real EvalF();               // Evaluate f 
  real EvalF(ColumnVector& x); // Evaluate f and return the value

};

#endif














