#define WANT_STREAM
#include "opt.h"
extern double alpha;
/* Example file to demonstrate the calling sequence to a 
 * simple NLF1 function
 */
void init_rosen (int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    printf ("init_rosen: ndim != 2, ndim = %d\n", ndim);
    exit (1);
  }
  x(1) = -1.2;
  x(2) =  1.0;
}
void rosen(int mode, int n, ColumnVector& x, double& fx, ColumnVector& g)
{ // Rosenbrock's function
  double f1, f2, x1, x2;
  
  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  if (mode == 0 || mode == 3) {
    fx  = 100.* f1*f1 + f2*f2;
  }
  if (mode == 1 || mode == 3) {
    g(1) = -400.*f1*x1 - 2.*f2; 
    g(2) = 200.*f1;
  }
}
void rosen0(int n, ColumnVector& x, double& fx)
{ // Rosenbrock's function
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  fx  = 100.* f1*f1 + f2*f2;
}
void rosen2(int mode, int n, ColumnVector& x, double& fx, ColumnVector& g, 
	   SymmetricMatrix& H)
// Rosenbrock's function, n = 2 with first and second order derivatives

{ 
  int dummy = n;
  double f1, f2, x1, x2;
//  cout << "calling rosen2: mode = " << mode << "\n";

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  if (mode == 0 || mode == 3) {
    fx  = 100.* f1*f1 + f2*f2;
  }
  if (mode == 1 || mode == 3) {
    g(1) = -400.*f1*x1 - 2.*f2; 
    g(2) = 200.*f1;
  }
  
  if (mode == 2 || mode == 3) {

    f1 = (x2 - 3.0*x1*x1);
    
    H(1,1) = -400.0*f1 + 2.0;
    H(2,1) = -400.0*x1;
    H(2,2) = 200.0;
//    Print(x);
//    Print(H);
  }
}
