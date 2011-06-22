//$$ precisio.h                          floating point constants

#ifndef NEWMAT_PRECISIO_H
#define NEWMAT_PRECISIO_H

#include "include.h"
#include <limits>

#ifdef HAVE_STD
#include <cfloat>
#include <cmath>
#include <cstdlib>
#endif


namespace NEWMAT {


#ifdef HAVE_STD                 // standard library available
	
class FloatingPointPrecision
{
public:

   typedef RBD_COMMON::Real Real;

   static int Dig()              // number of decimal digits or precision
     { return std::numeric_limits<Real>::digits10 ; }

   static Real Epsilon()         // smallest number such that 1+Eps!=Eps
      { return std::numeric_limits<Real>::epsilon(); }

   static int Mantissa()         // bits in mantisa
      { return std::numeric_limits<Real>::digits; }

   static Real Maximum()         // maximum value
      { return std::numeric_limits<Real>::max(); }

   static int MaximumDecimalExponent()  // maximum decimal exponent
      { return std::numeric_limits<Real>::max_exponent10; }

   static int MaximumExponent()  // maximum binary exponent
      { return std::numeric_limits<Real>::max_exponent; }

   static Real LnMaximum()       // natural log of maximum
      { return (Real)std::log(Maximum()); }

   static Real Minimum()         // minimum positive value
      { return std::numeric_limits<Real>::min(); } 

   static int MinimumDecimalExponent() // minimum decimal exponent
      { return std::numeric_limits<Real>::min_exponent10; }

   static int MinimumExponent()  // minimum binary exponent
      { return std::numeric_limits<Real>::min_exponent; }

   static Real LnMinimum()       // natural log of minimum
      { return (Real)std::log(Minimum()); }

   static int Radix()            // exponent radix
      { return std::numeric_limits<Real>::radix; }

   static int Rounds()           // addition rounding (1 = does round)
   {
	  return std::numeric_limits<Real>::round_style ==
		 std::round_to_nearest ? 1 : 0;
   }

};


#else                              // HAVE_STD not defined

#ifndef SystemV                    // if there is float.h


#ifdef USING_FLOAT


class FloatingPointPrecision
{
public:
   static int Dig()
      { return FLT_DIG; }        // number of decimal digits or precision

   static Real Epsilon()
      { return FLT_EPSILON; }    // smallest number such that 1+Eps!=Eps

   static int Mantissa()
      { return FLT_MANT_DIG; }   // bits in mantisa

   static Real Maximum()
      { return FLT_MAX; }        // maximum value

   static int MaximumDecimalExponent()
      { return FLT_MAX_10_EXP; } // maximum decimal exponent

   static int MaximumExponent()
      { return FLT_MAX_EXP; }    // maximum binary exponent

   static Real LnMaximum()
      { return (Real)std::log(Maximum()); } // natural log of maximum

   static Real Minimum()
      { return FLT_MIN; }        // minimum positive value

   static int MinimumDecimalExponent()
      { return FLT_MIN_10_EXP; } // minimum decimal exponent

   static int MinimumExponent()
      { return FLT_MIN_EXP; }    // minimum binary exponent

   static Real LnMinimum()
      { return (Real)std::log(Minimum()); } // natural log of minimum

   static int Radix()
      { return FLT_RADIX; }      // exponent radix

   static int Rounds()
      { return FLT_ROUNDS; }     // addition rounding (1 = does round)

};

#endif                           // USING_FLOAT


#ifdef USING_DOUBLE

class FloatingPointPrecision
{
public:

   static int Dig()
      { return DBL_DIG; }        // number of decimal digits or precision

   static Real Epsilon()
      { return DBL_EPSILON; }    // smallest number such that 1+Eps!=Eps

   static int Mantissa()
      { return DBL_MANT_DIG; }   // bits in mantisa

   static Real Maximum()
      { return DBL_MAX; }        // maximum value

   static int MaximumDecimalExponent()
      { return DBL_MAX_10_EXP; } // maximum decimal exponent

   static int MaximumExponent()
      { return DBL_MAX_EXP; }    // maximum binary exponent

   static Real LnMaximum()
      { return (Real)std::log(Maximum()); } // natural log of maximum

   static Real Minimum()
   {
//#ifdef __BCPLUSPLUS__
//       return 2.225074e-308;     // minimum positive value
//#else
       return DBL_MIN;
//#endif
   }

   static int MinimumDecimalExponent()
      { return DBL_MIN_10_EXP; } // minimum decimal exponent

   static int MinimumExponent()
      { return DBL_MIN_EXP; }    // minimum binary exponent

   static Real LnMinimum()
      { return (Real)std::log(Minimum()); } // natural log of minimum


   static int Radix()
      { return FLT_RADIX; }      // exponent radix

   static int Rounds()
      { return FLT_ROUNDS; }     // addition rounding (1 = does round)

};

#endif                             // USING_DOUBLE

#else                              // if there is no float.h

#ifdef USING_FLOAT

class FloatingPointPrecision
{
public:

   static Real Epsilon()
      { return std::pow(2.0,(int)(1-FSIGNIF)); }
                                   // smallest number such that 1+Eps!=Eps

   static Real Maximum()
      { return MAXFLOAT; }            // maximum value

   static Real LnMaximum()
      { return (Real)std::log(Maximum()); }  // natural log of maximum

   static Real Minimum()
      { return MINFLOAT; }             // minimum positive value

   static Real LnMinimum()
      { return (Real)std::log(Minimum()); }  // natural log of minimum

};

#endif                                  // USING_FLOAT


#ifdef USING_DOUBLE

class FloatingPointPrecision
{
public:

   static Real Epsilon()
      { return std::pow(2.0,(int)(1-DSIGNIF)); }
                                      // smallest number such that 1+Eps!=Eps

   static Real Maximum()
      { return MAXDOUBLE; }           // maximum value

   static Real LnMaximum()
      { return LN_MAXDOUBLE; }        // natural log of maximum

#ifdef IEEE
   static Real Minimum()
      { return MINDOUBLE; }
#else
   static Real Minimum()
      { return DBL_MIN; }
#endif

   static Real LnMinimum()
      { return LN_MINDOUBLE; }        // natural log of minimum
};

#endif                                // USING_DOUBLE

#endif                                // SystemV

#endif                                // HAVE_STD


}  // NEWMAT namespace



#endif                                // PRECISION_LIB
