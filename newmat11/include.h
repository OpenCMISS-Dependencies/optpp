//$$ include.h           include files required by various versions of C++

#ifndef NEWMAT_INCLUDE_H
#define NEWMAT_INCLUDE_H

#ifdef HAVE_OPTPP_CONFIG_H
#include "OPT++_config.h"
#endif

#define use_namespace                   // define name spaces

//#define SETUP_C_SUBSCRIPTS              // allow element access via A[i][j]

#define OPT_COMPATIBLE                  // for use with opt++

// Activate just one of the following 3 statements

//#define SimulateExceptions              // use simulated exceptions
//#define UseExceptions                   // use C++ exceptions
#define DisableExceptions               // do not use exceptions

#define TEMPS_DESTROYED_QUICKLY         // for compilers that delete
					// temporaries too quickly

//#define TEMPS_DESTROYED_QUICKLY_R       // the same thing but applied
					// to return from functions only

//#define DO_FREE_CHECK                   // check news and deletes balance

#define USING_DOUBLE                    // elements of type double
//#define USING_FLOAT                   // elements of type float

#define bool_LIB 0                      // for compatibility with my older libraries

//#define ios_format_flags ios::fmtflags  // for Gnu 3 and Intel for Linux

//#define use_float_h                   // use float.h for precision data

//#define HAS_INT64                     // if unsigned _int64 is recognised
                                        // used by newran03
                                        
// comment out next line if Exception causes a problem
#define TypeDefException

//*********************** end of options set by user ********************


#define ios_format_flags ios::fmtflags

#ifdef HAVE_STD                       // using standard library
   #include <cstdlib>
   #ifdef _MSC_VER
      #include <limits>                 // for VC++6
   #endif
   #ifdef WANT_STREAM
      #include <iostream>
      #include <iomanip>
   #endif
   #ifdef WANT_MATH
      #include <cmath>
   #endif
   #ifdef WANT_STRING
      #include <cstring>
   #endif
   #ifdef WANT_TIME
      #include <ctime>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream>
   #endif

#else

// BMA 11/18/2010: Not cleaning the following up since we expect to
// almost always be in the STD block

#define DEFAULT_HEADER                  // use AT&T style header
                                        // if no other compiler is recognised

#ifdef _MSC_VER                         // Microsoft
   #include <stdlib.h>

//   reactivate these statements to run under MSC version 7.0
//   typedef int jmp_buf[9];
//   extern "C"
//   {
//      int __cdecl setjmp(jmp_buf);
//      void __cdecl longjmp(jmp_buf, int);
//   }

   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif

#ifdef __ZTC__                          // Zortech
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.hpp>
      #include <iomanip.hpp>
      #define flush ""                  // not defined in iomanip?
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif

#if defined __BCPLUSPLUS__ || defined __TURBOC__  // Borland or Turbo
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>            // Borland has both float and values
                                    // but values.h returns +INF for
                                    // MAXDOUBLE in BC5
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif

#ifdef __GNUG__                         // Gnu C++
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif

#ifdef __PGI                         // Gnu C++
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif

#ifdef __WATCOMC__                      // Watcom C/C++
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif


#ifdef macintosh                        // MPW C++ on the Mac
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <float.h>
      #include <math.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif

#ifdef use_float_h                      // use float.h for precision values
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <float.h>
      #include <math.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
   #undef DEFAULT_HEADER
#endif


#ifdef DEFAULT_HEADER                   // for example AT&T
   #define ATandT
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #define SystemV                   // use System V
      #include <values.h>
   #endif
   #ifdef WANT_STRING
      #include <string.h>
   #endif
   #ifdef WANT_TIME
      #include <time.h>
   #endif
   #ifdef WANT_FSTREAM
      #include <fstream.h>
   #endif
#endif                                  // DEFAULT_HEADER

#endif                                  // HAVE_STD: else


namespace RBD_COMMON {

#ifdef USING_FLOAT                      // set precision type to float
typedef float Real;
typedef double long_Real;
#endif

#ifdef USING_DOUBLE                     // set precision type to double
typedef double Real;
typedef long double long_Real;
#endif

}


namespace RBD_COMMON {}
namespace RBD_LIBRARIES                 // access all my libraries
{
   using namespace RBD_COMMON;
}


#endif  // NEWMAT_INCLUDE_H
