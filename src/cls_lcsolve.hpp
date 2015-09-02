/* class lcsolve (definition)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * solver class for longcope solver designed to 
 * work on canvas class type "stack"
 *
 */

#ifndef CLS_LCSOLVE
#define CLS_LCSOLVE

#include "mpi.h"
#include "cls_stack.hpp"
//#include "cls_redhallmhd.hpp"
#include "nsp_constants.hpp"

#include <assert.h>
#include <fstream>
#include <complex>
#include <vector>
#include <cstddef>
#include<iomanip>

#ifdef HAVE_CUDA_H
#include "cls_lcsolve_cuda_ext.hpp"
// #else
// #include<fftw3.h>
#endif

using namespace constants;

class lcsolve
{

   friend class redhallmhd;
//   friend class fft;

   private:

#ifndef HAVE_CUDA_H

   ComplexArray tU0; /* ~ for holding predictor results?       ~ */
   ComplexArray tU1;
   ComplexArray tU2;
   ComplexArray tU3;

   RealArray    SE0; /* ~ same for different layers of U's      ~ */
   RealArray    SE1;
   RealArray    SE2;
   RealArray    SE3;

   RealArray    SI0; /* ~ same for different layers of U's      ~ */
   RealArray    SI1;
   RealArray    SI2;
   RealArray    SI3;

   ComplexArray B0; /* ~ different for different layers of U's ~ */
   ComplexArray B1;
   ComplexArray B2;
   ComplexArray B3;

   ComplexArray D0; /* ~ different for different layers of U's ~ */
   ComplexArray D1;
   ComplexArray D2;
   ComplexArray D3;

   ComplexArray A0; /* ~ different for different layers of U's ~ */
   ComplexArray A1;
   ComplexArray A2;
   ComplexArray A3;

   void createFields(stack& run );
   void destroyFields();

   void partialsInXandY(ComplexArray& U, RealArray& Ux, RealArray& Uy);
   void bracket( ComplexArray& BrKt, RealArray& dx1, RealArray& dy1, RealArray& dx2, RealArray& dy2);

   double maxdU(RealArray& dx, RealArray&  dy);
   void averageAcrossLayers( int shift_sign, RealArray& dx, RealArray&  dy);

#endif


   void Step( std::string str_step );
     
   public:

// lcsolve();                                    /* ~ Constructors                         ~ */
//
   lcsolve( stack& run );

   void Loop(      stack& run );                /* ~ stepping and such                     ~ */

   void passAdjacentLayers( std::string str_step, stack& run);

   ~lcsolve();                                   /* ~ Destructor                           ~ */

};

#endif
