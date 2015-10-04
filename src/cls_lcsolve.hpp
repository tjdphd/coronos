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
#include "cls_redhallmhd.hpp"
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

   private:

#ifndef HAVE_CUDA_H


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

   void setS(    std::string str_step,   stack& run, redhallmhd& physics );
   void setB(    std::string str_step,   stack& run, redhallmhd& physics );
   void setD(    std::string str_step,   stack& run, redhallmhd& physics );
   void setAi(                           stack& run, redhallmhd& physics );

   void partialsInXandY(stack& run, redhallmhd& physics, ComplexArray& U, RealArray& Ux, RealArray& Uy);
   void bracket( stack& run, redhallmhd& physics, ComplexArray& BrKt, RealArray& dx1, RealArray& dy1, RealArray& dx2, RealArray& dy2);

   double maxdU(RealArray& dx, RealArray&  dy);
   void averageAcrossLayers( stack& run, int shift_sign, RealArray& dx, RealArray&  dy);

#endif


//   void Step( std::string str_step, stack& run );

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

 void Step( std::string str_step, stack& run, redhallmhd& physics );

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
     
   public:

// lcsolve();                                    /* ~ Constructors                         ~ */
//
   lcsolve( stack& run );

   void Loop(      stack& run );                /* ~ stepping and such                     ~ */

   void passAdjacentLayers( std::string str_step, stack& run);

   ~lcsolve();                                   /* ~ Destructor                           ~ */

};

#endif
