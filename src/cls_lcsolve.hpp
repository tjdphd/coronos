/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 *
 * CORONOS||SONOROC - Version 0.1
 *
 * (S)ynthesized  (O)bject-based (N)umerical (O)bservatory for (R)HMHD [also RMHD and IRHMD] with (O)ptional (C)UDA-acceleration
 *
 * AUTHOR: Timothy J. Dennis
 *         tdennis10@alaska.edu
 *
 * CONTRIBUTORS:
 *
 *         C. S. Ng
 *         LiWei Lin
 *         Others to be included prior to public release
 *
 * copyright 2014-2016 
 *
 * Space Physics and Aeronomy
 * Geophysical Institute
 * University of Alaska, Fairbanks
 *
 * All Rights Reserved.
 *
 * This version of the code is pre-public release.
 * Please contact the author if you are not certain
 * you have an up-to-date working copy.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* 
 *        FILE: Definition of class "lcsolve"
 *
 * DESCRIPTION: solver class for longcope solver designed to work on canvas 
 *              class type "stack"
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

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

   RealVar maxdU(RealArray& dx, RealArray&  dy, int i_grid, int i_layers);
   void averageAcrossLayers( stack& run, int shift_sign, RealArray& dx, RealArray&  dy);


   void Step( std::string str_step, stack& run );

#endif

   public:

// lcsolve();                                    /* ~ Constructors                         ~ */

   lcsolve( stack& run );
   void Loop(      stack& run );                /* ~ stepping and such                     ~ */

#ifndef HAVE_CUDA_H

   void passAdjacentLayers( std::string str_step, stack& run);

#endif

   ~lcsolve();                                   /* ~ Destructor                           ~ */

};

#endif
