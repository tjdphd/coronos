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
 *        FILE: definition of class "stack"
 *
 * DESCRIPTION: child class of canvas contains conmplete information about 
 *              problem and solution in its parameter maps and defines the 
 *              primary data members 
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef CLS_STACK
#define CLS_STACK

#include "nsp_constants.hpp"
#include "cls_canvas.hpp"
#include "mpi.h"

#ifdef HAVE_CUDA_H
  #include "cls_stack_cuda_ext.hpp"
#endif

using namespace constants;

class stack : public canvas
{

  private:

  void init_stack_data();                     /* ~ gather/infer information to be
                                                   included in stack_data container        ~ */
  public:

  stack();                                    /* ~ Constructors                            ~ */
  stack(std::string coronos_in);

  parameter_map stack_data;

  InputOutputArray U;                         /* ~ raw input/output array                  ~ */
  InputOutputArray AUX;                       /* ~ raw output array for auxiliary fields   ~ */

  RealArray x;                                /* ~ For holding x-coordinates               ~ */
  RealArray y;                                /* ~ For holding y-coordinates               ~ */
  RealArray z;                                /* ~ For holding z-coordinates               ~ */

  void   allocUi();                           /* ~ Allocators/De-allocators                ~ */
  void deallocUi();

  void   allocAUX();                          /* ~ Allocators/De-allocators                ~ */
  void deallocAUX();

  void     zeroU();                           /* ~ a convenience function                  ~ */

#ifndef HAVE_CUDA_H

   ComplexArray U0;                           /* ~Fourier Space Field Arrays               ~ */
   ComplexArray U1;
   ComplexArray U2;
   ComplexArray U3;

   ComplexArray tU0;                          /* ~ for holding predictor results           ~ */
   ComplexArray tU1;
   ComplexArray tU2;
   ComplexArray tU3;

  RealArray kx;                               /* ~ Fourier Space Related                   ~ */
  RealArray ky;
  RealArray k2;
  RealArray inv_k2;

#endif

  void initAUX();                             /* ~ for containing auxiliary field data     ~ */

  void writeUData();                          /* ~ Input/Output                            ~ */


  std::string getLastDataFilename(int srun);
  std::string getNextDataFilename();
  void writeParameters();
  void writeParameters(int srun);

  void initxyz();                             /* ~ Calculate x-and y-coordinates of layers ~ */

  ~stack();                                   /* ~ Destructor                              ~ */

};

#endif
