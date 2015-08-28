/* class stack (definition)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * child class of canvas contains conmplete information about problem
 * and solution in its parameter maps and defines the primary data
 * members 
 *
 */

#ifndef CLS_STACK
#define CLS_STACK

#include "cls_canvas.hpp"
// #include "cls_lcsolve.hpp"
#include "nsp_constants.hpp"
//#include "nsp_utilities.hpp"

#include <cmath>
#include <fstream>
#include<new>

#include<fftw3.h>
#include <assert.h>

#include "mpi.h"

using namespace constants;

class stack : public canvas
{

  private:

  void init_stack_data();           /* ~ gather/infer information to be
                                         included in stack_data container        ~ */
  public:

  stack();                          /* ~ Constructors                            ~ */
  stack(std::string coronos_in);

  parameter_map stack_data;

  InputOutputArray U;               /* ~ raw input array                         ~ */

  RealArray x;                      /* ~ For holding x-coordinates               ~ */
  RealArray y;                      /* ~ For holding y-coordinates               ~ */
  RealArray z;                      /* ~ For holding z-coordinates               ~ */

  void   allocUi();                 /* ~ Allocators/De-allocators                ~ */
  void deallocUi();
  void     zeroU();                 /* ~ a convenience function                  ~ */

#ifndef HAVE_CUDA_H

  RealArray rt;                     /* ~ Fourier Transform Related               ~ */

  void rtInit();
  void rtFree();

  RealArray kx;
  RealArray ky;
  RealArray k2;
  RealArray inv_k2;

  void kInit( );                     /* ~ their initializers and "disposers"     ~ */
  void kFree( );

/* ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft  ~ */

//   double          * r_in;                       /* ~ input and output data arguments      ~ */
//   double          * r_out;                      /* ~ to FFT routines.                     ~ */
//
//   std::complex<double> * cplx_in;
//   std::complex<double> * cplx_out;
//
//   fftw_plan      p_lay_for;                     /* ~ For establishing plans for forward   ~ */
//   fftw_plan      p_lay_rev;                     /* ~ reverse FFT's of layers              ~ */
//
//   void fftwInitialize();                        /* ~ For allocating and deallocating "in" ~ */
//   void fftwFinalize();                          /* ~ and "out" arguments of FFT's, and    ~ */
//                                                 /* ~ for initializing and "destroying"    ~ */
//                                                 /* ~ FFT plans.                           ~ */
//
//   void fftwForwardAll( lcsolve& solve);         /* ~ Forward FFT all fields all layers    ~ */
//   void fftwReverseAll( lcsolve& solve);         /* ~ Reverst FFT all fields all layers    ~ */
//
//   void fftwForwardLayerofField ( lcsolve& solve, int layer, int field );
//   void fftwReverseLayerofField ( lcsolve& solve, int layer, int field );
//
//   void fftwForwardRaw( RealArray&    Rin, ComplexArray& Cout);
//   void fftwReverseRaw( ComplexArray& Cin, RealArray&    Rout);

/* ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft  ~ */

#endif

  std::string getLastDataFilename(int srun); /* ~ Input/Output                            ~ */
  std::string getNextDataFilename();

  void initxy();                     /* ~ Calculate x-and y-coordinates of layers ~ */
  void initz();                      /* ~ Calculate z-coordinates of layers       ~ */

  ~stack();                          /* ~ Destructor                              ~ */

};

#endif
