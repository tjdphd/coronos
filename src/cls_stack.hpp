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

#include "nsp_constants.hpp"
#include "cls_canvas.hpp"
#include "mpi.h"

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

  InputOutputArray U;                         /* ~ raw input array                         ~ */

  RealArray x;                                /* ~ For holding x-coordinates               ~ */
  RealArray y;                                /* ~ For holding y-coordinates               ~ */
  RealArray z;                                /* ~ For holding z-coordinates               ~ */

  void   allocUi();                           /* ~ Allocators/De-allocators                ~ */
  void deallocUi();
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

  void writeUData();                          /* ~ Input/Output                            ~ */
  std::string getLastDataFilename(int srun);
  std::string getNextDataFilename();
  void writeParameters();
  void writeParameters(int srun);

  void initxyz();                             /* ~ Calculate x-and y-coordinates of layers ~ */

  ~stack();                                   /* ~ Destructor                              ~ */

};

#endif
