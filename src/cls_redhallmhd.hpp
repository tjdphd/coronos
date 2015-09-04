/* class redhallmhd (definition)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * For defining and implementing the reduced mhd Hall physics for coronos.
 * this class is responsible for "filling" lcsolve's data structures with the 
 * appropriate values - based on its "knowledge" of the physical model of
 * the plasma. These values are needed by lcsolve so that lcsolve can update 
 * its stack.
 *
 */

#ifndef CLS_REDHALLMHD
#define CLS_REDHALLMHD

#include "mpi.h"
#include "nsp_constants.hpp"
#include "cls_parameter_map.hpp"
#include "cls_stack.hpp"
#include "cls_fft.hpp"
#include <assert.h>
#include <fstream>
#include <complex>
#include <vector>
#include <cstddef>
#include<iomanip>

class redhallmhd
{

  private:






  ComplexArray roldlb;
  ComplexArray rnewlb;

  ComplexArray roldub;
  ComplexArray rnewub;

  void init_physics_data(               stack& run);

  void initTimeInc(                     stack& run);
  void initU(                           stack& run);                  /* ~ U initialization functions                ~ */
  void computeFourierU(                 stack& run);
  void computeRealU(                    stack& run);
  void readUData(                       stack& run);

  void pLinzEnv(                        stack& run);
  void initBoundaries(                  stack& run);
  void countModes(                      stack& run);
  void initFootPointDriving(            stack& run);
  void initNoDrive(                     stack& run);

  void initialize(                      stack& run );
  void OfromP(                          stack& run );                 /* ~ Obtain vorticity from P                   ~ */
  void HfromA(                          stack& run );                 /* ~ Obtain H from A                           ~ */

  void PfromO(                          stack& run );                 /* ~ Obtain P from vorticity                   ~ */
  void AfromH(                          stack& run );                 /* ~ Obtain A from H                           ~ */

  void applyFootPointDrivingBC(         stack& run );                 /* ~ "pevol"                                   ~ */
  void applyLineTiedBC(                 stack& run );                 /* ~ pbot and p(:,n3) set to zero              ~ */

//  void finalizeFootPointDriving(        stack& run, lcsolve& solve);  /* ~                                         ~ */

  public:

  parameter_map physics_data;

  fft fftw;

  ComplexArray P;                                                     /* ~ Stream function Phi in Fourier Space ~ */
  ComplexArray A;                                                     /* ~ flux function A in Fourier Space     ~ */
  ComplexArray J;                                                     /* ~ current density in Fourier Space     ~ */

  RealArray maxU;                                                     /* ~ for time-step determination          ~ */

  void updatePAJ( std::string str_step, stack& run );
  void applyBC(                         stack& run );                 /* ~ Apply Boundary Conditions at current step ~ */
  void updateTimeInc(                   stack& run                );
  void writeUData(                      stack& run);

//  void finalize(                        stack& run, lcsolve& solve);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  redhallmhd();                                                       /* ~ Constructor (default)                     ~ */
  redhallmhd(                           stack& run );                 /* ~ Constructor                               ~ */

  ~redhallmhd();                                                      /* ~ Destructor                                ~ */

};

#endif
