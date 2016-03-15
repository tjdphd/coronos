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
 * FILE: Definition of class "redhallmhd"
 *
 * DESCRIPTION: For defining and implementing the reduced mhd Hall physics for
 *              coronos. this class is responsible for "filling" lcsolve's data
 *              structures with the appropriate values - based on its
 *              "knowledge" of the physical model of * the plasma. These values
 *              are needed by lcsolve so that lcsolve can update its stack.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

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

#ifdef HAVE_CUDA_H
  #include "cls_redhallmhd_cuda_ext.hpp"
#endif

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

  void initIRMHD(                       stack& run);
 
  void initialize(                      stack& run );
  void OfromP(                          stack& run );                 /* ~ Obtain vorticity from P                   ~ */
  void HfromA(                          stack& run );                 /* ~ Obtain H from A                           ~ */
  void JfromA(                          stack& run );                 /* ~ Obtain J from A                           ~ */


  void applyFootPointDrivingBC(         stack& run );                 /* ~ "pevol"                                   ~ */
  void applyLineTiedBC( std::string str_step, stack& run );           /* ~ pbot and p(:,n3) set to zero              ~ */

//  void finalizeFootPointDriving(        stack& run, lcsolve& solve);  /* ~                                         ~ */

  public:

  parameter_map physics_data;

  fft fftw;

  ComplexArray P;                                                     /* ~ Stream function Phi in Fourier Space      ~ */
  ComplexArray A;                                                     /* ~ flux function A in Fourier Space          ~ */
  ComplexArray J;                                                     /* ~ current density in Fourier Space          ~ */

  RealArray valfven;                                                  /* ~ Needed for Inhomogeneous RMHD             ~ */
  RealArray dvalfdz;                                                  /* ~ Needed for Inhomogeneous RMHD             ~ */
  RealArray nofz;                                                     /* ~ Needed for Inhomogeneous RMHD             ~ */
  RealArray dndz;                                                     /* ~ Needed for Inhomogeneous RMHD             ~ */
  RealArray umean;

  RealArray Elln;
  RealArray EllA;
  RealArray EllB;

  RealArray h11;
  RealArray h12;
  RealArray h21;
  RealArray h22;

  RealArray maxU;                                                     /* ~ for time-step determination               ~ */

  void updatePAJ( std::string str_step, stack& run );
  void applyBC(   std::string str_step, stack& run );                 /* ~ Apply Boundary Conditions at current step ~ */
  void updateTimeInc(                   stack& run );

  void PfromO(                          stack& run );                 /* ~ Obtain P from vorticity                   ~ */
  void AfromH(                          stack& run );                 /* ~ Obtain A from H                           ~ */

  void evalElls(                        stack& run );                 /* ~ calculate l's and h's at each layer       ~ */
  void evalValf(                        stack& run );                 /* ~ calculate Va at each layer                ~ */
  void evalUmean(                       stack& run );                 /* ~ calculate Va at each layer                ~ */


  void trackEnergies(int l, int nw,     stack& run );                 /* ~ update energy quantities between steps    ~ */
  void evalTotalKineticEnergy(          stack& run, int i_pe);
  void evalTotalMagneticEnergy(         stack& run, int i_me);
  void evalTotalVorticitySqd(           stack& run, int i_oe);
  void evalTotalCurrentSqd(             stack& run, int i_ce);
  void evalTotalFootPointKE(            stack& run, int i_fp);        /* ~ Misnomer?                                 ~ */
  void evalTotalPoyntingFlux(           stack& run, int i_fe);        /* ~ Poynting Flux                             ~ */

  void physicsFinalize(                 stack& run );                 /* ~ end of subrun bookkeeping                 ~ */
                                                                      /* ~ this is just a stub right now, but will   ~ */
                                                                      /* ~ eventually make it possible to privatize  ~ */
                                                                      /* ~ some of the currently public procedures   ~ */

  redhallmhd();                                                       /* ~ Constructor (default)                     ~ */

  redhallmhd(                           stack& run );                 /* ~ Constructor                               ~ */

  ~redhallmhd();                                                      /* ~ Destructor                                ~ */

};

#endif
