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
 * FILE: definition of namespace "constants"
 *
 * DESCRIPTION: To be added prior to public-release
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef NSP_CONSTANTS
#define NSP_CONSTANTS

#include "../config.h"
#include <cmath>
#include<vector>
#include<complex>

namespace constants {

#ifdef LD_PRECISION_H

typedef long double ***                         InputOutputArray;

typedef std::vector<long double>                RealArray;
typedef std::vector< RealArray >                Real2DArray;
typedef std::vector<Real2DArray >               Real3DArray;
typedef std::vector<std::complex<long double> > ComplexArray;

typedef long double                             RealVar;
typedef std::complex<long double>               ComplexVar;

   static const long double pi                  =  3.14159265358979323846264338327950288L;
   static const long double zero                =  0.0L;
   static const long double one                 =  1.0L;
   static const long double two                 =  2.0L;
   static const long double three               =  3.0L;
   static const long double four                =  4.0L;
   static const long double eight               =  8.0L;
   static const long double nine                =  9.0L;
   static const long double ten                 = 10.0L;
   static const long double half                =  0.5L;
   static const long double two_pi              =  two*pi;
   static const long double two_thirds          =  two / three;
   static const long double smallish            =  1.0e-12L;
   static const long double biggish             =  1.0e+16L;
   static const long double tiny                =  1.0e-16L;
   static const long double teensy              =  1.0e-30L;
   static const long double huge                =  1.0e+30L;

   static const std::complex<long double> iunit = std::complex<long double>(zero, one );
   static const std::complex<long double> czero = std::complex<long double>(zero, zero);
   static const std::complex<long double> chuge = std::complex<long double>(huge, huge);
   static const std::complex<long double> cone  = std::complex<long double>(one,  zero);
   static const std::complex<long double> ctwo  = std::complex<long double>(two,  zero);
   static const std::complex<long double> chalf = std::complex<long double>(half, zero);

#elif defined OD_PRECISION_H

typedef double ***                              InputOutputArray;
typedef std::vector<double>                     RealArray;
typedef std::vector< RealArray >                Real2DArray;
typedef std::vector< Real2DArray >              Real3DArray;
typedef std::vector<std::complex<double> >      ComplexArray;

typedef double                                  RealVar;
typedef std::complex<double>                    ComplexVar;

   static const double pi                  =  3.14159265358979323846264338327950288L;
   static const double zero                =  0.0L;
   static const double one                 =  1.0L;
   static const double two                 =  2.0L;
   static const double three               =  3.0L;
   static const double four                =  4.0L;
   static const double eight               =  8.0L;
   static const double nine                =  9.0L;
   static const double ten                 = 10.0L;
   static const double half                =  0.5L;
   static const double two_pi              =  two*pi;
   static const double two_thirds          =  two / three;
   static const double smallish            =  1.0e-12L;
   static const double biggish             =  1.0e+16L;
   static const double tiny                =  1.0e-16L;
   static const double teensy              =  1.0e-30L;
   static const double huge                =  1.0e+30L;

   static const std::complex<double> iunit = std::complex<double>(zero, one );
   static const std::complex<double> czero = std::complex<double>(zero, zero);
   static const std::complex<double> chuge = std::complex<double>(huge, huge);
   static const std::complex<double> cone  = std::complex<double>(one,  zero);
   static const std::complex<double> ctwo  = std::complex<double>(two,  zero);
   static const std::complex<double> chalf = std::complex<double>(half, zero);

#endif

}

#endif
