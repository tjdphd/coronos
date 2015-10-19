#ifndef NSP_CONSTANTS
#define NSP_CONSTANTS

#include<config.h>
#include<vector>
#include<complex>

namespace constants {

#ifdef LD_PRECISION_H

typedef long double ***                         InputOutputArray;
typedef std::vector<long double>                RealArray;
typedef std::vector<std::complex<long double> > ComplexArray;

typedef long double                             RealVar;
typedef std::complex<long double>               ComplexVar;

   static const long double pi                  = 3.14159265358979323846264338327950288L;
   static const long double zero                = 0.0L;
   static const long double one                 = 1.0L;
   static const long double two                 = 2.0L;
   static const long double three               = 3.0L;
   static const long double four                = 4.0L;
   static const long double half                = 0.5L;
   static const long double two_pi              = two*pi;
   static const long double two_thirds          = two / three;
   static const long double teensy              = 1.0e-100L;

   static const std::complex<long double> iunit = std::complex<long double>(zero, one );
   static const std::complex<long double> czero = std::complex<long double>(zero, zero);
   static const std::complex<long double> cone  = std::complex<long double>(one,  zero);
   static const std::complex<long double> ctwo  = std::complex<long double>(two,  zero);
   static const std::complex<long double> chalf = std::complex<long double>(half, zero);

#elif defined OD_PRECISION_H

typedef double ***                         InputOutputArray;
typedef std::vector<double>                RealArray;
typedef std::vector<std::complex<double> > ComplexArray;

typedef double                             RealVar;
typedef std::complex<double>               ComplexVar;

   static const double pi                  = 3.14159265358979323846264338327950288L;
   static const double zero                = 0.0L;
   static const double one                 = 1.0L;
   static const double two                 = 2.0L;
   static const double three               = 3.0L;
   static const double four                = 4.0L;
   static const double half                = 0.5L;
   static const double two_pi              = two*pi;
   static const double two_thirds          = two / three;
   static const double teensy              = 1.0e-100L;

   static const std::complex<double> iunit = std::complex<double>(zero, one );
   static const std::complex<double> czero = std::complex<double>(zero, zero);
   static const std::complex<double> cone  = std::complex<double>(one,  zero);
   static const std::complex<double> ctwo  = std::complex<double>(two,  zero);
   static const std::complex<double> chalf = std::complex<double>(half, zero);

#endif

}

#endif
