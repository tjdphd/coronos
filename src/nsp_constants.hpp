#ifndef NSP_CONSTANTS
#define NSP_CONSTANTS

#include<vector>
#include<complex>

namespace constants {

typedef double ***                         InputOutputArray;
typedef std::vector<double>                RealArray;
typedef std::vector<std::complex<double> > ComplexArray;

   static const double pi                  = 3.14159265358979323846264338327950288L;
   static const double zero                = 0.0L;
   static const double one                 = 1.0L;
   static const double two                 = 2.0L;
   static const double four                = 4.0L;
   static const double half                = 0.5L;
   static const double two_pi              = two*pi;
   static const double two_thirds          = 0.66666666666666666666666666666666666L;
   static const double teensy              = 1.0e-100L;

   static const std::complex<double> iunit = std::complex<double>(zero, one );
   static const std::complex<double> czero = std::complex<double>(zero, zero);
   static const std::complex<double> cone  = std::complex<double>(one,  zero);
   static const std::complex<double> ctwo  = std::complex<double>(two,  zero);
   static const std::complex<double> chalf = std::complex<double>(half, zero);

}

#endif
