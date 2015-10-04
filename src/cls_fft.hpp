#ifndef UTIL_FFT
#define UTIL_FFT

#include "nsp_constants.hpp"
#include "cls_stack.hpp"
#include <assert.h>
//#include "cls_lcsolve.hpp"
//#include "cls_redhallmhd.hpp"

#ifdef HAVE_CUDA_H
#include "cls_lcsolve_cuda_ext.hpp"
#else
#include<fftw3.h>
#endif

class fft {

public:

    double          * r_in;                                /* ~ input and output data arguments      ~ */
    double          * r_out;                               /* ~ to FFT routines.                     ~ */

    std::complex<double> * cplx_in;
    std::complex<double> * cplx_out;

    fftw_plan      p_lay_for;                              /* ~ For establishing plans for forward   ~ */
    fftw_plan      p_lay_rev;                              /* ~ reverse FFT's of layers              ~ */

    void fftwInitialize( stack& run);                      /* ~ For allocating and deallocating "in" ~ */
    void fftwFinalize();                                   /* ~ and "out" arguments of FFT's, and    ~ */
                                                           /* ~ for initializing and "destroying"    ~ */
                                                           /* ~ FFT plans.                           ~ */

    void fftwKInit( stack& run );
    void fftwKFree( stack& run );

    RealArray rt;                                          /* ~ Fourier Transform Related               ~ */

    void fftwrtInit( stack& run );
    void fftwrtFree( );

    void fftwForwardAll( stack& run);                      /* ~ Forward FFT all fields all layers    ~ */
    void fftwReverseAll( stack& run);                      /* ~ Reverst FFT all fields all layers    ~ */

    void fftwForwardLayerofField ( stack& run, int layer, int field );
    void fftwReverseLayerofField ( stack& run, int layer, int field );

    void fftwForwardRaw( stack& run, RealArray&    Rin, ComplexArray& Cout);
    void fftwReverseRaw( stack& run, ComplexArray& Cin, RealArray&    Rout);

    void fftwForwardIC(RealArray&    Rin, ComplexArray& Cout);
    void fftwReverseIC(ComplexArray& Cin, RealArray&    Rout);

};

#endif
