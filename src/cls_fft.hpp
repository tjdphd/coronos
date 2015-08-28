#ifndef UTIL_FFT
#define UTIL_FFT

#include "nsp_constants.hpp"
#include "cls_stack.hpp"
#include "cls_lcsolve.hpp"

class fft {

/* ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft  ~ */

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

    void fftwForwardAll( stack& run, lcsolve& solve);      /* ~ Forward FFT all fields all layers    ~ */
    void fftwReverseAll( stack& run, lcsolve& solve);      /* ~ Reverst FFT all fields all layers    ~ */

    void fftwForwardLayerofField ( stack& run, lcsolve& solve, int layer, int field );
    void fftwReverseLayerofField ( stack& run, lcsolve& solve, int layer, int field );

    void fftwForwardRaw( stack& run, RealArray&    Rin, ComplexArray& Cout);
    void fftwReverseRaw( stack& run, ComplexArray& Cin, RealArray&    Rout);


/* ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft ~ fft  ~ */

};

#endif
