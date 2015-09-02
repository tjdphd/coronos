#include "cls_fft.hpp"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~ non-CUDA Fourier-related ~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef HAVE_CUDA_H

void fft::fftwInitialize( stack& run ) {

  int n1; 
  run.stack_data.fetch("n1",   &n1);
  int n2;
  run.stack_data.fetch("n2",   &n2);
  int nr_in; 
  run.stack_data.fetch("n1n2", &nr_in);
  int nc_out;

  nc_out    = n1 * (((int)(0.5*n2)) + 1);

  /* ~ Forward field/layer transforms: ~ */

  r_in      = (double *)               fftw_malloc(sizeof(double)               * nr_in  );
  cplx_out  = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * nc_out );

    p_lay_for = fftw_plan_dft_r2c_2d(n1, n2, r_in, reinterpret_cast<fftw_complex*>(cplx_out), FFTW_MEASURE);

  /* ~ Reverse field/layer transforms: ~ */

  cplx_in   = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * nc_out );
  r_out     = (double *)               fftw_malloc(sizeof(double)               * nr_in  );
  
  p_lay_rev = fftw_plan_dft_c2r_2d(n1, n2, reinterpret_cast<fftw_complex*>(cplx_in), r_out, FFTW_MEASURE);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwFinalize() {

  fftw_destroy_plan(p_lay_for);
  fftw_destroy_plan(p_lay_rev);

  fftw_free(r_in);
  fftw_free(r_out);

  fftw_free(cplx_in); 
  fftw_free(cplx_out); 

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardAll( stack& run ) {

  InputOutputArray& U = run.U;               /* ~ raw input array                         ~ */

  ComplexArray& U0 = run.U0;
  ComplexArray& U1 = run.U1;
  ComplexArray& U2 = run.U2;
  ComplexArray& U3 = run.U3;

  int n1; 
  run.stack_data.fetch("n1"   , &n1      );
  int n2; 
  run.stack_data.fetch("n2"   , &n2      );
  int n1n2c;
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2; 
  run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; 
  run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;
  run.stack_data.fetch("iu3"  , &n_flds  );

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx;

  for (int i_f = 0; i_f < n_flds; i_f++)  {
    for ( int i_l = 0; i_l < n_layers; i_l++) {

      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = zero; }
      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = U[k][i_l][i_f]; }
      for (unsigned k = 0; k < nc;   k++) { cplx_out[k] = (std::complex<double>) zero; }

      fftw_execute(p_lay_for);

      strt_idx = i_l * nc;

      switch(i_f) {
      case(0) :
        for (unsigned k = 0; k < nc;   k++) { U0[(strt_idx + k)] = cplx_out[k]; }
        break;
      case(1) :  
        for (unsigned k = 0; k < nc;   k++) { U1[(strt_idx + k)] = cplx_out[k]; }
        break;
      case(2) :  
        for (unsigned k = 0; k < nc;   k++) { U2[(strt_idx + k)] = cplx_out[k]; }
        break;
      case(3) :  
        for (unsigned k = 0; k < nc;   k++) { U3[(strt_idx + k)] = cplx_out[k]; }
        break;
      }
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseAll( stack& run ) {

  InputOutputArray& U = run.U;               /* ~ raw input array                         ~ */

  ComplexArray& U0 = run.U0;
  ComplexArray& U1 = run.U1;
  ComplexArray& U2 = run.U2;
  ComplexArray& U3 = run.U3;

  int n1; 
  run.stack_data.fetch("n1"   , &n1      );
  int n2; 
  run.stack_data.fetch("n2"   , &n2      );
  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2; 
  run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; 
  run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;
  run.stack_data.fetch("iu3"  , &n_flds  );

  ComplexArray::size_type nc = n1n2c;

  double scale = (double) one/((double) (n1n2));
  
  unsigned strt_idx;

  for (int i_f = 0; i_f < n_flds; i_f++)  {
    for ( int i_l = 0; i_l < n_layers; i_l++) {
      
      strt_idx = (i_l * nc);

      for (unsigned   k = 0; k < nc; k++) { cplx_in[k] = (std::complex<double>) zero; }
      switch(i_f) {
      case(0) :
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U0[strt_idx + k];  }
        break;
      case(1) :
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U1[strt_idx + k];  }
        break;
      case(2) :
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U2[strt_idx + k];  }
        break;
      case(3) :
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U3[strt_idx + k];  }
        break;
      }

      for (int k = 0; k < n1n2; k++) { r_out[k] = zero; }
      fftw_execute(p_lay_rev);

      for (int k = 0; k < n1n2; k++) { U[k][i_l][i_f] = (scale * r_out[k]); }

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardLayerofField ( stack& run, int i_l, int i_f ) {

  InputOutputArray& U = run.U;               /* ~ raw input array                         ~ */

  ComplexArray& U0 = run.U0;
  ComplexArray& U1 = run.U1;
  ComplexArray& U2 = run.U2;
  ComplexArray& U3 = run.U3;

  int n1; 
  run.stack_data.fetch("n1"   , &n1      );
  int n2; 
  run.stack_data.fetch("n2"   , &n2      );
  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2; 
  run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; 
  run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;
  run.stack_data.fetch("iu3"  , &n_flds  );

  ComplexArray::size_type nc = n1n2c;

  for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = zero; }
  for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = U[k][i_l][i_f]; }
  for (unsigned k = 0; k < nc;   k++) { cplx_out[k] = (std::complex<double>) zero; }

  fftw_execute(p_lay_for);

  unsigned strt_idx = i_l * nc;

  switch(i_f) {
  case(0) :
    for (unsigned k = 0; k < nc;   k++) { U0[(strt_idx + k)] = cplx_out[k]; }
    break;
  case(1) :  
    for (unsigned k = 0; k < nc;   k++) { U1[(strt_idx + k)] = cplx_out[k]; }
    break;
  case(2) :  
    for (unsigned k = 0; k < nc;   k++) { U2[(strt_idx + k)] = cplx_out[k]; }
    break;
  case(3) :  
    for (unsigned k = 0; k < nc;   k++) { U3[(strt_idx + k)] = cplx_out[k]; }
    break;
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseLayerofField ( stack& run, int i_l, int i_f) {

  InputOutputArray& U = run.U;               /* ~ raw input array                         ~ */

  ComplexArray& U0 = run.U0;
  ComplexArray& U1 = run.U1;
  ComplexArray& U2 = run.U2;
  ComplexArray& U3 = run.U3;

  int n1; 
  run.stack_data.fetch("n1"   , &n1      );
  int n2; 
  run.stack_data.fetch("n2"   , &n2      );
  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2; 
  run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; 
  run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;
  run.stack_data.fetch("iu3"  , &n_flds  );

  ComplexArray::size_type nc = n1n2c;

  double scale      = (double) one/((double) (n1n2));
  
  unsigned strt_idx = (i_l * nc);

  for (unsigned   k = 0; k < nc; k++) { cplx_in[k] = czero; }
  switch(i_f) {
  case(0) :
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U0[strt_idx + k];  }
    break;
  case(1) :
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U1[strt_idx + k];  }
    break;
  case(2) :
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U2[strt_idx + k];  }
    break;
  case(3) :
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U3[strt_idx + k];  }
    break;
  }

  for (int k = 0; k < n1n2; k++) { r_out[k] = zero; }
  fftw_execute(p_lay_rev);

  for (int k = 0; k < n1n2; k++) { U[k][i_l][i_f] = (scale * r_out[k]); }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardRaw( stack& run, RealArray& Rin, ComplexArray& Cout) {

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2; 
  run.stack_data.fetch("n1n2" , &n1n2    );
  int iu2;
  run.stack_data.fetch("iu2"  , &iu2     );

  unsigned c_strt_idx = 0;
  unsigned r_strt_idx = 0;

  for (unsigned i_l   = 0; i_l < iu2; i_l++) {

    c_strt_idx        = i_l * n1n2c;
    r_strt_idx        = i_l * n1n2;

    for (unsigned k   = 0 ; k < n1n2 ; k++) { r_in[k]              =  zero;               }
    for (unsigned k   = 0 ; k < n1n2c; k++) { cplx_out[k]          = czero;               }
    for (unsigned k   = 0 ; k < n1n2 ; k++) { r_in[k]              = Rin[r_strt_idx + k]; }
    fftw_execute(p_lay_for);
    for (unsigned k   = 0 ; k < n1n2c ; k++){ Cout[c_strt_idx + k] = cplx_out[k];         }

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseRaw( stack& run, ComplexArray& Cin, RealArray& Rout) {

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2; 
  run.stack_data.fetch("n1n2" , &n1n2    );
  int iu2;
  run.stack_data.fetch("iu2"  , &iu2     );

  assert(n1n2        != 0 );

  double scale        = (double) one/((double) (n1n2));
  
  unsigned c_strt_idx = 0;
  unsigned r_strt_idx = 0;

  for (unsigned i_l   = 0; i_l < iu2; i_l++) {

    c_strt_idx        = i_l * n1n2c;
    r_strt_idx        = i_l * n1n2;

  for (unsigned k     = 0 ; k < n1n2c; k++) { cplx_in[k]           = czero;               }
  for (unsigned k     = 0 ; k < n1n2 ; k++) { r_out[k]             =  zero;               }
  for (unsigned k     = 0 ; k < n1n2c; k++) { cplx_in[k]           = Cin[c_strt_idx + k]; }
  fftw_execute(p_lay_rev);
  for (unsigned k     = 0 ; k < n1n2 ; k++) { Rout[r_strt_idx + k] = (scale * r_out[k]);  }
  
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseIC(ComplexArray& Cin, RealArray& Rout ) {

 int rank;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2c;
  n1n2c           = Cin.size();
  int n1n2;
  n1n2            = Rout.size();

  assert(n1n2    != 0);

  double scale    = (double) one/((double) (n1n2));

  for (unsigned k = 0 ; k < n1n2c; k++) {cplx_in[k]  = czero; }
  for (unsigned k = 0 ; k < n1n2 ; k++) {r_out[k]    =  zero; }
  for (unsigned k = 0 ; k < n1n2c; k++) {cplx_in[k]  = Cin[k];}

  fftw_execute(p_lay_rev);
  for (unsigned k = 0 ; k < n1n2 ; k++) {Rout[k]     = (scale * r_out[k]);}

}
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardIC( RealArray& Rin, ComplexArray& Cout) {

  int n1n2c;
  n1n2c           = Cout.size();
  int n1n2;
  n1n2            = Rin.size();

  for (unsigned k = 0 ; k < n1n2 ; k++) {r_in[k]     =  zero;       }
  for (unsigned k = 0 ; k < n1n2c; k++) {cplx_out[k] = czero;       }
  for (unsigned k = 0 ; k < n1n2 ; k++) {r_in[k]     = Rin[k];      }
  fftw_execute(p_lay_for);
  for (unsigned k = 0 ; k < n1n2c; k++) {Cout[k]     = cplx_out[k]; }

}

#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// void fft::test_run( stack& run) {
// 
//   int i_t;
// 
// 
// }
// void fft::test_solve(lcsolve& solve ) {
// 
//   int i_t;
// 
// }
// void fft::test_both(stack& run, lcsolve& solve ) {
// 
