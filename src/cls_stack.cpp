/* class stack (implementation)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * Longcope-type staggered stack of slabs
 *
 */

#include "cls_stack.hpp"

/* ~~~~~~~~~~~~~~~~ */
/* ~ Constructors ~ */
/* ~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

stack::stack() : canvas::canvas() {

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

stack::stack(std::string coronos_in) : canvas::canvas(coronos_in) {

  init_stack_data();
  allocUi();
  initxyz();
//initz();
//  kInit();
//  rtInit();

}

/* ~~~~~~~~~~~~~~~~ */
/* ~ initializers ~ */
/* ~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::init_stack_data() {                  /* ~ gather/infer information to be           ~ */
                                                 /* ~ included in stack_data container         ~ */

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* ~ incoming parameters from palette         ~ */
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  std::string model;                             /* ~ reduced mhd or hall-mrhd                 ~ */
  palette.fetch("model",&model);
  int p1;                                        /* ~ power of 2 specifying resolution in x    ~ */
  palette.fetch("p1"   ,&p1   );
  int p2;                                        /* ~ power of 2 specifying resolution in y    ~ */
  palette.fetch("p2"   ,&p2   );
  int p3;                                        /* ~ total number of layers in z              ~ */
  palette.fetch("p3"   ,&p3   );
  int np;                                        /* ~ np number of processes                   ~ */
  palette.fetch("np"   ,&np   );
  double zl;                                     /* ~ length in z of computational domain      ~ */
  palette.fetch("zl"   ,&zl   );

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* ~ to be made parameters of stack_data      ~ */
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  std::string resolution;                        /* ~ full 3-d resolution string                ~ */

  int n1    = (int) pow(2.0, p1);                /* ~ number of x-coordinates in a layer        ~ */
  int n2    = (int) pow(2.0, p2);                /* ~ number of y-coordinates in a layer        ~ */
  int n3    =                p3 ;                /* ~ number of interior layers per process     ~ */

  int n1n2  = n1*n2;                             /* ~ total number points on a (real) layer     ~ */
  int n1n2c = n1 * (((int)(half * n2)) + 1);     /* ~ total number points on a (Fourier) layer  ~ */

  int iu1   = n1n2;                              /* ~ dimension of first index of U             ~ */
  int iu2   = n3+2;                              /* ~ dimension of second index of U            ~ */

                                                 /* ~ Ui's should have dimensions:              ~ */
                                                 /* ~ Ui[n1n2c * iu2]                           ~ */
                                                 /* ~                                           ~ */
                                                 /* ~ NOTE: relative to old codes:              ~ */
                                                 /* ~ pbot(:) = U0[0:(n1n2c - 1)]               ~ */
                                                 /* ~ atop(:) = U1[n1n2c*(iu2-1):(n1ncc*iu2)-1] ~ */

  int iu3;                                       /* ~ dimension of third index of U             ~ */
                                                 /* ~ counts number of fields in plasma model   ~ */
  if (model.compare("rmhd") == 0) iu3 = 2;       /* ~ fix number of field variables             ~ */
  if (model.compare("hall") == 0) iu3 = 4;

  double dz        = zl/(n3*np);                 /* ~ layer separation in z                     ~ */

  int izres        = (int) (n3 * np)/zl;         /* ~ integer effective resolution in z         ~ */

  std::string xres = static_cast<std::ostringstream*>( &(std::ostringstream() << n1   ) ) -> str();
  std::string yres = static_cast<std::ostringstream*>( &(std::ostringstream() << n2   ) ) -> str();
  std::string zres = static_cast<std::ostringstream*>( &(std::ostringstream() << izres) ) -> str();

  if (xres.compare(yres) == 0 ) resolution.assign(xres + "_" + zres);
  else             resolution.assign(xres + "_" + yres + "_" + zres);

  std::string pname;                             /* ~ for containing parameter names            ~ */
  std::string padjust;                           /* ~ for specifying adjustability              ~ */
  
  padjust.assign("rfx"  );                       /* ~ assigning parameters that are run-fixed   ~ */

  pname.assign("n1"     );
  stack_data.emplace(pname, n1,         padjust);
  pname.assign("n2"     );
  stack_data.emplace(pname, n2,         padjust);
  pname.assign("n3"     );
  stack_data.emplace(pname, n3,         padjust);
  pname.assign("n1n2"   );
  stack_data.emplace(pname, n1n2,       padjust);
  pname.assign("n1n2c"  );
  stack_data.emplace(pname, n1n2c,      padjust);
  pname.assign("iu1"    );
  stack_data.emplace(pname, iu1,        padjust);
  pname.assign("iu2"    );
  stack_data.emplace(pname, iu2,        padjust);
  pname.assign("iu3"    );
  stack_data.emplace(pname, iu3,        padjust);
  pname.assign("dz"     );
  stack_data.emplace(pname, dz,         padjust);
  pname.assign("res_str");
  stack_data.emplace(pname, resolution, padjust);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~ Allocators / De-allocators ~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::allocUi() {                          /* ~ U is the input/output array for the fields ~ */
                                                 /* ~ NOTE: U holds the real-space fields        ~ */
  int iu1, iu2, iu3;                             /* ~       and is read and written to dataframe ~ */
                                                 /* ~       files for each process               ~ */

  stack_data.fetch("iu1", &iu1);
  stack_data.fetch("iu2", &iu2);
  stack_data.fetch("iu3", &iu3);

  std::cout << "iu1 = " << iu1 << std::endl;
  std::cout << "iu2 = " << iu2 << std::endl;
  std::cout << "iu3 = " << iu3 << std::endl;


  U = new double**[iu1];                         /* ~ allocate U dynamically using dimensions    ~ */

  for (int i = 0; i< iu1; ++i) {                 /* ~ determined in init_stack_data              ~ */
    U[i] = new double*[iu2];
    for (int j = 0; j < iu2; ++j) {
      U[i][j] = new double[iu3];

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::zeroU() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  int iu1, iu2, iu3;

  stack_data.fetch("iu1", &iu1);
  stack_data.fetch("iu2", &iu2);
  stack_data.fetch("iu3", &iu3);

  int i, j, k;

  for(k = 0; k < iu3; ++k) {
    for(j = 0; j < iu2; ++j) {
      for(i = 0; i < iu1; ++i) {

          U[i][j][k] = zero;

      }
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::deallocUi() {

  int iu1, iu2;

  stack_data.fetch("iu1", &iu1);
  stack_data.fetch("iu2", &iu2);

  for (int i = 0; i< iu1; ++i) {
    for (int j = 0; j < iu2; ++j) delete [] U[i][j];
    delete [] U[i];
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::rtInit( ) {                          /* ~ initialize de-aliasing array              ~ */

  int n1;                                        /* ~ number of x-coordinates                   ~ */
  stack_data.fetch("n1",    &n1);
  int n2;
  stack_data.fetch("n2",    &n2);                /* ~ number of y-coordinates                   ~ */
  int n1n2c;
  stack_data.fetch("n1n2c", &n1n2c);             /* ~ number of Fourier space points per layer  ~ */

  RealArray::size_type nc;                       /* ~ a vector size-type version of n1n2c       ~ */
  nc  = n1n2c;

  int n1h = (((int)(half*n1)) + 1);
  int n2h = (((int)(half*n2)) + 1);

  double threshold = two_thirds * sqrt(two * ((( half * n1 - one) * two_pi) * ((half * n1 - one) * two_pi)));

  rt.reserve(nc);                                /* ~ create space in rt                        ~ */

  RealArray::size_type ndx;

  for (unsigned i = 0; i < n1; i++) {            /* ~ initialization                            ~ */
    for (unsigned j = 0; j < n1h; j++) {

      ndx = i * n2h + j;

      if ( sqrt(std::abs(k2[ndx])) < threshold ) rt[ndx] = one;
      else                                       rt[ndx] = zero;

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::rtFree() {   /* ~ this is a bit excessive no? maybe just put in destructor or something ? ~ */

 rt.resize(0);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::kInit() {                            /* ~ initialize the wave number arrays kx, ky, k2, and inv_k2 ~ */

  int n1;                                        /* ~ number of coordinates in x                               ~ */
  stack_data.fetch("n1",    &n1);
  int n2; 
  stack_data.fetch("n2",    &n2);                /* ~ number of coordinates in y                               ~ */
  int nc;
  stack_data.fetch("n1n2c", &nc);                /* ~ number of Fourier space points per layer                 ~ */

  int n1h    = (((int)(half*n1)) + 1);
  int n2h    = (((int)(half*n2)) + 1);

  RealArray kp;                                  /* ~ temporary definition of convenience                      ~ */

  kp.reserve(n1);

      kx.reserve(nc);                            /* ~ these are all lcstack members and carry the wave-number  ~ */
      ky.reserve(nc);                            /* ~ information                                              ~ */
      k2.reserve(nc);
  inv_k2.reserve(nc);

                                                 /* ~ the initialization as done in gpu port of old code       ~ */
  int ndx;

  for (int i = 0;   i < n1h; i++ ) kp[i] = ((double) i) * two_pi;
  for (int i = n1h; i < n1;  i++ ) kp[i] = -kp[n1 - i];

  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n1h; j++) {
   
      ndx     = (i * n2h) + j;
      kx[ndx] = kp[j];
      ky[ndx] = kp[i];
    }
  }

  kp.resize(0);

  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n1h; j++) {

      ndx     = (i * n2h) + j;
      k2[ndx] = (kx[ndx] * kx[ndx]) +  (ky[ndx] * ky[ndx]);

      if (std::abs(k2[ndx]) > teensy) inv_k2[ndx] = one/k2[ndx];
      else inv_k2[ndx] = zero;

    }
  }

/* ~ save for debugging ~ */

/*   int rank;
 *  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 * 
 *  if (rank == 0) {
 * 
 *    std::cout << "ky: " << std::endl << std::endl;
 * 
 *    for (int i = 0; i < n1; ++i) {
 *      for (int j = 0; j < n1/2 + 1; ++j) {
 * 
 *        ndx = (i * n2h) + j;
 *        std::cout << std::setw(10) << std::right << std::setprecision(4) << std::scientific << ky[ndx] << "  ";
 * 
 *      }
 *        std::cout << std::endl;
 *    }
 *  }
 */  

/* ~ save for debugging ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::kFree() {                            /* ~ also possibly a bit excessive ~ */

    kx.resize(0);
    ky.resize(0);
    k2.resize(0);
    inv_k2.resize(0);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~ non-CUDA Fourier-related ~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//#ifndef HAVE_CUDA_H

//void stack::fftwInitialize( ) {
//
//  int n1; 
//  stack_data.fetch("n1",   &n1);
//  int n2;
//  stack_data.fetch("n2",   &n2);
//  int nr_in; 
//  stack_data.fetch("n1n2", &nr_in);
//  int nc_out;
//
//  nc_out    = n1 * (((int)(0.5*n2)) + 1);
//
//  /* ~ Forward field/layer transforms: ~ */
//
//  r_in      = (double *)               fftw_malloc(sizeof(double)               * nr_in  );
//  cplx_out  = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * nc_out );
//
//  p_lay_for = fftw_plan_dft_r2c_1d(nr_in, r_in, reinterpret_cast<fftw_complex*>(cplx_out), FFTW_MEASURE);
//
//  /* ~ Reverse field/layer transforms: ~ */
//
//  cplx_in   = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * nc_out );
//  r_out     = (double *)          fftw_malloc(sizeof(double)          * nr_in  );
//  
//  p_lay_rev = fftw_plan_dft_c2r_1d(nr_in, reinterpret_cast<fftw_complex*>(cplx_in), r_out, FFTW_MEASURE);
//
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::fftwFinalize() {
//
//  fftw_destroy_plan(p_lay_for);
//  fftw_destroy_plan(p_lay_rev);
//
//  fftw_free(r_in);
//  fftw_free(r_out);
//
//  fftw_free(cplx_in); 
//  fftw_free(cplx_out); 
//
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::fftwForwardAll( lcsolve& solve) {
//
//  ComplexArray& U0 = solve.U0;
//  ComplexArray& U1 = solve.U1;
//  ComplexArray& U2 = solve.U2;
//  ComplexArray& U3 = solve.U3;
//
//  int n1; 
//  stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  stack_data.fetch("n2"   , &n2      );
//  int n1n2c;
//  stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  stack_data.fetch("iu3"  , &n_flds  );
//
//  ComplexArray::size_type nc = n1n2c;
//
//  unsigned strt_idx;
//
//  for (int i_f = 0; i_f < n_flds; i_f++)  {
//    for ( int i_l = 0; i_l < n_layers; i_l++) {
//
//      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = zero; }
//      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = U[k][i_l][i_f]; }
//      for (unsigned k = 0; k < nc;   k++) { cplx_out[k] = (std::complex<double>) zero; }
//
//      fftw_execute(p_lay_for);
//
//      strt_idx = i_l * nc;
//
//      switch(i_f) {
//      case(0) :
//        for (unsigned k = 0; k < nc;   k++) { U0[(strt_idx + k)] = cplx_out[k]; }
//        break;
//      case(1) :  
//        for (unsigned k = 0; k < nc;   k++) { U1[(strt_idx + k)] = cplx_out[k]; }
//        break;
//      case(2) :  
//        for (unsigned k = 0; k < nc;   k++) { U2[(strt_idx + k)] = cplx_out[k]; }
//        break;
//      case(3) :  
//        for (unsigned k = 0; k < nc;   k++) { U3[(strt_idx + k)] = cplx_out[k]; }
//        break;
//      }
//    }
//  }
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::fftwReverseAll( lcsolve& solve ) {
//
//  ComplexArray& U0 = solve.U0;
//  ComplexArray& U1 = solve.U1;
//  ComplexArray& U2 = solve.U2;
//  ComplexArray& U3 = solve.U3;
//
//  int n1; 
//  stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  stack_data.fetch("n2"   , &n2      );
//  int n1n2c; 
//  stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  stack_data.fetch("iu3"  , &n_flds  );
//
//  ComplexArray::size_type nc = n1n2c;
//
//  double scale = (double) one/((double) (n1n2));
//  
//  unsigned strt_idx;
//
//  for (int i_f = 0; i_f < n_flds; i_f++)  {
//    for ( int i_l = 0; i_l < n_layers; i_l++) {
//      
//      strt_idx = (i_l * nc);
//
//      for (unsigned   k = 0; k < nc; k++) { cplx_in[k] = (std::complex<double>) zero; }
//      switch(i_f) {
//      case(0) :
//        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U0[strt_idx + k];  }
//        break;
//      case(1) :
//        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U1[strt_idx + k];  }
//        break;
//      case(2) :
//        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U2[strt_idx + k];  }
//        break;
//      case(3) :
//        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U3[strt_idx + k];  }
//        break;
//      }
//
//      for (int k = 0; k < n1n2; k++) { r_out[k] = zero; }
//      fftw_execute(p_lay_rev);
//
//      for (int k = 0; k < n1n2; k++) { U[k][i_l][i_f] = (scale * r_out[k]); }
//
//    }
//  }
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::fftwForwardLayerofField ( lcsolve& solve, int i_l, int i_f ) {
//
//  ComplexArray& U0 = solve.U0;
//  ComplexArray& U1 = solve.U1;
//  ComplexArray& U2 = solve.U2;
//  ComplexArray& U3 = solve.U3;
//
//  int n1; 
//  stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  stack_data.fetch("n2"   , &n2      );
//  int n1n2c; 
//  stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  stack_data.fetch("iu3"  , &n_flds  );
//
//  ComplexArray::size_type nc = n1n2c;
//
//  for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = zero; }
//  for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = U[k][i_l][i_f]; }
//  for (unsigned k = 0; k < nc;   k++) { cplx_out[k] = (std::complex<double>) zero; }
//
//  fftw_execute(p_lay_for);
//
//  unsigned strt_idx = i_l * nc;
//
//  switch(i_f) {
//  case(0) :
//    for (unsigned k = 0; k < nc;   k++) { U0[(strt_idx + k)] = cplx_out[k]; }
//    break;
//  case(1) :  
//    for (unsigned k = 0; k < nc;   k++) { U1[(strt_idx + k)] = cplx_out[k]; }
//    break;
//  case(2) :  
//    for (unsigned k = 0; k < nc;   k++) { U2[(strt_idx + k)] = cplx_out[k]; }
//    break;
//  case(3) :  
//    for (unsigned k = 0; k < nc;   k++) { U3[(strt_idx + k)] = cplx_out[k]; }
//    break;
//  }
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::fftwReverseLayerofField ( lcsolve& solve, int i_l, int i_f) {
//
//  ComplexArray& U0 = solve.U0;
//  ComplexArray& U1 = solve.U1;
//  ComplexArray& U2 = solve.U2;
//  ComplexArray& U3 = solve.U3;
//
//  int n1; 
//  stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  stack_data.fetch("n2"   , &n2      );
//  int n1n2c; 
//  stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  stack_data.fetch("iu3"  , &n_flds  );
//
//  ComplexArray::size_type nc = n1n2c;
//
//  double scale      = (double) one/((double) (n1n2));
//  
//  unsigned strt_idx = (i_l * nc);
//
//  for (unsigned   k = 0; k < nc; k++) { cplx_in[k] = czero; }
//  switch(i_f) {
//  case(0) :
//    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U0[strt_idx + k];  }
//    break;
//  case(1) :
//    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U1[strt_idx + k];  }
//    break;
//  case(2) :
//    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U2[strt_idx + k];  }
//    break;
//  case(3) :
//    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U3[strt_idx + k];  }
//    break;
//  }
//
//  for (int k = 0; k < n1n2; k++) { r_out[k] = zero; }
//  fftw_execute(p_lay_rev);
//
//  for (int k = 0; k < n1n2; k++) { U[k][i_l][i_f] = (scale * r_out[k]); }
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::fftwForwardRaw( RealArray& Rin, ComplexArray& Cout) {
//
//  int n1n2c; 
//  stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  stack_data.fetch("n1n2" , &n1n2    );
//  int iu2;
//  stack_data.fetch("iu2"  , &iu2     );
//
//  unsigned c_strt_idx = 0;
//  unsigned r_strt_idx = 0;
//
//  for (unsigned i_l   = 0; i_l < iu2; i_l++) {
//
//    c_strt_idx        = i_l * n1n2c;
//    r_strt_idx        = i_l * n1n2;
//
//    for (unsigned k   = 0 ; k < n1n2 ; k++) { r_in[k]              =  zero;               }
//    for (unsigned k   = 0 ; k < n1n2c; k++) { cplx_out[k]          = czero;               }
//    for (unsigned k   = 0 ; k < n1n2 ; k++) { r_in[k]              = Rin[r_strt_idx + k]; }
//    fftw_execute(p_lay_for);
//    for (unsigned k   = 0 ; k < n1n2c ; k++){ Cout[c_strt_idx + k] = cplx_out[k];         }
//
//  }
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::fftwReverseRaw( ComplexArray& Cin, RealArray& Rout) {
//
//  int n1n2c; 
//  stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  stack_data.fetch("n1n2" , &n1n2    );
//  int iu2;
//  stack_data.fetch("iu2"  , &iu2     );
//
//  assert(n1n2        != 0 );
//
//  double scale        = (double) one/((double) (n1n2));
//  
//  unsigned c_strt_idx = 0;
//  unsigned r_strt_idx = 0;
//
//  for (unsigned i_l   = 0; i_l < iu2; i_l++) {
//
//    c_strt_idx        = i_l * n1n2c;
//    r_strt_idx        = i_l * n1n2;
//
//  for (unsigned k     = 0 ; k < n1n2c; k++) { cplx_in[k]           = czero;               }
//  for (unsigned k     = 0 ; k < n1n2 ; k++) { r_out[k]             =  zero;               }
//  for (unsigned k     = 0 ; k < n1n2c; k++) { cplx_in[k]           = Cin[c_strt_idx + k]; }
//  fftw_execute(p_lay_rev);
//  for (unsigned k     = 0 ; k < n1n2 ; k++) { Rout[r_strt_idx + k] = (scale * r_out[k]);  }
//  
//  }
//}
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//#endif
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* ~ Initialization            ~ */
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

std::string stack::getLastDataFilename(int srun) {

  int rank;
  run_data.fetch("rank",      &rank   );

  std::string prefix;
  palette.fetch("prefix",     &prefix );
  std::string res_str;
  stack_data.fetch("res_str", &res_str);

  std::string data_file = "not_this_one";

  std::string rnk_str;
  std::string srn_str;

  rnk_str               = static_cast<std::ostringstream*>( &(std::ostringstream() << rank ) ) -> str();
  srn_str               = static_cast<std::ostringstream*>( &(std::ostringstream() << srun ) ) -> str();

  int rnk_len           = rnk_str.length();

  switch(rnk_len) {

  case(1) : rnk_str     = "00" + rnk_str;
            break;
  case(2) : rnk_str     =  "0" + rnk_str;
            break;
  default : std::cout << "this can't be right " << std::endl;

  }

  data_file             = prefix + "_" + res_str + "." + rnk_str + ".ots" + srn_str;

  return data_file;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  std::string stack::getNextDataFilename() {

  std::string data_file = "not_this_one";

  std::string prefix;
  palette.fetch(   "prefix" , &prefix );

  int rank; 
  run_data.fetch(  "rank"   , &rank   );
  std::string rnk_str;
  rnk_str     = static_cast<std::ostringstream*>( &(std::ostringstream() << rank ) ) -> str();

  int srun;
  palette.fetch(   "srun"   , &srun   );
  std::string srn_str;
  srn_str     = static_cast<std::ostringstream*>( &(std::ostringstream() << srun ) ) -> str();

  std::string res_str;
  stack_data.fetch("res_str", &res_str);

  int rnk_len = rnk_str.length();

  switch(rnk_len) {

  case(1) : rnk_str = "00" + rnk_str;
            break;
  case(2) : rnk_str =  "0" + rnk_str;
            break;
  default : std::cout << "this can't be right " << std::endl;

  }

  data_file = prefix + "_" + res_str + "." + rnk_str + ".ots" + srn_str;

  return data_file;

  }

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::initxyz() {                     /* ~ Calculate x- and y-coordinates of layers ~ */

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1; 
  stack_data.fetch("n1",    &n1 );
  x.reserve(n1);

  double dx = one / ((double) n1);


  double next_x;
  for (int i = 0; i < n1; ++i) {

    next_x = ( (double) i) * dx; 
    x.push_back(next_x);

//  x[i] = next_x;                           /* ~ seems to work but I'll worry about this later ~ */

  }

  int n2;
  stack_data.fetch("n2",    &n2 );
  y.reserve(n2);

  double dy = one / ((double) n2);

  double next_y;

  for (int j = 0; j < n2; ++j) {

    next_y = ( (double) j) * dy;
    y.push_back(next_y);

//  y[j] = next_y;                           /* ~ seems to work but I'll worry about this later ~ */
//

  }

  int np;
  palette.fetch(   "np"  , &np );

  int    iu2;
  stack_data.fetch("iu2" , &iu2 );
  double dz;
  stack_data.fetch("dz"  , &dz );
  double zl;
  palette.fetch(   "zl"  , &zl );

  z.reserve(iu2);

  double next_z;

  for (int i = 0; i < iu2 - 1; ++i) {

    next_z =  ((double) (rank)) * (zl / ((double) (np))) + (((double) i) * dz);
    z.push_back(next_z);
 
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::initz() {

  int rank; 
  run_data.fetch(  "rank", &rank);
  int np;
  palette.fetch(   "np"  , &np );
  int iu2;
  stack_data.fetch("iu2" , &iu2 );
  double zl;
  palette.fetch(   "zl"  , &zl );
  double dz;
  stack_data.fetch("dz"  , &dz );

  z.reserve(iu2);

  double next_z;

  for (int i = 0; i < iu2 - 1; ++i) {

    next_z =  ((double) (rank)) * (zl / ((double) (np))) + (((double) i) * dz) ;
    z.push_back(next_z);
 
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~ */
/* ~ Destructor  ~ */
/* ~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

stack::~stack() {

  x.resize(0);
  y.resize(0);
  z.resize(0);

  rtFree();
  kFree();
//  deallocUi();

}
