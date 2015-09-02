/* class lcsolve (implementation)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * solver class for longcope solver designed to 
 * work on canvas class type "stack"
 *
 */

#include "cls_lcsolve.hpp"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~ */
/* ~ Constructors ~ */
/* ~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

lcsolve::lcsolve( stack& run ) {

  createFields( run );

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::Loop( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2;
  run.stack_data.fetch("n1n2", &n1n2);

  int l, ndt,    iptest;
  double tstart, t_cur, dt;

  run.palette.fetch("ndt",    &ndt   );
  run.palette.fetch("tstart", &t_cur );

  run.palette.fetch("iptest", &iptest);

  for (l = 0; l < ndt;l++) {

    /* ~ iptest conditional goes here              ~ */
    /* ~ mv, mb, etc.... initialization goes here  ~ */

    passAdjacentLayers ( "predict", run        );
//    updatePAJ(                 "predict", run, solve );   /* ~ P, A, and J contain un-updated/corrector-updated values ~ */
//    applyBC(                              run, solve );
//    setS(                      "predict", run, solve );   /* ~ set predictor S's                                       ~ */
//    setB(                      "predict", run, solve );   /* ~ set predictor Brackets                                  ~ */
//    setD(                      "predict", run, solve );   /* ~ set predictor finite differences                        ~ */
//    setAi(                                run, solve );   /* ~ set predictor A's                                       ~ */
//    solve.Step(                "predict", run        );   /* ~ execute predictor update                                ~ */
//
//    solve.passAdjacentLayers ( "correct", run        );
//    updatePAJ(                 "correct", run, solve );   /* ~ P, A, and J now contain predictor-updated values        ~ */
//    setS(                      "correct", run, solve );   /* ~ set corrector S's                                       ~ */
//    setB(                      "correct", run, solve );   /* ~ set corrector Brackets                                  ~ */
//    setD(                      "correct", run, solve );   /* ~ set corrector finite differences                        ~ */
//    setAi(                                run, solve );   /* ~ set corrector A's                                       ~ */
//    solve.Step(                "correct", run        );   /* ~ execute corrector update                                ~ */
//
//    run.palette.fetch("dt", &dt);
//    t_cur       = t_cur + dt;
//
//    updateTimeInc(                        run        );
//
///* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
//
////  int rank;
////  int mp,  n1n2c, n_layers;
////  unsigned strt_idx, stop_idx;
////
////  MPI_Comm_rank(MPI_COMM_WORLD,   &rank);
////
////  run.stack_data.fetch("n1n2c", &n1n2c   );
////  run.stack_data.fetch("n3"   , &n_layers);
////  run.palette.fetch("mp"    ,   &mp      );
////
////  ComplexArray& O            = solve.tU1;
////  ComplexArray::size_type nc = n1n2c;
////
////  if (rank == 1) {
////
////    strt_idx = nc;
////
////  }
////  else if (rank == 0) {
////
////    strt_idx = (n_layers + 1) * nc;
////
////  }
////  stop_idx   = strt_idx + n1n2c;
////
////  if (rank == 0) {
////
////    for (unsigned k = strt_idx; k < stop_idx; k++) {
////
////      std::cout << std::setw(24) << std::right << std::setprecision(16) << std::scientific << O[k].real() << std::endl;
////
////    }
////  }
//
///* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
//
//    /* ~ bookkeeping goes here                     ~ */
//    /* ~ time-step determination goes here         ~ */
//
//       
  }
//
//  updatePAJ("predict", run, solve );   /* ~ P, A, and J contain final corrector-updated values ~ */
//
//  run.palette.reset(   "tstart", t_cur    );
//
//  std::string prefix;
//  run.palette.fetch(   "prefix", &prefix  );
//
//  std::string res_str;
//  run.stack_data.fetch("res_str", &res_str);
//
//  int srun;
//  run.palette.fetch(   "srun"   , &srun   );
//
//  if (rank == 0) { run.palette.report(prefix + '_' +  res_str, srun); }
//

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::passAdjacentLayers( std::string str_step, stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Status * status = 0;

  ComplexArray& O     = run.U0;                    /* ~ for predictor case                          ~ */
  ComplexArray& H     = run.U1;                    /* ~ un-updated values are transferred           ~ */
  ComplexArray& Z     = run.U2;
  ComplexArray& V     = run.U3;

  ComplexArray& tO    = tU0;                      /* ~ for corrector case                          ~ */
  ComplexArray& tH    = tU1;                      /* ~ results from predictor step are transferred ~ */
  ComplexArray& tZ    = tU2;
  ComplexArray& tV    = tU3;

  int np;
  run.palette.fetch(   "np"   , &np      );       /* ~ number of processes                         ~ */

  int n1n2c;
  run.stack_data.fetch("n1n2c", &n1n2c   );       /* ~ number of complex elements in a layer       ~ */

  int n_layers;
  run.stack_data.fetch("n3"   , &n_layers);       /* ~ number of layers in stack                   ~ */

  unsigned n3_idx     =   n_layers       * n1n2c; /* ~ starting index for n3'th layer              ~ */
  unsigned atop_idx   = ( n_layers + 1 ) * n1n2c; /* ~ starting index for top boundary layer       ~ */
  unsigned abot_idx   =                    n1n2c; /* ~ starting index for first layer              ~ */


  std::string model;
  run.palette.fetch("model", &model);

  if (str_step.compare("predict") == 0) {         /* ~ predictor case                              ~ */

    if (rank != 0 ) {
      if (rank != np - 1) {

          MPI_Send(  &O[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank + 1     , MPI_COMM_WORLD        ); // send " p(:,n3)"
          if (model.compare("hall") == 0 ) {
            MPI_Send(&Z[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank         , MPI_COMM_WORLD        ); // send "bz(:,n3)"
          }

      }
          MPI_Recv(  &O.front(),   n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, rank         , MPI_COMM_WORLD, status); // receive " pbot"
          if (model.compare("hall") == 0 ) {
            MPI_Recv(&Z.front(),   n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, rank - 1     , MPI_COMM_WORLD, status); // receive "bzbot"
          }

          MPI_Send(  &H[abot_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, np + rank - 1, MPI_COMM_WORLD        ); // send " a(:,1)"
          if (model.compare("hall") == 0 ) {
            MPI_Send(&V[abot_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, np + rank    , MPI_COMM_WORLD        ); // send "vz(:,1)"
          }

        if ( rank != np - 1) {

          MPI_Recv(  &H[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank    , MPI_COMM_WORLD, status); // receive " atop"
          if (model.compare("hall") == 0 ) {
            MPI_Recv(&V[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank + 1, MPI_COMM_WORLD, status); // receive "vztop"
          }

        }
    }
    else {

          MPI_Send(  &O[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank + 1     , MPI_COMM_WORLD        ); // send " p(:,n3)"
          if (model.compare("hall") == 0 ) {
            MPI_Send(&Z[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank         , MPI_COMM_WORLD        ); // send "bz(:,n3)"
          }

          MPI_Recv(  &H[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank    , MPI_COMM_WORLD, status); // receive " atop"
          if (model.compare("hall") == 0 ) {
            MPI_Recv(&V[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank + 1, MPI_COMM_WORLD, status); // receive "vztop"
          }
    }
}
  else if(str_step.compare("correct") == 0) {     /* ~ corrector case                              ~ */

    if (rank != 0 ) {
      if (rank != np - 1) {

          MPI_Send(  &tO[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank + 1     , MPI_COMM_WORLD        ); // send " tp(:,n3)"
          if (model.compare("hall") == 0 ) {
            MPI_Send(&tZ[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank         , MPI_COMM_WORLD        ); // send "tbz(:,n3)"
          }

      }
          MPI_Recv(  &tO.front(),   n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, rank         , MPI_COMM_WORLD, status); // receive " pbot"
          if (model.compare("hall") == 0 ) {
            MPI_Recv(&tZ.front(),   n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, rank - 1     , MPI_COMM_WORLD, status); // receive "bzbot"
          }

          MPI_Send(  &tH[abot_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, np + rank - 1, MPI_COMM_WORLD        ); // send " ta(:,1)"
          if (model.compare("hall") == 0 ) {
            MPI_Send(&tV[abot_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank - 1, np + rank    , MPI_COMM_WORLD        ); // send "tvz(:,1)"
          }

        if ( rank != np - 1) {

          MPI_Recv(  &tH[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank    , MPI_COMM_WORLD, status); // receive " atop"
          if (model.compare("hall") == 0 ) {
            MPI_Recv(&tV[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank + 1, MPI_COMM_WORLD, status); // receive "vztop"
          }

        }
    }
    else {

          MPI_Send(  &tO[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank + 1     , MPI_COMM_WORLD        ); // send " tp(:,n3)"
          if (model.compare("hall") == 0 ) {
            MPI_Send(&tZ[n3_idx],   n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, rank         , MPI_COMM_WORLD        ); // send "tbz(:,n3)"
          }

          MPI_Recv(  &tH[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank    , MPI_COMM_WORLD, status); // receive " atop"
          if (model.compare("hall") == 0 ) {
            MPI_Recv(&tV[atop_idx], n1n2c, MPI::DOUBLE_COMPLEX, rank + 1, np + rank + 1, MPI_COMM_WORLD, status); // receive "vztop"
          }
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~ non-CUDA vector representations for fields ~ */

void lcsolve::createFields( stack& run ) {

    int n1n2c;
    run.stack_data.fetch("n1n2c", &n1n2c );        /* ~ number of complex elements per layer ~ */
  
    int iu2;
    run.stack_data.fetch("iu2"  , &iu2   );        /* ~ number of layers in stack            ~ */
  
    std::string model;
    run.palette.fetch("model"   , &model );
  
    tU0.reserve(n1n2c * iu2);                      /* ~ for predictor-step results           ~ */
    tU1.reserve(n1n2c * iu2);                      /* ~ Note: U0, U1, U2, & U3 are defined   ~ */
                                                   /* ~       on the stack.                  ~ */
  
    if (model.compare("hall") == 0 ) {
      tU2.reserve(n1n2c * iu2);
      tU3.reserve(n1n2c * iu2);
    }
  
    SE0.reserve(n1n2c);                            /* ~ the S -arrays. see documentation     ~ */
    SE1.reserve(n1n2c);
  
    if (model.compare("hall") == 0 ) {
      SE2.reserve(n1n2c);
      SE3.reserve(n1n2c);
    }
  
    SI0.reserve(n1n2c);
    SI1.reserve(n1n2c);
  
    if (model.compare("hall") == 0 ) {
      SI2.reserve(n1n2c);
      SI3.reserve(n1n2c);
    }
  
    B0.reserve(n1n2c * iu2);                      /* bracket terms                           ~ */
    B1.reserve(n1n2c * iu2);
  
    if (model.compare("hall") == 0 ) {
      B2.reserve(n1n2c * iu2);
      B3.reserve(n1n2c * iu2);
    }
  
    D0.reserve(n1n2c * iu2);                      /* finite differences in z                 ~ */
    D1.reserve(n1n2c * iu2);
  
    if (model.compare("hall") == 0 ) {
      D2.reserve(n1n2c * iu2);
      D3.reserve(n1n2c * iu2);
    }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::destroyFields() {

  tU0.resize(0);
  tU1.resize(0);
  tU2.resize(0);
  tU3.resize(0);

  SE0.resize(0);
  SE1.resize(0);
  SE2.resize(0);
  SE3.resize(0);

  SI0.resize(0);
  SI1.resize(0);
  SI2.resize(0);
  SI3.resize(0);

  D0.resize(0);
  D1.resize(0);
  D2.resize(0);
  D3.resize(0);

}

///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
///* ~ non-CUDA Fourier-related ~ */
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//#ifndef HAVE_CUDA_H
//
//void lcsolve::fftwInitialize( stack& run ) {
//
//  int n1; 
//  run.stack_data.fetch("n1",   &n1);
//  int n2;
//  run.stack_data.fetch("n2",   &n2);
//  int nr_in; 
//  run.stack_data.fetch("n1n2", &nr_in);
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
//void lcsolve::fftwFinalize() {
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
//void lcsolve::fftwForwardAll( stack& run ) {
//
//  InputOutputArray& U = run.U;
//
//  int n1; 
//  run.stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  run.stack_data.fetch("n2"   , &n2      );
//  int n1n2c;
//  run.stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  run.stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  run.stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  run.stack_data.fetch("iu3"  , &n_flds  );
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
//void lcsolve::fftwReverseAll( stack& run ) {
//
//  InputOutputArray& U = run.U;
//
//  int n1; 
//  run.stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  run.stack_data.fetch("n2"   , &n2      );
//  int n1n2c; 
//  run.stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  run.stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  run.stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  run.stack_data.fetch("iu3"  , &n_flds  );
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
//void lcsolve::fftwForwardLayerofField ( int i_l, int i_f, stack& run) {
//
//  InputOutputArray& U = run.U;
//
//  int n1; 
//  run.stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  run.stack_data.fetch("n2"   , &n2      );
//  int n1n2c; 
//  run.stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  run.stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  run.stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  run.stack_data.fetch("iu3"  , &n_flds  );
//
//  ComplexArray::size_type nc = n1n2c;
//
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
//void lcsolve::fftwReverseLayerofField ( int i_l, int i_f, stack& run) {
//
//  InputOutputArray& U = run.U;
//
//  int n1; 
//  run.stack_data.fetch("n1"   , &n1      );
//  int n2; 
//  run.stack_data.fetch("n2"   , &n2      );
//  int n1n2c; 
//  run.stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  run.stack_data.fetch("n1n2" , &n1n2    );
//  int n_layers; 
//  run.stack_data.fetch("iu2"  , &n_layers);
//  int n_flds;
//  run.stack_data.fetch("iu3"  , &n_flds  );
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
//void lcsolve::fftwForwardRaw( stack& run, RealArray& Rin, ComplexArray& Cout) {
//
//  int n1n2c; 
//  run.stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  run.stack_data.fetch("n1n2" , &n1n2    );
//  int iu2;
//  run.stack_data.fetch("iu2"  , &iu2     );
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
//void lcsolve::fftwReverseRaw( stack& run, ComplexArray& Cin, RealArray& Rout) {
//
//  int n1n2c; 
//  run.stack_data.fetch("n1n2c", &n1n2c   );
//  int n1n2; 
//  run.stack_data.fetch("n1n2" , &n1n2    );
//  int iu2;
//  run.stack_data.fetch("iu2"  , &iu2     );
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

void lcsolve::partialsInXandY( ComplexArray& U, RealArray& Ux, RealArray& Uy) {

//  unsigned usize  = U.capacity();
//
//  RealArray& kx   = run.kx;
//  RealArray& ky   = run.ky;
//  
//  unsigned kx_size = kx.capacity();
//
//  ComplexArray U_tmp(usize, czero);
//
//  int      n1n2c;
//  unsigned idk;
//
//  run.stack_data.fetch("n1n2c", &n1n2c   );              /* ~ number of complex elements in a layer       ~ */
//
//  for (unsigned k = 0; k < usize; k++) { 
//
//    if ( k % n1n2c == 0 ) { idk = 0; }
//    U_tmp[k] =  iunit * kx[idk] * U[k];
//    ++idk;
//
//  } /* ~ dU/dx -> ik_x U         ~ */
//
//  fftwReverseRaw( run, U_tmp,  Ux);                                          /* ~ transform to real space ~ */
//  for (unsigned k = 0; k < usize; k++) { U_tmp[k] =  iunit * ky[k] * U[k]; } /* ~ dU/dy -> ik_y U ~ */
//  fftwReverseRaw( run, U_tmp,  Uy);                                          /* ~ transform to real space ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::bracket( ComplexArray& BrKt, RealArray& d1x, RealArray& d1y, RealArray& d2x, RealArray& d2y) {

//  unsigned dsize = d1x.capacity();
//
//  RealArray B_tmp(dsize, zero);
//
//  for ( unsigned k = 0; k < dsize; k++ ) { B_tmp[k] = (d1y[k] * d2x[k]) - (d1x[k] * d2y[k]); }
//
//  fftwForwardRaw( run, B_tmp, BrKt);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double lcsolve::maxdU( RealArray& dx, RealArray& dy) {

  unsigned dsize  = dx.capacity();

  double test     = 0;
  double max      = 0;

   for (unsigned k = 0; k< dsize; k++ ) {

     test = ((dx[k] * dx[k]) + (dy[k] * dy[k]));
     if ( test > max ) { max = test; }
  }

  return max;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::averageAcrossLayers( int shift_sign, RealArray& dx, RealArray& dy) {

//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  int n1n2;
//  run.stack_data.fetch("n1n2",  &n1n2 );       /* ~ number of real elements in a layer       ~ */
//  int iu2;
//  run.stack_data.fetch("iu2"  , &iu2  );       /* ~ number of layers in stack                ~ */
//
//  int dsize = dx.capacity();
//
//  assert(dsize == n1n2 * iu2);
//
//  RealArray d_tmp_x(dsize, zero);
//  RealArray d_tmp_y(dsize, zero);
//
//  unsigned kstart    = n1n2;
//  unsigned kstop     = n1n2 * (iu2 - 1);
//  unsigned kshift    = shift_sign * n1n2;
//
//  unsigned idx       = 0;
//
//    for (unsigned k = kstart; k < kstop; k++) {
//
//          idx        =  k + kshift;
//
//        d_tmp_x[k]   = half * (dx[k] + dx[idx]) ;
//        d_tmp_y[k]   = half * (dy[k] + dy[idx]) ;
//
//    }
//
//    for (unsigned k = kstart; k < kstop; k++) {
//
//      dx[k] = d_tmp_x[k];
//      dy[k] = d_tmp_y[k];
//
//    }

}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::Step( std::string str_step ) {

//  std::string model;
//  run.palette.fetch(   "model", &model );          /* ~ either "hall" or "rmhd"              ~ */
//
//  int n1n2c;
//  run.stack_data.fetch("n1n2c" ,&n1n2c );          /* ~ number of complex elements per layer ~ */
//  int iu2;
//  run.stack_data.fetch("iu2"   ,&iu2   );          /* ~ number layers                        ~ */
//
//  double dt; 
//  run.stack_data.fetch("dt"    ,&dt    );          /* ~ current time increment               ~ */
//  double pfrac; 
//  run.stack_data.fetch("pfrac" ,&pfrac );          /* ~ fraction of dt to use for predictor- ~ */
//                                                   /* ~ step                                 ~ */
//  if (str_step.compare("predict") == 0 ) {         /* ~ use partial step in predictor case   ~ */
//    dt             = pfrac * dt;                   /* ~ dt is local so no problem here       ~ */
//  }
//
//  int kstart       = n1n2c;                        /* ~ stepping is only done for layers 1   ~ */
//  int kstop        = n1n2c * (iu2 - 1);            /* ~ through iu2 - 2                      ~ */
//                                                   /* ~ layers 0 and iu2 - 1 are for the     ~ */
//                                                   /* ~ boundaries and overlaps              ~ */
//  int idx;                                         /* ~ an index for the S's which are       ~ */
//                                                   /* ~ field-independent and thus the same  ~ */
//                                                   /* ~ across layers                        ~ */
//
//  for (unsigned k  = kstart; k < kstop; k++) {
//
//    if (k % kstart == 0 ) { idx = 0; }             /* ~ reset idx when starting new layer    ~ */
//
//      if (     str_step.compare("predict") == 0) { /* ~ the predictor case                   ~ */
//
//        tU0[k]     = (SE0[idx] * U0[k] + (dt * (B0[k] + D0[k] + A0[k]))) * SI0[idx];
//        tU1[k]     = (SE1[idx] * U1[k] + (dt * (B1[k] + D1[k] + A1[k]))) * SI1[idx];
//
//        if ( model.compare("hall") == 0 ) {
//
//          tU2[k]   = (SE2[idx] * U2[k] + (dt * (B2[k] + D2[k] + A2[k]))) * SI2[idx];
//          tU3[k]   = (SE3[idx] * U3[k] + (dt * (B3[k] + D3[k] + A3[k]))) * SI3[idx];
//
//        }
//      }
//      else if (str_step.compare("correct") == 0) { /* ~ the corrector case                   ~ */
//
//         U0[k]     = (SE0[idx] * U0[k] + (dt * (B0[k] + D0[k] + A0[k]))) * SI0[idx];
//         U1[k]     = (SE1[idx] * U1[k] + (dt * (B1[k] + D1[k] + A1[k]))) * SI1[idx];
//
//        if ( model.compare("hall") == 0 ) {
//
//           U2[k]   = (SE2[idx] * U2[k] + (dt * (B2[k] + D2[k] + A2[k]))) * SI2[idx];
//           U3[k]   = (SE3[idx] * U3[k] + (dt * (B3[k] + D3[k] + A3[k]))) * SI3[idx];
//
//        }
//      }
//
//    ++idx;
//
//  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

   /* ~ Destructor  ~ */

  lcsolve::~lcsolve( ) {

//    fftwFinalize();
//     destroyFields();

  }
