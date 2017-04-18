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
 *        FILE: Implementation of class "redhallmhd"
 *
 * DESCRIPTION: For defining and implementing the reduced mhd Hall physics for
 *              coronos. this class is responsible for "filling" lcsolve's data
 *              structures with the appropriate values - based on its
 *              "knowledge" of the physical model of * the plasma. These values
 *              are needed by lcsolve so that lcsolve can update its stack.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "cls_redhallmhd.hpp"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~ */
/* ~ Constructors ~ */
/* ~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

redhallmhd::redhallmhd() {

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

redhallmhd::redhallmhd(stack& run ) {

#ifndef HAVE_CUDA_H

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD,          &rank );
  int srun;  run.palette.fetch(         "srun",     &srun );
  int bdrys; run.palette.fetch(        "bdrys",    &bdrys );
  std::string model; run.palette.fetch("model",    &model );
  std::string init; run.palette.fetch("initMode",  &init  );

  init_physics_data(     run   );                          /* ~ physics - specific parameters               ~ */
  initU(                 run   );                          /* ~ initialization of layers 1 - n3 of U        ~ */
                                                           /* ~ for srun > 1 AUX has real space O and J now ~ */
  initTimeInc(           run   );
  fftw.fftwForwardAll(   run   );                          /* ~ puts real-space fields into Fourier space   ~ */

  initBoundaries(        run   );                          /* ~ initialization of quantities needed for     ~ */

  if (bdrys == 0) { applyBC("predict", run); }
  
  if ((srun > 1) || (init.compare("from_data") == 0)) { fftw.fftwForwardAll(   run, O, J);}

  applyBC( "predict",    run   );
  OfromP(                run   );                          /* ~ the data files contain both stream func. P  ~ */
  HfromA(                run   );

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

//if (srun == 1) { checkState(4,      run, "c"); }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

  evalValf(              run   );
  evalUmean(             run   );
  evalElls(              run   );

  if (srun == 1) { fftw.fftwReverseAll(   run, O, J);
                   run.writeUData  (               ); }    /* ~ the real space values of O and J as well as ~ */
                                                           /* ~ P and A in a "subrun zero" data file        ~ */
}

#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef HAVE_CUDA_H

void redhallmhd::initTimeInc( stack& run ){

 int n_flds;
 run.stack_data.fetch("iu3" , &n_flds);

 maxU.assign(n_flds, zero);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initU( stack& run ) {

  fftw.fftwInitialize( run );

  MPI_Barrier(MPI_COMM_WORLD);

  int srun;             run.palette.fetch("srun",     &srun    );
  std::string scenario; run.palette.fetch("scenario", &scenario);

  if (scenario.compare("reconnection")  == 0) {

    if (srun == 1) {

      std::string init; run.palette.fetch("initMode", &init);

      if (init.compare("fourierspace") == 0) computeFourierU( run );
      if (init.compare("realspace")    == 0) computeRealU(    run );
      if (init.compare("from_data" )   == 0) readUData(       run );

      int ilnr; run.palette.fetch("ilnr", &ilnr);

      if (ilnr != 0 && init.compare("from_data") !=0 ) pLinzEnv( run );

    }
    else           { readUData(   run ); }
  }
  else if (scenario.compare("parker") == 0 ) {

    if (srun == 1) {

      std::string init; run.palette.fetch("initMode", &init);
      if (init.compare("from_data" )  == 0 ) readUData(       run );
      else { run.zeroU(); }

    }
    else { readUData( run ); }

  }
  else if (scenario.compare("gauss")  == 0 ) {

    if (srun == 1) { initGauss( run ); }
    else           { readUData( run ); }

  }

  int n3;      run.palette.fetch("p3",       &n3    );
  int iu2;     run.stack_data.fetch("iu2",   &iu2   );
  int np;      run.palette.fetch("np",       &np    );
  int n1n2c;   run.stack_data.fetch("n1n2c", &n1n2c );
  

  int usize = iu2*n1n2c;

  P.assign(usize,czero);
  A.assign(usize,czero);
  O.assign(usize,czero);
  J.assign(usize,czero);

  int calcqvz; run.palette.fetch("calcqvz", &calcqvz);
  int calcsvz; run.palette.fetch("calcsvz", &calcsvz);

  if (calcqvz == 1) {

    RealArray Qrecord; Qrecord.assign(26,zero); QtyVsZ.assign(np*n3, Qrecord);

  }

  if (calcsvz == 1) {

    int isp;    physics_data.fetch(    "isp", &isp);
    int n1;     run.stack_data.fetch(   "n1", &n1 );
    int n2;     run.stack_data.fetch(  "n2",  &n2 );
    RealVar kc; run.palette.fetch(     "kc",  &kc );

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank );

    RealArray& k2=run.k2;

    int size_k2=k2.size();
    for (unsigned k = 0; k< size_k2-1;++k ) { if (k > 0 && std::sqrt(k2[k]) >= kc) { kb = std::sqrt(k2[k]); break;}}


//  dk     = sqrt(two)*(pi * n1) / isp;
//
    dk     = (pi*n1) / ((RealVar) isp);

    dk_m1  = one / dk;
    kf     = kb + (isp * dk);

//  kb     = two * pi;

    ikb    = 1 + ((int) (kb * dk_m1));
    ikf    = 1 + ((int) (kf * dk_m1));
//  ikb    = ((int) (kb * dk_m1));
//  ikf    = ((int) (kf * dk_m1));

    nk     = ikf - ikb + 1;

//  std::cout << "initU: ikb = " << ikb << std::endl;
//  std::cout << "initU: ikf = " << ikf << std::endl;
//  std::cout << "initU: isp = " << isp << std::endl;
//  std::cout << "initU: nk  = " << nk  << std::endl;
//  std::cout << "initU: kb  = " << kb  << std::endl;
//  std::cout << "initU: kf  = " << kf  << std::endl;
//  std::cout << "initU: dk  = " << dk  << std::endl;

    RealArray   SpcRecord; SpcRecord.assign(11,zero);
    Real2DArray SpcLayer;  SpcLayer.assign(n3,SpcRecord);
    SpcVsZ.assign(isp+1,SpcLayer);
    ke.assign(isp+1,zero);

    for (unsigned k = 0; k < nk; ++k) { ke[k] = log((k+ikb)*dk); }

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::init_physics_data( stack& run ) {

  std::string model;  run.palette.fetch("model",   &model  ); 
  int         bdrys;  run.palette.fetch("bdrys",   &bdrys  );
  RealVar     tstart; run.palette.fetch("tstart" , &tstart );
  RealVar     nu;     run.palette.fetch("nu"     , &nu     );
  RealVar     eta;    run.palette.fetch("eta"    , &eta    );

  RealVar qs0   = nu;
  RealVar qs1   = eta;

  RealVar delta;  
  RealVar epratio;
  RealVar beta;   
  RealVar kappa;  
  RealVar qs2;
  RealVar qs3;
  RealVar ssqd;
  RealVar rho;

  /* ~ hall - related ~ */

  if ( model.compare("hall") == 0 )  {           /* ~ initialize only if needed  ~ */

    run.palette.fetch("delta"  , &delta   );
    run.palette.fetch("epratio", &epratio );
    run.palette.fetch("beta"   , &beta    );
    run.palette.fetch("kappa"  , &kappa   );

    qs2         = kappa + (half * beta * eta);
    qs3         = nu;
    ssqd        = two * delta * sqrt(epratio);
    rho         = sqrt(beta) * delta;

  }

  std::string pname;
  std::string padjust;

  padjust.assign("adj"  );

  pname.assign(  "t_cur");   physics_data.emplace(pname, tstart,  padjust);
  pname.assign(  "dtvb" );   physics_data.emplace(pname, zero,    padjust);

  int calcsvz; run.palette.fetch("calcsvz", &calcsvz);

  if (calcsvz == 1 ) {

    int     n1;  run.stack_data.fetch("n1", &n1);
    int     n2;  run.stack_data.fetch("n2", &n2);
    int isp = n1/2;
       

    padjust.assign("rfx");
    pname.assign("isp"); physics_data.emplace(pname, isp, padjust);

  }

  if (bdrys > 0) {

    int brcount = 0;
    int trcount = 0;

    padjust.assign("adj");

    pname.assign("brcount"); physics_data.emplace(pname, brcount, padjust);
    pname.assign("trcount"); physics_data.emplace(pname, trcount, padjust);

  }

  padjust.assign("rfx");

  pname.assign( "qs0" );     physics_data.emplace(pname, qs0,     padjust);
  pname.assign( "qs1" );     physics_data.emplace(pname, qs1,     padjust);

  if ( model.compare("hall") == 0 )  {
  
    pname.assign( "qs2" );   physics_data.emplace(pname, qs2,     padjust);
    pname.assign( "qs3" );   physics_data.emplace(pname, qs3,     padjust);

    pname.assign( "ssqd");   physics_data.emplace(pname, ssqd,    padjust);
    pname.assign( "rho");    physics_data.emplace(pname, rho,     padjust);

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::computeRealU( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1;     run.stack_data.fetch( "n1",    &n1    );
  int n2;     run.stack_data.fetch( "n2",    &n2    );
  int n3;     run.stack_data.fetch( "n3",    &n3    );
  int n1n2;   run.stack_data.fetch( "n1n2",  &n1n2  );
  int n_flds; run.stack_data.fetch("iu3",    &n_flds);

  int n_lyrs              = n3;

  InputOutputArray& U     = run.U;
  RealArray&        x     = run.x;
  RealArray&        y     = run.y;

  int idx                 = 0;

  for (int i_f = 0; i_f < n_flds; ++i_f) {

    switch(i_f) {

    case(0) :
      for (int i_x=0;i_x < n1; ++i_x) {
        for (int j_y=0;j_y < n1; ++j_y) {

          idx             = (i_x * n1) + j_y;

          U[idx][n3][i_f] = - 0.0032L * ( cos(two_pi*x[i_x]) - cos(two_pi * y[j_y]) );

        }
      }
      break;
    case(1) :
      for (int i_x=0;i_x < n1; ++i_x) {
        for (int j_y=0;j_y < n1; ++j_y) {
 
          idx             = (i_x * n1) + j_y;
          U[idx][n3][i_f] = four * 0.1L * sin(two_pi *x[i_x]) * sin(two_pi * y[j_y]);
 
        }
      }

      break;
    case(2) :

      break;
    case(3) :

      break;

    }
  }

  for (  int i_f = 0; i_f < n_flds;     ++i_f) {
    for (int i_l = 1; i_l < n_lyrs + 1; ++i_l) {
    
    for (int   k = 0; k< n1n2; ++k) { U[k][i_l][i_f] = U[k][n3][i_f]; }

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initGauss( stack& run ) {

  int               rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int               n1;      run.stack_data.fetch( "n1"     , &n1     );
  int               n2;      run.stack_data.fetch( "n2"     , &n2     );
  int               n3;      run.stack_data.fetch( "n3"     , &n3     );
  int               n1n2;    run.stack_data.fetch( "n1n2"   , &n1n2   );
  int               n_flds;  run.stack_data.fetch( "iu3"    , &n_flds );

  RealVar           zl;      run.palette.fetch(    "zl"     , &zl     );
  RealVar           r_x;     run.palette.fetch(    "r_x"    , &r_x    );
  RealVar           r_y;     run.palette.fetch(    "r_y"    , &r_y    );

  RealVar           sigma_x; run.palette.fetch(    "sigma_x", &sigma_x);
  RealVar           sigma_y; run.palette.fetch(    "sigma_y", &sigma_y);
  RealVar           sigma_z; run.palette.fetch(    "sigma_z", &sigma_z);

  int               n_lyrs       = n3;

  InputOutputArray& U            = run.U;

  RealArray&        x            = run.x;
  RealArray&        y            = run.y;
  RealArray&        z            = run.z;

  int               idx          = 0;

  RealVar kappa_x;    kappa_x    = one / (two*sigma_x*sigma_x);
  RealVar kappa_y;    kappa_y    = one / (two*sigma_y*sigma_y);
  RealVar kappa_z;    kappa_z    = one / (two*sigma_z*sigma_z);

  RealVar mu_x;       mu_x       = kappa_x / four;
  RealVar mu_y;       mu_y       = kappa_y / four;

  RealVar varsigma_x; varsigma_x = r_x * kappa_x;
  RealVar varsigma_y; varsigma_y = r_y * kappa_y;

  RealVar psi_xzero; psi_xzero   = acos( kappa_x / (sqrt(pow(varsigma_x,2)+pow(kappa_x,2))));
  RealVar psi_yzero; psi_yzero   = acos( kappa_y / (sqrt(pow(varsigma_y,2)+pow(kappa_y,2))));

  RealVar xfld_zero; xfld_zero   = exp(-mu_x)*sin((varsigma_x / four) + psi_xzero);
  RealVar yfld_zero; yfld_zero   = exp(-mu_y)*sin((varsigma_y / four) + psi_yzero);

  RealVar eenv_zero; run.palette.fetch("eenv_zero", &eenv_zero);

  RealVar xfld;
  RealVar yfld;
  RealVar zenv;
  RealVar phi_norm;

  phi_norm                       = one / sqrt((pow(varsigma_x,2) + pow(kappa_x,2))*(pow(varsigma_y,2) + pow(kappa_y,2)));

  if (rank == 0) {

    std::cout << "sigma_x    = " << sigma_x    << std::endl;
    std::cout << "sigma_y    = " << sigma_y    << std::endl;
    std::cout << "" << std::endl;
    std::cout << "kappa_x    = " << kappa_x    << std::endl;
    std::cout << "kappa_y    = " << kappa_y    << std::endl;
    std::cout << "" << std::endl;
    std::cout << "varsigma_x = " << varsigma_x << std::endl;
    std::cout << "varsigma_y = " << varsigma_y << std::endl;
    std::cout << "" << std::endl;
    std::cout << "mu_x       = " << mu_x       << std::endl;
    std::cout << "mu_y       = " << mu_y       << std::endl;
    std::cout << "" << std::endl;
    std::cout << "r_x        = " << r_x        << std::endl;
    std::cout << "r_y        = " << r_y        << std::endl;
    std::cout << "" << std::endl;
    std::cout << "psi_xzero  = " << psi_xzero  << std::endl;
    std::cout << "psi_yzero  = " << psi_yzero  << std::endl;
    std::cout << "" << std::endl;
    std::cout << "phi_norm   = " << phi_norm   << std::endl;

  }

  RealVar next_u;
  for (int i_f = 0; i_f < n_flds; ++i_f) {

    switch(i_f) {

    case(0) :
    for (unsigned k_z=0; k_z < n3; ++k_z) {
      zenv                     = eenv_zero * pow(sin((pi/zl) * (z[k_z] - zl)),2) * exp(-kappa_z * pow((z[k_z] - (half*zl)),2));
      for (unsigned i_x=0; i_x < n1; ++i_x) {
        xfld                   = exp(-kappa_x * pow(( x[i_x] - half),2) ) * sin(varsigma_x * pow((x[i_x] - half),2) + psi_xzero ); 
        for (unsigned j_y=0; j_y < n2; ++j_y) {
          yfld                 = exp(-kappa_y * pow(( y[j_y] - half),2) ) * sin(varsigma_y * pow((y[j_y] - half),2) + psi_yzero ); 
          idx                  = (i_x * n1) + j_y;
          next_u               = phi_norm * zenv * (xfld - xfld_zero) * (yfld - yfld_zero);
          if (std::abs(next_u) >= teensy) {
            U[idx][k_z+1][i_f] = next_u;
          }
          else {
            U[idx][k_z+1][i_f] = zero;
          }
        }
      }
    }

      break;
    case(1) :
      break;
    case(2) :

      break;
    case(3) :
      break;
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::computeFourierU( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1;    run.stack_data.fetch("n1",    &n1    );
  int n2;    run.stack_data.fetch("n2",    &n2    );
  int n3;    run.stack_data.fetch("n3",    &n3    );
  int n1n2;  run.stack_data.fetch("n1n2",  &n1n2  );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int iu1;   run.stack_data.fetch("iu1",   &iu1   );
  int iu2;   run.stack_data.fetch("iu2",   &iu2   );
  int iu3;   run.stack_data.fetch("iu3",   &iu3   );

  int n_flds          = iu3;
  int n_lyrs          = n3;

  ComplexArray Cin(n1n2c, czero);
  RealArray    Rout(n1n2, zero);
  
  InputOutputArray& U = run.U;

  RealVar               real_part;
  RealVar               imag_part;

  ComplexVar tuple;

  unsigned idx        = 0;

  for (int i_f = 0; i_f < n_flds; ++i_f) {

     for (unsigned k  = 0; k < n1n2c; ++k) { Cin[k]    = czero; }
     for (unsigned k  = 0; k < n1n2;  ++k) { Rout[k]   =  zero; }
     
     switch(i_f) {

     case(0) :
       idx            =        0 * (n2/2 + 1) + 1;
       real_part      =  1.0e-03L;
       imag_part      =  0.0L;
       tuple          = ComplexVar(real_part, imag_part);
       Cin[idx]       = tuple;

       idx            =        1 * (n2/2 + 1);
       real_part      = -1.0e-03L;
       imag_part      =  0.0L;
       tuple          = ComplexVar(real_part, imag_part);
       Cin[idx]       = tuple;

       idx            = (n1-1) * (n2/2 + 1);
       Cin[idx]       = tuple;

       break;
     case(1) :
       idx            =       1 * (n2/2 + 1) + 1;
       real_part      =    -1.0e-01L;
       imag_part      =     0.0L;
       tuple          = ComplexVar(real_part, imag_part);
       Cin[idx]       = tuple;

       idx            =  (n1-1) * (n2/2 + 1) + 1;
       real_part      =     1.0e-01L;
       imag_part      =     0.0L;
       tuple          = ComplexVar(real_part, imag_part);
       Cin[idx]       = tuple;

       break;
     case(2) :
       /* ~ edit for non-zero initial bz ~ */
       break;
     case(3) :
       /* ~ edit for non-zero initial vz ~ */
       break;

    }

    fftw.fftwReverseIC(Cin, Rout);
    for (unsigned k   = 0; k < n1n2; ++k) { U[k][n3][i_f] = Rout[k]; }

  }

  for (     int i_f = 0; i_f < n_flds; ++i_f) {
    for (   int i_l = 1; i_l < n3;     ++i_l) {
      for ( int k   = 0;   k < n1n2;   ++k)   { U[k][i_l][i_f] = U[k][n3][i_f]; }
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::readUData( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int srun; run.palette.fetch("srun", &srun);
  std::string data_dir; run.palette.fetch("data_dir", &data_dir);

  std::string data_file;
  data_file                = "./" + data_dir + "/" + run.getLastDataFilename(srun-1);

  const char *c_data_file                = data_file.c_str();

  std::ifstream ifs;
  ifs.open( c_data_file, std::ios::in );

  if ( ifs.good() ) {

    InputOutputArray& U                  = run.U;
    InputOutputArray& AUX                = run.AUX;
    
    int iu3;  run.stack_data.fetch("iu3",  &iu3 );
    int n1;   run.stack_data.fetch("n1",   &n1  );
    int n2;   run.stack_data.fetch("n2",   &n2  );
    int n3;   run.stack_data.fetch("n3",   &n3  );
    int n1n2; run.stack_data.fetch("n1n2", &n1n2);

    int n_slab_points                    = n1n2 * iu3;
    int point_count                      = 0;
    int slab_index                       = 1;
    int from_col_maj_idx                 = 0;
    int to_row_maj_idx                   = 0;
    int i                                = 0;
    int j                                = 0;

    RealVar next_p; 
    RealVar next_a; 
    RealVar next_bz; 
    RealVar next_vz;

    RealVar next_o;
    RealVar next_j;
    
    while ( !ifs.eof() ) {

      if (slab_index > n3) break; 

      ifs >> next_p;
      ++point_count;
      ifs >> next_a;
      ++point_count;
//    if (next_p < teensy) { next_p = zero;}
//    if (next_a < teensy) { next_a = zero;}

      U[to_row_maj_idx][slab_index][0]   = next_p;
      U[to_row_maj_idx][slab_index][1]   = next_a;

      if(iu3 <= 2) {
        ifs >> next_o;
        ifs >> next_j;
//      if (next_o < teensy) { next_o = zero;}
//      if (next_j < teensy) { next_j = zero;}

        AUX[to_row_maj_idx][slab_index][0] = next_o;
        AUX[to_row_maj_idx][slab_index][1] = next_j;

      }
      else if(iu3 > 2) {

        ifs >> next_bz;
        ++point_count;
        ifs >> next_vz;
        ++point_count;

        ifs >> next_o;
        ifs >> next_j;

        U[to_row_maj_idx][slab_index][2] = next_bz;
        U[to_row_maj_idx][slab_index][3] = next_vz;

        AUX[to_row_maj_idx][slab_index][0] = next_o;
        AUX[to_row_maj_idx][slab_index][1] = next_j;

      }

      if (from_col_maj_idx < n1n2) {

        ++from_col_maj_idx;
        if (from_col_maj_idx % n2 != 0) ++j;
        else {
                                       j = 0;
                                       ++i;
        }
      }

      if (to_row_maj_idx < n1n2 - 1) to_row_maj_idx = i + (j*n1);
      else to_row_maj_idx                = 0;

      if (from_col_maj_idx == n1n2) {

        from_col_maj_idx                 = 0;
        i                                = 0;
        j                                = 0;
      }

      if(point_count == n_slab_points) {

        point_count                      = 0;
        ++slab_index;
      }
    }

    ifs.close();
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::pLinzEnv( stack& run ) {

  int    n1n2; run.stack_data.fetch("n1n2", &n1n2);
  int    n3;   run.stack_data.fetch("n3",   &n3  );
  RealVar zl;  run.palette.fetch(   "zl",   &zl  );

  InputOutputArray&  U  = run.U;
  RealArray&         z  = run.z;

  int i, j;

  for (i = 0; i < n1n2; ++i) {
    for (j = 1; j < n3+1; ++j) {

      U[i][j][0]        = U[i][j][0] * (1.0 - (std::abs(z[j] - (0.5*zl))/(0.5*zl)));

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initBoundaries( stack& run) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);  

  int np; run.palette.fetch("np", &np);
  int srun; run.palette.fetch("srun", &srun);

  if ( rank == 0 || rank == np - 1 ) {
    
    int bdrys; run.palette.fetch("bdrys", &bdrys);

    if (bdrys > 0 ) { initFootPointDriving( run ); }
    else            { initNoDrive(          run ); }

  }

  /* ~ possibly want to put an MPI_Barrier here so the stacks 
   *   with exterior layers don't get out of sync with the others ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::finalizeBoundaries( stack& run) {

    int bdrys; run.palette.fetch("bdrys", &bdrys);

    if (bdrys > 0) { finalizeFootPointDriving(   run ); }
    else           { finalizeLineTiedBoundaries( run ); }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::countModes( stack& run ) {

  RealVar kc; run.palette.fetch("kc", &kc);

  bool   l_reset_success = false;

  int m                  = 0;
  while ( ((RealVar) m) <= kc ) { ++m; } 

  RealVar arg;
  int nf                 = 0;
  for (   int ix = -m; ix < (m + 1);  ++ix) {
    for ( int iy = -m; iy < (m + 1);  ++iy) {

    arg                  = sqrt( pow((two_pi * ((RealVar) ix)),2) + pow((two_pi * ((RealVar) iy)),2) );

    if (arg != zero && arg <= kc) {   ++nf; }
    
    }
  }
 
  l_reset_success        = run.palette.reset("nf", nf);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initFootPointDriving( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  MPI_Status status;

  int tag_old           = 0;
  int tag_new           = 1;
  int tag_ptop          = 2;

  countModes( run );

  int       nf; run.palette.fetch("nf"  , &nf  );
  RealVar tauC; run.palette.fetch("tauC", &tauC);
  RealVar tauE; run.palette.fetch("tauE", &tauE);

  RealVar dtau          = ((RealVar) (1.383)) * tauC;
  RealVar qfp           = zero;
  RealVar ffp           = (two_pi / tauE) * (one / sqrt( (RealVar) nf));

  bool l_reset_success  = false;

  l_reset_success       = run.palette.reset("dtau", dtau);
  l_reset_success       = run.palette.reset("qfp" , qfp);
  l_reset_success       = run.palette.reset("ffp" , ffp);

  const int    seed     = 1234567;
  srand(seed);

  int rcount; run.palette.fetch( "rcount",  &rcount);

  RealVar dummy;
  if (rcount > 0) { for (int l = 0; l < rcount; ++l) { dummy = (double) rand() / RAND_MAX; }}

  int srun;   run.palette.fetch(     "srun", &srun  );
  int n1;     run.stack_data.fetch(   "n1",  &n1    );
  int n2;     run.stack_data.fetch(   "n2",  &n2    );
  int n1n2c;  run.stack_data.fetch( "n1n2c", &n1n2c );
  RealVar kc; run.palette.fetch(       "kc", &kc    );
  int np;     run.palette.fetch(       "np", &np    );
  int bdrys;  run.palette.fetch(    "bdrys", &bdrys );
  int n3;     run.palette.fetch(       "p3", &n3    );

  int i_ua;
  int i_cpy;

  ComplexArray Ptop(n1n2c, czero);

  RealArray& k2         = run.k2;

  RealVar    next_real; 
  RealVar    next_imag;
  ComplexVar tuple;

  i_ua = 0;

  if (rank == 0 || rank == np - 1 ) {

    roldlb.assign(n1n2c,czero);
    rnewlb.assign(n1n2c,czero);
    roldub.assign(n1n2c,czero);
    rnewub.assign(n1n2c,czero);

  }

  if (rank == 0) {

    int brcount; physics_data.fetch("brcount", &brcount);

    if (srun == 1) {

      for (int l = 0; l < n1n2c; ++l) {

/* ~ fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  */

        if (l  < (n1n2c - (n2/2)-1)) {
          if ((l == 0) || (l %((n2/2)+1) !=0) ){

            if ( sqrt(k2[l]) < kc ) {

                next_real   = ffp * ( ((double) rand() / RAND_MAX ) * two - one); ++brcount;
                dummy       =         ((double) rand() / RAND_MAX ) * two - one;  ++brcount;
                next_imag   = ffp * ( ((double) rand() / RAND_MAX ) * two - one); ++brcount;
                dummy       =         ((double) rand() / RAND_MAX ) * two - one;  ++brcount;

                tuple       = ComplexVar(next_real, next_imag);
                rnewlb[l]   = tuple;

            }
            if ( l == 0) { rnewlb[l] = czero; }
            ++i_ua;
          } // not copying
          else { i_cpy    = (l/((n2/2)+1));
//          if (i_cpy < 0 || i_cpy > n2/2-1) {
//            std::cout << "initFoot: WARNING ONE - i_cpy = " << i_cpy << " for l = " << l << std::endl;
//          } // i_cpy out of range
//          else { old
              if ( i_cpy < n2/2 ) {
                rnewlb[l] = rnewlb[i_cpy];
//              rnewlb[l] =  std::conj(rnewlb[i_cpy]);
              } // really copying
              else {
                if ( sqrt(k2[l]) < kc ) {

                  next_real   = ffp * ( ((double) rand() / RAND_MAX ) * two - one); ++brcount;
                  dummy       =         ((double) rand() / RAND_MAX ) * two - one;  ++brcount;
                  next_imag   = ffp * ( ((double) rand() / RAND_MAX ) * two - one); ++brcount;
                  dummy       =         ((double) rand() / RAND_MAX ) * two - one;  ++brcount;

                  tuple       = ComplexVar(next_real, next_imag);
                  rnewlb[l]   = tuple;

                }
                ++i_ua;
              } // not really copying
//          }   // i_cpy in range
          }     // copying
        }       // l is <  (n1n2c - (n2/2) - 1
        else {

          i_cpy = (l - n1n2c + n2/2+2);
          if (i_cpy < 0 || i_cpy > n2/2+1) {
            std::cout << "initFoot: WARNING TWO - i_cpy = "             << i_cpy << " for l = " << l <<  std::endl;
          }
          else {
//          std::cout << "initFoot: l >= n1n2c - n2/2 <-> i_cpy = " << i_cpy << " for l = " << l << std::endl;
            if ( l != (n1n2c - 1)) {
              rnewlb[l] = std::conj(rnewlb[i_cpy]);
            }
            else { rnewlb[l] = czero;}
          }
        }       // l is >= (n1n2c - (n2/2) - 1

/* ~ fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  */

      }        // for loop in l

      std::cout << "initFoot: i_ua = " << i_ua << " for rank = " << rank << std::endl; 

    l_reset_success     = physics_data.reset("brcount", brcount);

   }           // srun = 1
   else {
      
      std::string data_dir;  run.palette.fetch(    "data_dir",   &data_dir  );
      std::string prefix;    run.palette.fetch(    "prefix",     &prefix    );
      std::string run_label; run.palette.fetch(    "run_label",  &run_label );
      std::string res_str;   run.stack_data.fetch( "res_str",    &res_str   );

      std::string srn_str              = static_cast<std::ostringstream*>( &(std::ostringstream() << ( srun - 1 ) ) ) -> str();
      std::string boundary_data_file   = data_dir + "/" + prefix + '_' + res_str + "r" + srn_str;

      const char *c_boundary_data_file = boundary_data_file.c_str();

      std::ifstream ifs;
      ifs.open( c_boundary_data_file, std::ios::in );

      if ( ifs.good() ) {

        int fld_count                  = 0;
        int line_count                 = 0;

        while ( !ifs.eof() && line_count < n1n2c ){

          if (fld_count == 6) { ++line_count; fld_count = 0; }

          ifs >> next_real;
          ifs >> next_imag;

          if (ifs.eof()) { break;}

          tuple=ComplexVar(next_real, next_imag);
          switch(fld_count) {

          case(0) : roldlb[line_count] = tuple; ++fld_count; break;
          case(1) : rnewlb[line_count] = tuple; ++fld_count; break;
          case(2) : roldub[line_count] = tuple; ++fld_count; break;
          case(3) : rnewub[line_count] = tuple; ++fld_count; break;
          case(4) : P[line_count]      = tuple; ++fld_count; break;
          case(5) : Ptop[line_count]   = tuple; ++fld_count; break;

          default : std::cout << "initFootPointDriving: WARNING - fld_count maximum exceeded." << std::endl;

          }
        }

        ifs.close();

       if ( bdrys == 2 ) {
         MPI_Send( &roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np - 1, tag_old,  MPI_COMM_WORLD );
         MPI_Send( &rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np - 1, tag_new,  MPI_COMM_WORLD );
         MPI_Send( &Ptop.front(),   n1n2c, MPI::DOUBLE_COMPLEX, np - 1, tag_ptop, MPI_COMM_WORLD );
       }
     }
     else { std::cout << "initFootPointDriving: WARNING - could not open file " << boundary_data_file << std::endl; }

    } // srun > 1
  }   //rank is zero

  if ( rank == np - 1 ) {
    if (bdrys  == 2) {

      int trcount; physics_data.fetch("trcount", &trcount);

      if (srun == 1) {
        for (int l = 0; l < n1n2c; ++l) {

          if (l  < (n1n2c - (n2/2)-1)) {
            if ((l == 0) || (l %((n2/2)+1) !=0) ) {

/* ~ fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  */

              if ( sqrt(k2[l]) < kc ) {

                  dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                  dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                  dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                  next_real   = ffp * (((double) rand() / RAND_MAX ) * two - one); ++trcount;
                  dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                  next_imag   = ffp * (((double) rand() / RAND_MAX ) * two - one); ++trcount;
                  dummy       =        ((double) rand() / RAND_MAX );              ++trcount;

                  tuple       = ComplexVar(next_real, next_imag);
                  rnewub[l]   = tuple;

              }
              if ( l == 0) { rnewub[l] = czero; }
              ++i_ua;
            } // not copying
            else { i_cpy    = (l/((n2/2)+1));
//            if (i_cpy < 0 || i_cpy > n2/2-1) {
//              std::cout << "initFoot: WARNING THREE - i_cpy = " << i_cpy << " for l = " << l << std::endl;
//            } // i_cpy out of range
//            else {
                if ( i_cpy < n2/2 ) {
                  rnewub[l] = rnewub[i_cpy];
                } // really copying
                else {

                  if ( sqrt(k2[l]) < kc ) {

                      dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                      dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                      dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                      next_real   = ffp * (((double) rand() / RAND_MAX ) * two - one); ++trcount;
                      dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
                      next_imag   = ffp * (((double) rand() / RAND_MAX ) * two - one); ++trcount;
                      dummy       =        ((double) rand() / RAND_MAX );              ++trcount;

                      tuple       = ComplexVar(next_real, next_imag);
                      rnewub[l]   = tuple;

                  }

                  ++i_ua;
                }  // not really copying
//            }    // i_cpy in range
            }      // copying
          } // l  < (n1n2c - (n2/2)-1)
          else {

            i_cpy = (l - n1n2c + n2/2+2);
            if (i_cpy < 0 || i_cpy > n2/2+1) {
              std::cout << "initFoot: WARNING FOUR - i_cpy = "             << i_cpy << " for l = " << l <<  std::endl;
            }
            else {
//            std::cout << "initFoot: l >= n1n2c - n2/2 <-> i_cpy = " << i_cpy << " for l = " << l << std::endl;
              if (l != (n1n2c - 1)) {
              rnewub[l] = std::conj(rnewub[i_cpy]);
              }
              else { rnewub[l] = czero; }
            }
          } //l  >= (n1n2c - (n2/2)-1)

/* ~ fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  */

        } // end for loop in l

//    assert (n1n2c == i_ua + n1);
        std::cout << "initFoot: i_ua = " << i_ua << " for rank = " << rank << std::endl; 
      l_reset_success   = physics_data.reset("trcount", trcount);

      } // srun = 1 new
      else {
        for (int l = 0; l < n1n2c; ++l) {
          if ( sqrt(k2[l]) < kc ) {
  
              dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
              dummy       =        ((double) rand() / RAND_MAX );              ++trcount;

          } // k < kc for srun > 1
        }

         MPI_Recv(&roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_old,  MPI_COMM_WORLD, &status);
         MPI_Recv(&rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_new,  MPI_COMM_WORLD, &status);
         MPI_Recv(&Ptop.front(),   n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_ptop, MPI_COMM_WORLD, &status);

         for (unsigned k = 0; k < n1n2c; ++k) { P[(n3*n1n2c) + k] = Ptop[k]; }

      }  // srun > 1  (else) new
    }    // bdrys = 2
  }      // rank is np - 1
}        // end initFoot

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::finalizeFootPointDriving( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status   status;

  int  rcount;  run.palette.fetch(  "rcount",  &rcount  );
  int brcount;  physics_data.fetch( "brcount", &brcount );
  int trcount;  physics_data.fetch( "trcount", &trcount );

  rcount              = rcount + ( (brcount) > (trcount) ? (brcount) : (trcount) );

  run.palette.reset("rcount", rcount);

  int n3;     run.palette.fetch(    "p3",        &n3    );
  int np;     run.palette.fetch(    "np",      &np      );
  int n1n2c;  run.stack_data.fetch( "n1n2c",   &n1n2c   );
  int srun;   run.palette.fetch(    "srun",    &srun    );

  ComplexArray Ptop(n1n2c, czero);

  int tag_old         = 0;
  int tag_new         = 1;
  int tag_ptop        = 2;
  
  if (rank            == 0)      {
     
    MPI_Recv(&roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np -1, tag_old,  MPI_COMM_WORLD, &status);
    MPI_Recv(&rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np -1, tag_new,  MPI_COMM_WORLD, &status);
    MPI_Recv(&Ptop.front(),   n1n2c, MPI::DOUBLE_COMPLEX, np -1, tag_ptop, MPI_COMM_WORLD, &status);

    std::string data_dir;  run.palette.fetch(    "data_dir",   &data_dir  );
    std::string prefix;    run.palette.fetch(    "prefix",     &prefix    );
    std::string run_label; run.palette.fetch(    "run_label",  &run_label );
    std::string res_str;   run.stack_data.fetch( "res_str",    &res_str   );

    std::string srn_str              = static_cast<std::ostringstream*>( &(std::ostringstream() << ( srun ) ) ) -> str();
    std::string boundary_data_file   = data_dir + "/" + prefix + '_' + res_str + "r" + srn_str;
    const char *c_boundary_data_file = boundary_data_file.c_str();

    std::ofstream ofs;
    ofs.open( c_boundary_data_file, std::ios::out | std::ios::trunc );

    if ( ofs.good() ) {
      for (unsigned k = 0; k < n1n2c; k++) {               /* ~ order: roldlb rnewlb roldub rnewub pbot ~ */

        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << roldlb[k].real() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << roldlb[k].imag() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << rnewlb[k].real() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << rnewlb[k].imag() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << roldub[k].real() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << roldub[k].imag() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << rnewub[k].real() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << rnewub[k].imag() << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << P[k].real()      << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << P[k].imag()      << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << Ptop[k].real()   << " ";
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << Ptop[k].imag()   << " ";
        ofs << std::endl;

      }
      ofs.close();
    }
    else { std::cout << "finalizeFootPointDriving: WARNING - could not open file " << boundary_data_file << std::endl; }
  }
  else if (rank       == np - 1) {

    assert (roldub.size() == n1n2c);
    assert (rnewub.size() == n1n2c);

    for (unsigned k = 0; k < n1n2c; ++k) { Ptop[k] = P[(n3*n1n2c) + k]; }

    MPI_Send( &roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_old,  MPI_COMM_WORLD );
    MPI_Send( &rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_new,  MPI_COMM_WORLD );
    MPI_Send( &Ptop.front(),   n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_ptop, MPI_COMM_WORLD );

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::finalizeLineTiedBoundaries( stack& run ) {

   applyBC("finalize", run);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initNoDrive( stack& run) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  int srun; run.palette.fetch(   "srun" ,  &srun );
  int np;   run.palette.fetch(    "np"  ,  &np   );      /* ~ number of processes ~ */
  int n1n2; run.stack_data.fetch( "n1n2",  &n1n2 );
  int n3;   run.stack_data.fetch( "n3"  ,  &n3   );      /* ~ number of layers    ~ */
  int iu2;  run.stack_data.fetch("iu2",    &iu2   );     /* ~ n3 + 2              ~ */
  int iu3;  run.stack_data.fetch("iu3"  ,  &iu3   );     /* ~ number of fields    ~ */

  InputOutputArray& U = run.U;

  if ( rank == 0 || rank == np - 1) {

    if ( rank == 0 ) {
      for (int k   = 0; k< n1n2; ++k) { U[k][0][0]      = zero; }
      for (int k   = 0; k< n1n2; ++k) { U[k][n3+1][1]   = zero; }  /* ~ atop is zero on rank 0                   ~ */

      if ( iu3  > 2) {
        for (int k = 0; k< n1n2; ++k) { U[k][0][2]      = zero; }
        for (int k = 0; k< n1n2; ++k) { U[k][n3+1][3]   = zero; }  /* ~ vztop is zero on rank 0                  ~ */

      }
    }

    if ( rank == np - 1 ) {
      for (int k   = 0; k< n1n2; ++k) { U[k][n3][0]     = zero; }
      for (int k   = 0; k< n1n2; ++k) { U[k][n3+1][1]   = zero; }  /* ~ atop is zero on rank 0                   ~ */

      if ( iu3  > 2) {
        for (int k = 0; k< n1n2; ++k) { U[k][n3][2]     = zero; }
        for (int k = 0; k< n1n2; ++k) { U[k][n3+1][3]   = zero; }  /* ~ vztop is zero on rank 0                  ~ */

      }
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::retrieveOJ( stack& run )  {

  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c    ); /* ~ number of complex elements per layer       ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers ); /* ~ number of layers in stack                  ~ */

  ComplexArray::size_type usize;
  usize = n1n2c *n_layers;
  fftw.fftwForwardAll( run, O, J); /* ~ resulting O and J should compare well with OfromP and HfromA output below ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::OfromP( stack& run )  {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD,     &rank     ); /* ~ possibly needed for debugging                   ~ */
  int srun;     run.palette.fetch(   "srun",  &srun     ); /* ~ for srun > 1 we are restarting.                 ~ */
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c    ); /* ~ number of complex elements per layer            ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers ); /* ~ number of layers in stack                       ~ */
  int np;       run.palette.fetch(    "np"  , &np       ); /* ~ number of processes ~ */
  std::string init; run.palette.fetch( "initMode", &init);
  
  RealArray&    k2 = run.k2;                               /* ~ square magnitudes of k-space vectors            ~ */
  ComplexArray& U0 = run.U0;                               /* ~ holds phi (i.e. P ) at this point               ~ */ 
  
  ComplexArray::size_type usize;
  usize            = U0.capacity();                        /* ~ current capacity of U0 - should be known        ~ */
  
  assert(usize     == (n1n2c * n_layers));                 /* ~ test usize                                      ~ */
  assert(P.size()  == usize);
  assert(O.size()  == usize);
  
  int kstart;
  int kstop;

  if ( rank == 0 )          { kstart = n1n2c;               /* ~ layer 0 of P and O handled by BC's for process 0     ~ */
                              kstop  = usize;
  }
  else if (rank == (np-1))  { kstart = 0;                   /* ~ layer n3 of P and O handled by BC's for process np-1 ~ */
                              kstop  = (n_layers-2)*n1n2c;
  } 
  else {                      kstart = 0;                   /* ~ interior layers/processes                            ~ */
                              kstop  = usize;
  }

  for (unsigned k = kstart; k <kstop; k++) { P[k] = U0[k];} /* ~ copy the stream function and place it in P           ~ */
     
  if (!(init.compare("from_data") == 0 )) {
    unsigned  idx    = 0;                                     /* ~ index for k2                                         ~ */
    for (unsigned k  = 0; k < usize; k++) {

      if (k % n1n2c  == 0 ) { idx = 0; }                     /* ~ reset idx when starting new layer                     ~ */
    
      U0[k] = k2[idx] * P[k];                                /* ~ Omega = - delperp^2 P                                 ~ */
       O[k] = U0[k];
      ++idx;
    }

  }
  else {
  for (unsigned k = kstart; k <kstop; k++) {U0[k] = O[k];}  /* ~ O already initialized just needs to be placed in U0  ~ */
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::HfromA( stack& run )  {

  int rank;     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int srun;     run.palette.fetch(   "srun",  &srun     );
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c    ); /* ~ number of complex elements per layer       ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers ); /* ~ number of layers in stack                  ~ */
  std::string init; run.palette.fetch( "initMode", &init);
  
  RealArray&    k2 = run.k2;                               /* ~ square magnitude of k-space vectors        ~ */
  ComplexArray& U1 = run.U1;                               /* ~ holds A (flux function) at this point      ~ */

  ComplexArray::size_type usize;
  usize            = U1.capacity();                        /* ~ current capacity of U1 - should be known   ~ */
  assert(usize     == (n1n2c * n_layers));                 /* ~ test usize                                 ~ */
  assert(A.size()  == usize);
  assert(J.size()  == usize);

  RealVar ssqd;      run.palette.fetch(  "ssqd",  &ssqd ); /* ~ parameter sigma^2 needed for H             ~ */
  std::string model; run.palette.fetch(  "model", &model); 

  int kstart       = n1n2c;
  int kstop        = n1n2c*(n_layers - 1);

  if (!(init.compare("from_data") == 0)) {
    int idx          = 0;
    for (unsigned k  = kstart; k < kstop; k++) { 
      if (k % n1n2c == 0) {idx = 0;}
      A[k]           = U1[k]; 
      J[k]           = k2[idx]*A[k];
      if (model.compare("hall") == 0 ) {
        U1[k]        = A[k] + ssqd * J[k];                   /* ~ H = A + sigma^2 J                          ~ */
      }

      ++idx;
    }                                                         /* ~ preserve flux function in A                ~ */
//  std::cout << "this should not print" << std::endl;
  }
  else {
    for (unsigned k  = kstart; k < kstop; k++) { 
      A[k]           = U1[k]; 
    }                                                         /* ~ preserve flux function in A                ~ */
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::PfromO( stack& run )  {
 
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   );    /* ~ number of complex elements per layer        ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers);    /* ~ number of layers in stack                   ~ */

  RealArray&    k2     = run.k2;                             /* ~ square magnitude of k-space vectors         ~ */
  RealArray&    inv_k2 = run.inv_k2;                         /* ~ inverse of square magnitude of k-space vectors~ */
  ComplexArray& U0     = run.U0;                             /* ~ holds Omega (vorticity) at this point       ~ */

  ComplexArray::size_type usize;
  usize                = U0.capacity();                      /* ~ current capacity of U0 - should be known    ~ */

  assert(usize         == (n1n2c * n_layers));               /* ~ test usize                                  ~ */

  unsigned idx         = 0;                                  /* ~ index for inv_k2                            ~ */

  for (unsigned k = 0; k < usize; k++) {

    if ( k % n1n2c == 0) { idx = 0; }                        /* ~ reset idx when starting new layer           ~ */

    if ( k2[idx] != zero) { P[k]  = U0[k] / k2[idx];}        /* ~ Omega = - delperp^2 P                       ~ */
    else                  { U0[k] = czero;          }

    O[k]               = U0[k];                              /* ~ preserve vorticity for output               ~ */
    U0[k]              = P[k];

    ++idx;

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::AfromH( stack& run )  {

  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   ); /* ~ number of complex elements per layer       ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers); /* ~ number of layers in stack                  ~ */

  RealArray& k2    = run.k2;                              /* ~ square magnitude of k-space vectors        ~ */
  ComplexArray& U1 = run.U1;                              /* ~ holds H function at this point             ~ */

  ComplexArray::size_type usize;
  usize            = U1.capacity();                       /* ~ current capacity of U1 - should be known   ~ */

  assert(usize     == (n1n2c * n_layers));                /* ~ test usize                                 ~ */

  RealVar ssqd;      run.palette.fetch("ssqd",   &ssqd ); /* ~ parameter sigma^2 needed for A             ~ */
  std::string model; run.palette.fetch( "model", &model);

  unsigned idx     = 0;                                   /* ~ index for k2                               ~ */
  for (unsigned k  = 0; k < usize; k++) {

    if (k % n1n2c  == 0) { idx = 0; }                     /* ~ reset idx when starting new layer          ~ */

    if (model.compare("hall") == 0 ) {
      A[k]         = U1[k] / (one + ssqd*k2[idx]);
      U1[k]        = A[k];
    }
    else if( model.compare("rmhd") == 0 || model.compare("inhm") == 0 ) {

      A[k]         = U1[k];

    }
      J[k]         = k2[idx] * A[k];                       /* ~ J = -delperp A                             ~ */
    ++idx;
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalElls(    stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string model   ; run.palette.fetch("model",    &model   );
  std::string nprofile; run.palette.fetch("nprofile", &nprofile);
  int p3;               run.palette.fetch("p3",       &p3      );

  Elln.assign(p3+1,zero);
  EllA.assign(p3+1,zero);
  EllB.assign(p3+1,zero);

  kpm.assign( p3+1,zero);
  kpp.assign( p3+1,zero);
  kmm.assign( p3+1,zero);
  kmp.assign( p3+1,zero);

  if (model.compare("inhm") == 0) {
    int i_profile;

    if      (nprofile.compare("flat")   == 0) { i_profile =  0; }
    else if (nprofile.compare("torus")  == 0) { i_profile =  1; }
    else if (nprofile.compare("cloop")  == 0) { i_profile =  2; }
    else                                      { i_profile = -1; }

    switch(i_profile) {

    case(0) : Elln.assign(p3+1, zero); 
              EllA.assign(p3+1, zero);
              EllB.assign(p3+1, zero);
              break;
    case(1) : EllB.assign(p3+1, zero);
              for ( int l = 0; l <  p3+1;  ++l) {

                EllA[l] = std::abs( dvalfdz[l] / (two  * valfven[l]) );
                Elln[l] = std::abs( dndz[l]    / (four * nofz[l]   ) );

              }
              EllB.assign(p3+1, zero);
              break;
    case(2) : Elln.assign(p3+1, zero);
              EllA.assign(p3+1, zero);
              EllB.assign(p3+1, zero);
              break;

    default :

      std::cout << "evalN: WARNING - the profile " << nprofile << " is not implemented. Assuming a flat profile." << std::endl;

      Elln.assign(p3+1, zero);
      EllA.assign(p3+1, zero);
      EllB.assign(p3+1, zero);

    }

    for ( int l = 0; l <  p3+1;  ++l) {

      kpm[l] =  EllB[l] + EllA[l] - Elln[l] ;
      kpp[l] =  EllB[l] + EllA[l] + Elln[l] ;

      kmm[l] =  EllB[l] - EllA[l] - Elln[l] ;
      kmp[l] =  EllB[l] - EllA[l] + Elln[l] ;

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalValf( stack& run ) {
  
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string model; run.palette.fetch("model",    &model   ); /* ~ specification for density profile    ~ */
  int         p3;    run.palette.fetch("p3",       &p3      ); /* ~ layers per stack not incl. ghosts    ~ */

  if (model.compare("rmhd") == 0) {

    valfven.assign(p3+1, one  );
    dvalfdz.assign(p3+1, zero );
    nofz.assign(   p3+1, one  );
    dndz.assign(   p3+1, zero );

  }

  else {

    std::string nprofile; run.palette.fetch("nprofile", &nprofile); /* ~ specification for density profile    ~ */
    double zl;            run.palette.fetch("zl",       &zl      ); /* ~ size of domain along z               ~ */
    double ninf;          run.palette.fetch("ninf",     &ninf    ); /* ~ density at z = +/- infinity          ~ */
    double n0;            run.palette.fetch("n0",       &n0      ); /* ~ density at z_0 (i.e. z = zl / 2)     ~ */
    double valfmax;       run.palette.fetch("valfmax",  &valfmax );
    double valfmin;
    double H0;                                                      /* ~ density scale-height constant        ~ */

    RealArray& z = run.z;

    int    i_profile;
    if      (nprofile.compare("flat")   == 0) { i_profile =  0; }   /* ~ switch statements don't like strings ~ */
    else if (nprofile.compare("torus")  == 0) { i_profile =  1; }
    else if (nprofile.compare("cloop")  == 0) { i_profile =  2; }
    else                                      { i_profile = -1; }

    switch(i_profile) {

    case(0) : valfven.assign(p3+1, one );
              dvalfdz.assign(p3+1, zero);
                 nofz.assign(p3+1, one );
                 dndz.assign(p3+1, zero);
              break;
    case(1) : H0            = zl / sqrt( - eight * log( (one - ninf) / (n0 - ninf) ));
              valfmin       = valfmax / sqrt(n0);

              if (rank == 0){ std::cout << "evalValf: H0      = " << H0      << std::endl; }

                 nofz.assign(p3+1,one );
                 dndz.assign(p3+1,zero);
              valfven.assign(p3+1,one );
              dvalfdz.assign(p3+1,zero);

              if (rank == 0){ std::cout << "evalValf: test one "       << std::endl; }

              for (unsigned k = 0; k < p3+1; ++k) {

                 nofz[k]    = ninf + ((n0 - ninf) * (exp(-half*pow((z[k] - (half*zl))/H0 ,two))));
                 dndz[k]    = -(n0 - ninf) * exp( -half * pow(((z[k]-(half*zl))/H0),2)) * (z[k] - (half*zl)) / ( pow(H0,2) );
                 valfven[k] = valfmin * sqrt( n0 / nofz[k] );
                 dvalfdz[k] = -half * valfven[k] * dndz[k] / nofz[k];

              }

              if (rank == 0){ std::cout << "evalValf: test two "       << std::endl; }

              break;
    case(2) : valfven.assign(p3+1, one ); 
              dvalfdz.assign(p3+1, zero);
                 nofz.assign(p3+1, one );
                 dndz.assign(p3+1, zero);
              break;

    default : 

      std::cout << "evalValf: WARNING - the profile " << nprofile << " is not implemented. Assuming a flat profile." << std::endl;
      valfven.assign(p3+1, one );
      dvalfdz.assign(p3+1, zero);

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalUmean( stack& run ) {
  
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string model; run.palette.fetch("model",    &model   ); /* ~ specification for density profile    ~ */
  int         p3;    run.palette.fetch("p3",       &p3      );

  umean.assign(p3+1, zero); 
  dudz.assign(p3+1,  zero);

  if (model.compare("rmhd") != 0)  {

    std::string uprofile; run.palette.fetch("uprofile", &uprofile);

    int        i_profile;

    if      (uprofile.compare("noflow")   == 0) { i_profile =  0; }
    else if (uprofile.compare("torus")    == 0) { i_profile =  1; }
    else if (uprofile.compare("whoknows") == 0) { i_profile =  2; }
    else                                        { i_profile = -1; }

    switch(i_profile) {

    case(0) : umean.assign(p3+1, zero); 
              dudz.assign(p3+1,  zero);
              break;
    case(1) : umean.assign(p3+1, zero);
              for (unsigned k = 0; k < p3 + 1; ++k){umean[k] = one / nofz[k];}
              for (unsigned k = 0; k < p3 + 1; ++k){dudz[k]  = (-one / (nofz[k]*nofz[k]))*(dndz[k]);}
              break;
    case(2) : umean.assign(p3+1, zero);
              dudz.assign(p3+1,  zero); break;

    default : 

      std::cout << "evalValf: WARNING - the profile " << uprofile << " is not implemented. Assuming a flat profile." << std::endl;

      umean.assign(p3+1,zero);
      dudz.assign(p3+1, zero);

    }

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::trackEnergies( RealVar t_cur, stack& run ) {


  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  RealArray& EnergyQs    = run.EnergyQs;

  static const int i_tcr =  0; double tcr;  /* ~ Current time                                     ~ */
  static const int i_pe  =  1; double pe ;  /* ~ Total perpendicular kinetic energy               ~ */
  static const int i_ae  =  2; double ae ;  /* ~ Total perpendicular perturbed magnetic energy    ~ */
  static const int i_mo  =  3; double mo ;  /* ~ Maximum vorticity                                ~ */
  static const int i_imo =  4; double imo;  /* ~ i index location of maximum vorticity            ~ */
  static const int i_jmo =  5; double jmo;  /* ~ j index location of maximum vorticity            ~ */
  static const int i_oe  =  6; double oe ;  /* ~ Total square magnitude of vorticity              ~ */
  static const int i_mj  =  7; double mj ;  /* ~ Maximum current                                  ~ */
  static const int i_imj =  8; double imj;  /* ~ i index location of maximum current              ~ */
  static const int i_jmj =  9; double jmj;  /* ~ j index location of maximum current              ~ */
  static const int i_ce  = 10; double ce ;  /* ~ Total square magnitude of current                ~ */
  static const int i_noe = 11; double noe;  /* ~ viscous dissipation                              ~ */
  static const int i_ece = 12; double ece;  /* ~ resistive dissipation                            ~ */
  static const int i_fe  = 13; double fe ;  /* ~ Poynting Flux                                    ~ */
  static const int i_ftp = 14; double ftp;  /* ~ Time-averaged Poynting Flux                      ~ */
  static const int i_eds = 15; double eds;  /* ~ Time average of total dissipated energy          ~ */
  static const int i_dng = 16; double dng;  /* ~ Average rate-of-change of total energy           ~ */
  static const int i_irc = 17; double irc;  /* ~ Instantaneous rate-of-change of total energy     ~ */
  static const int i_gml = 18; double gml;  /* ~ time-average of energy gains minus losses        ~ */
  static const int i_cns = 19; double cns;  /* ~ Check on energy conservation                     ~ */
  static const int i_ttc = 20; double ttc;  /* ~ Current time in units of correlation time        ~ */
  static const int i_dt  = 21; double dt ;  /* ~ current value of time increment                  ~ */
  static const int i_dtv = 22; double dtv;  /* ~ time increment adjustment parameter              ~ */
  static const int i_coc = 23; double coc;  /* ~ root mean square J^2 / Del J^2                   ~ */
  static const int i_vkt = 24; double vkt;  /* ~ average rate-of-change of j^2 / magnetic energy  ~ */ 
  static const int i_avm = 25; double avm;  /* ~ Time averaged magnetic field strength            ~ */
  static const int i_avp = 26; double avp;  /* ~ Time-averaged footpoint velocity                 ~ */
  static const int i_fp  = 27; double fp ;  /* ~ Total footpoint kinetic energy at z = 0          ~ */                  
  static const int i_he  = 28; double he ;  /* ~ Total internal energy due to inhomogeneities     ~ */

  RealVar dtvb;   physics_data.fetch("dtvb",   &dtvb   );
  RealVar AVEz;   run.palette.fetch("AVEz",    &AVEz   );
  RealVar AVEpv;  run.palette.fetch("AVEpv",   &AVEpv  );
  RealVar AVEpn;  run.palette.fetch("AVEpn",   &AVEpn  );
  RealVar AVEpe;  run.palette.fetch("AVEpe",   &AVEpe  );
  RealVar tauC;   run.palette.fetch( "tauC",   &tauC   );
  RealVar eta;    run.palette.fetch( "eta",    &eta    );
  RealVar nu;     run.palette.fetch( "nu",     &nu     );
  RealVar tstart; run.palette.fetch( "tstart", &tstart );

  int     srun;   run.palette.fetch("srun",    &srun   );
  int     ndt;    run.palette.fetch("ndt",     &ndt    );

  RealVar cee;
  RealVar t_old;
  RealVar dt_old;
  RealVar aeold;
  RealVar peold;
  RealVar heold;

  pe           = evalTotalKineticEnergy(  run );
  ae           = evalTotalMagneticEnergy( run );
  oe           = evalTotalVorticitySqd(   run );
  ce           = evalTotalCurrentSqd(     run );
  cee          = evalTotalGradCurrentSqd( run );
  fp           = evalTotalFootPointKE(    run );
  fe           = evalTotalPoyntingFlux(   run );

  std::string model; run.palette.fetch("model", &model);

  if (     model.compare("rmhd") == 0) { he = zero;                           }
  else if (model.compare("inhm") == 0) { he = evalTotalHelicalEnergy(  run ); }

  if ( EnergyQs.size() == 0) {

    if ( srun == 1 ) {
      if (t_cur == zero) {

        EnergyQs.assign(i_he + 1,zero);

        t_old  = t_cur;
        run.palette.fetch("dt", &dt_old);
        aeold  = ae;
        peold  = pe;
        heold  = he;
        noe    = zero;
        ece    = zero;
        ftp    = zero;
        eds    = zero;
        dng    = zero;
        vkt    = zero;
        avm    = zero;
        avp    = zero;

      }
      else {std::cout << "trackEnergies: WARNING - EnergyQs is unallocated for tstart = zero during first sub-run." << std::endl;}
    }
    else { /* ~ there will be a file ~ */

      if (rank == 0) {

        std::string prefix;    run.palette.fetch(   "prefix",    &prefix   );
        std::string run_label; run.palette.fetch(   "run_label", &run_label);
        std::string res_str;   run.stack_data.fetch("res_str",   &res_str  );

        std::string data_dir;  run.palette.fetch(   "data_dir",  &data_dir );

        std::string energy_data_file    = "./" + data_dir + "/" + prefix + "_" + res_str + ".o" + run_label;
        const char *c_data_file         = energy_data_file.c_str();

        std::ifstream ifs;
        ifs.open( c_data_file, std::ios::in );

        if (ifs.good()) {

          unsigned esize                = i_he + 1;

          std::streampos begin, end;
          std::streamoff bytes_per_line = (esize*24 + 28);
          RealVar nextE;
          begin                         = ifs.tellg();
          ifs.seekg(-bytes_per_line, std::ios::end);

          for (unsigned k = 0; k < esize; k++) {
            ifs >> nextE;
            EnergyQs.push_back(nextE);
          }
          ifs.close();

        }   // ifs is good
        else { std::cout << "trackEnergyQs: Warning - could not open file " << energy_data_file << std::endl; }

        t_old                           = EnergyQs[i_tcr];
        dt_old                          = EnergyQs[i_dt ];

        aeold                           = EnergyQs[i_ae ];
        peold                           = EnergyQs[i_pe ];
        heold                           = EnergyQs[i_he ];

        noe                             = nu  * oe;
        ece                             = eta * ce;

        ftp                             = ((EnergyQs[i_ftp]*t_old) + fe          * dt_old)        / t_cur;
        eds                             = ((EnergyQs[i_eds]*t_old) + (noe + ece) * dt_old)        / t_cur;
        dng                             = ( (EnergyQs[i_dng] * t_old) + (ae - aeold + pe - peold + he - heold) ) / t_cur;
        avm                             = ( EnergyQs[i_avm] * EnergyQs[i_avm] ) * t_old;
        avm                             = sqrt( (avm + ae * dt_old)        / t_cur );
        avp                             = sqrt(two * (AVEpe + fp * dt_old) / t_cur );

      }     // rank is zero
    }       // at beginning of srun > 1
  }         // EnergyQs not yet allocated
  else {    // EnergyQs is already allocated

    t_old                               = EnergyQs[i_tcr];
    dt_old                              = EnergyQs[i_dt ];
    aeold                               = EnergyQs[i_ae ];
    peold                               = EnergyQs[i_pe ];
    heold                               = EnergyQs[i_he ];

    noe                                 = nu  * oe;
    ece                                 = eta * ce;

    ftp                                 = ((EnergyQs[i_ftp]*t_old) + fe          * dt_old)                       / t_cur;
    eds                                 = ((EnergyQs[i_eds]*t_old) + (noe + ece) * dt_old)                       / t_cur;
    dng                                 = ( (EnergyQs[i_dng] * t_old) + (ae - aeold + pe - peold + he - heold) ) / t_cur;

    avm                                 = ( EnergyQs[i_avm] * EnergyQs[i_avm] ) * t_old;
    avm                                 = sqrt( (avm + ae * dt_old)        / t_cur );
    avp                                 = sqrt(two * (AVEpe + fp * dt_old) / t_cur );

  }         // EnergyQs is already allocated

/* ~ now calculate all the stuff that relies on above evualations ~ */

  if (rank == 0 ) {

    cns                                 = std::abs(((fe - noe - ece )-((ae - aeold + pe - peold + he - heold)/dt_old))*dt_old);

    std::cout << " " <<std::endl;
    std::cout << "trackEn: cns     =  " << cns     << std::endl;
    std::cout << "trackEn: fe      =  " << fe      << std::endl;
    std::cout << "trackEn: noe     =  " << noe     << std::endl;
    std::cout << "trackEn: ece     =  " << ece     << std::endl;
    std::cout << "trackEn: ae      =  " << ae      << std::endl;
    std::cout << "trackEn: aeold   =  " << aeold   << std::endl;
    std::cout << "trackEn: pe      =  " << pe      << std::endl;
    std::cout << "trackEn: peold   =  " << peold   << std::endl;
    std::cout << "trackEn: he      =  " << he      << std::endl;
    std::cout << "trackEn: heold   =  " << heold   << std::endl;
    std::cout << "trackEn: dt_old  =  " << dt_old  << std::endl;
    std::cout << "trackEn: t_cur   =  " << t_cur   << std::endl;
    std::cout << " " <<std::endl;

    mo                                  = zero;
    imo                                 = zero;
    jmo                                 = zero;

    mj                                  = zero;
    imj                                 = zero;
    jmj                                 = zero;

    AVEz                                = AVEz  + cns * dt_old;
    AVEpv                               = AVEpv + noe * dt_old;
    AVEpn                               = AVEpn + ece * dt_old;
    AVEpe                               = AVEpe + fp  * dt_old;

    run.palette.reset("AVEz",  AVEz  );
    run.palette.reset("AVEpv", AVEpv );
    run.palette.reset("AVEpn", AVEpn );
    run.palette.reset("AVEpe", AVEpe );

    ttc                                 = t_cur / tauC;
    run.palette.fetch("dt", &dt);
    dtv                                 = dtvb;
    if (cee == zero){ coc               = zero;         }
    else            { coc               = sqrt(ce/cee); }
    if (ae != zero ) {
      vkt                               = EnergyQs[i_vkt] * t_old;
      vkt                               = (vkt + (ce/ae)  * dt_old) / t_cur;
    }
    else { vkt                          = zero;}


    RealVar AVEz;  run.palette.fetch("AVEz",  &AVEz );
    RealVar AVEpv; run.palette.fetch("AVEpv", &AVEpv);
    RealVar AVEpn; run.palette.fetch("AVEpn", &AVEpn);
    RealVar AVEpe; run.palette.fetch("AVEpe", &AVEpe);

    irc                                 = (ae - aeold + pe - peold + he - heold) / dt_old;
    gml                                 = ftp - eds;

    AVEz                                = AVEz  + cns * dt_old;
    AVEpv                               = AVEpv + noe * dt_old;
    AVEpn                               = AVEpn + ece * dt_old;
    AVEpe                               = AVEpe + fp  * dt_old;

    run.palette.reset("AVEz",  AVEz );

    EnergyQs[ i_tcr ]                   = t_cur;
    EnergyQs[ i_pe  ]                   = pe ;
    EnergyQs[ i_ae  ]                   = ae ;
    EnergyQs[ i_mo  ]                   = mo ;
    EnergyQs[ i_imo ]                   = imo;
    EnergyQs[ i_jmo ]                   = jmo;
    EnergyQs[ i_oe  ]                   = oe ;
    EnergyQs[ i_mj  ]                   = mj ;
    EnergyQs[ i_imj ]                   = imj;
    EnergyQs[ i_jmj ]                   = jmj;
    EnergyQs[ i_ce  ]                   = ce ;
    EnergyQs[ i_noe ]                   = noe;      
    EnergyQs[ i_ece ]                   = ece;
    EnergyQs[ i_fe  ]                   = fe ;
    EnergyQs[ i_ftp ]                   = ftp;
    EnergyQs[ i_eds ]                   = eds;
    EnergyQs[ i_dng ]                   = dng;
    EnergyQs[ i_irc ]                   = irc;
    EnergyQs[ i_gml ]                   = gml;
    EnergyQs[ i_cns ]                   = cns;
    EnergyQs[ i_ttc ]                   = ttc;
    EnergyQs[ i_dt  ]                   = dt ;
    EnergyQs[ i_dtv ]                   = dtv;
    EnergyQs[ i_coc ]                   = coc;
    EnergyQs[ i_vkt ]                   = vkt;
    EnergyQs[ i_avm ]                   = avm;
    EnergyQs[ i_avp ]                   = avp;
    EnergyQs[ i_fp  ]                   =  fp;
    EnergyQs[ i_he  ]                   =  he;

  }
}          // end function

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::trackPowerSpectra( RealVar t_cur, stack& run ) {

  int calcsvz; run.palette.fetch("calcsvz", &calcsvz);
  if (calcsvz == 1) {

     int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     int   n3; run.palette.fetch("p3",   &n3    );

     static const int i_k     =  0;
     static const int i_spe   =  1;
     static const int i_sae   =  2;
     static const int i_ts    =  3;
     static const int i_szp   =  4;
     static const int i_szm   =  5;
     static const int i_tz    =  6;
     static const int i_pf    =  7;
     static const int i_af    =  8;
     static const int i_zpf   =  9;
     static const int i_zmf   = 10;

     int n1;    run.stack_data.fetch("n1",    &n1);
     int n2;    run.stack_data.fetch("n2",    &n2);
     int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c);

     RealArray& k2            = run.k2;
     int ksize                = k2.capacity();

     RealArray spe;
     RealArray sae;
     RealArray szp;
     RealArray szm;

//   RealArray pef;
//   RealArray aef;
//   RealArray zpf;
//   RealArray zmf;

     assert(n1n2c == ksize);

     int isp;physics_data.fetch("isp", &isp);

     for (unsigned l = 1; l < n3 + 1; ++l ) {

       spe.assign(isp+2,zero);
       sae.assign(isp+2,zero);
       szp.assign(isp+2,zero);
       szm.assign(isp+2,zero);

//     pef.assign(isp+1,zero);
//     aef.assign(isp+1,zero);
//     zpf.assign(isp+1,zero);
//     zmf.assign(isp+1,zero);

       

       int m;
       int idx;

       for (unsigned k = 0; k < n1n2c; ++k) {

         if (sqrt(k2[k]) >= kb && sqrt(k2[k]) <= kf) {
           m       = sqrt(k2[k]) * dk_m1;
           idx     = (l*n1n2c) + k;

           if ( m <= isp+1 ) {

             spe[m]  = spe[m] + k2[k] * ( pow(std::norm(P[idx]),         2) );
             sae[m]  = sae[m] + k2[k] * ( pow(std::norm(A[idx]),         2) );
             szp[m]  = szp[m] + k2[k] * ( pow(std::norm(P[idx] + A[idx]),2) );
             szm[m]  = szm[m] + k2[k] * ( pow(std::norm(P[idx] - A[idx]),2) );

           }

           else { std::cout << "trackPowerSpectra: WARNING - index m = " << m << " this exceeds upper limit isp + 1 = " << isp + 1 << std::endl;}        

         }
       }
  
         for (unsigned j = 0; j < isp + 1; ++j) {
  
           spe[j] = two * spe[j] * dk_m1;
           sae[j] = two * sae[j] * dk_m1;
           szp[j] = two * szp[j] * dk_m1;
           szm[j] = two * szm[j] * dk_m1;
  
           SpcVsZ[j][l-1][i_k  ] = j * dk;
           SpcVsZ[j][l-1][i_spe] = SpcVsZ[j][l-1][i_spe] + spe[j];
           SpcVsZ[j][l-1][i_sae] = SpcVsZ[j][l-1][i_sae] + sae[j];
           SpcVsZ[j][l-1][i_ts ] = SpcVsZ[j][l-1][i_ts ] + spe[j] + sae[j];
           SpcVsZ[j][l-1][i_szp] = SpcVsZ[j][l-1][i_szp] + szp[j];
           SpcVsZ[j][l-1][i_szm] = SpcVsZ[j][l-1][i_szm] + szm[j];
           SpcVsZ[j][l-1][i_tz ] = SpcVsZ[j][l-1][i_tz ] + szp[j] + szm[j];
  
         }

         for (unsigned j = 0; j < nk - 1;++j) {

           if (spe[ikb+j] == zero) { SpcVsZ[j][l-1][i_pf ] = SpcVsZ[j][l-1][i_pf ] + ((RealVar) (-300.0)   ); }
           else                    { SpcVsZ[j][l-1][i_pf ] = SpcVsZ[j][l-1][i_pf ] + std::log10(spe[ikb+j] ); }
           if (sae[ikb+j] == zero) { SpcVsZ[j][l-1][i_af ] = SpcVsZ[j][l-1][i_af ] + ((RealVar) (-300.0)   ); }
           else                    { SpcVsZ[j][l-1][i_af ] = SpcVsZ[j][l-1][i_af ] + std::log10(sae[ikb+j] ); }
           if (szp[ikb+j] == zero) { SpcVsZ[j][l-1][i_zpf] = SpcVsZ[j][l-1][i_zpf] + ((RealVar) (-300.0)   ); }
           else                    { SpcVsZ[j][l-1][i_zpf] = SpcVsZ[j][l-1][i_zpf] + std::log10(szp[ikb+j] ); }
           if (szm[ikb+j] == zero) { SpcVsZ[j][l-1][i_zmf] = SpcVsZ[j][l-1][i_zmf] + ((RealVar) (-300.0)   ); }
           else                    { SpcVsZ[j][l-1][i_zmf] = SpcVsZ[j][l-1][i_zmf] + std::log10(szm[ikb+j] ); }

         }
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::trackQtyVsZ(RealVar t_cur, stack& run ) {

  int calcqvz; run.palette.fetch("calcqvz", &calcqvz);

  if (calcqvz == 1) {

    RealVar eta;  run.palette.fetch(    "eta",     &eta );
    RealVar nu;   run.palette.fetch(    "nu",      &nu  );
    RealVar dz;   run.stack_data.fetch( "dz",      &dz  );
    int     n3;   run.palette.fetch(    "p3",      &n3  );
    int     np;   run.palette.fetch(    "np",      &np  );
    int     rank; MPI_Comm_rank(MPI_COMM_WORLD,    &rank);

    static const int i_z    =  0;
    static const int i_ae   =  1;
    static const int i_pe   =  2;
    static const int i_ch   =  3;
    static const int i_ep   =  4;
    static const int i_em   =  5;
    static const int i_ce   =  6;
    static const int i_oe   =  7;
    static const int i_zp   =  8;
    static const int i_zm   =  9;
    static const int i_nch  = 10;
    static const int i_lc   = 11;
    static const int i_lo   = 12;
    static const int i_lzp  = 13;
    static const int i_lzm  = 14;
    static const int i_tez  = 15;
    static const int i_rez  = 16;
    static const int i_nez  = 17;
    static const int i_rzp  = 18;
    static const int i_rzm  = 19;
    static const int i_zpm  = 20;
    static const int i_lpl  = 21;
    static const int i_lmi  = 22;
    static const int i_lpla = 23;
    static const int i_lmia = 24;
    static const int i_cnse = 25;

    RealArray& z            = run.z;

    RealArray aevsz;    aevsz.assign(   n3,zero); /* ~ magnetic energy by layer                                             ~ */
    RealArray pevsz;    pevsz.assign(   n3,zero); /* ~ kinetic  energy by layer                                             ~ */
    RealArray tevsz;    tevsz.assign(   n3,zero); /* ~ total energy by layer                                                ~ */
    RealArray revsz;    revsz.assign(   n3,zero); /* ~ residual energy by layer                                             ~ */
    RealArray nevsz;    nevsz.assign(   n3,zero); /* ~ normalized residual energy by layer                                  ~ */
    RealArray chvsz;    chvsz.assign(   n3,zero); /* ~ cross helicity by layer                                              ~ */
    RealArray cevsz;    cevsz.assign(   n3,zero); /* ~ Ohmic dissipation per unit magnetic diffusivity by layer             ~ */
    RealArray oevsz;    oevsz.assign(   n3,zero); /* ~ viscous dissipation per unit kinematic viscosity by layer            ~ */
    RealArray zpvsz;    zpvsz.assign(   n3,zero); /* ~ Z^+ elsasser dissipation per unit zeta by layer                      ~ */
    RealArray zmvsz;    zmvsz.assign(   n3,zero); /* ~ Z^- elsasser dissipation per unit zeta by layer                      ~ */
    RealArray zepvsz;   zepvsz.assign(  n3,zero); /* ~ Z^+ elsasser energy by layer                                         ~ */
    RealArray zemvsz;   zemvsz.assign(  n3,zero); /* ~ Z^- elsasser energy by layer                                         ~ */
    RealArray nchvsz;   nchvsz.assign(  n3,zero); /* ~ normalized cross helicity by layer                                   ~ */
    RealArray lcvsz;    lcvsz.assign(   n3,zero); /* ~ current density based scale length by layer                          ~ */
    RealArray lovsz;    lovsz.assign(   n3,zero); /* ~ vorticity based scale length by layer                                ~ */
    RealArray rzpvsz;   rzpvsz.assign(  n3,zero); /* ~ square root of Z^+ elsasser energy by layer                          ~ */
    RealArray rzmvsz;   rzmvsz.assign(  n3,zero); /* ~ square root of Z^- elsasser energy by layer                          ~ */
    RealArray zpmvsz;   zpmvsz.assign(  n3,zero); /* ~ product of rzpvsz and rzmvsz                                         ~ */
    RealArray lzpvsz;   lzpvsz.assign(  n3,zero); /* ~ Z^+ Elsasser vorticity length scale by layer                         ~ */
    RealArray lzmvsz;   lzmvsz.assign(  n3,zero); /* ~ Z^- Elsasser vorticity length scale by layer                         ~ */
    RealArray lplvsz;   lplvsz.assign(  n3,zero); /* ~ Z^+ Elsasser length scale from F.T. of Z^+                           ~ */
    RealArray lmivsz;   lmivsz.assign(  n3,zero); /* ~ Z^- Elsasser length scale from F.T. of Z^-                           ~ */
    RealArray lplavsz;  lplavsz.assign( n3,zero); /* ~ Z^+ Elsasser length scale from F.T. of Z^+ (alternative definition)  ~ */
    RealArray lmiavsz;  lmiavsz.assign( n3,zero); /* ~ Z^- Elsasser length scale from F.T. of Z^- (alternative definition)  ~ */
    RealArray consevsz; consevsz.assign(n3,zero); /* ~ Energy conservation check y layer                                    ~ */

    int kdx;
    RealArray& k2  = run.k2;
    int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );

    int capa       = A.size();

    assert(capa == (n3+2)*n1n2c);

    ComplexArray ZP; ZP.assign( capa, czero);
    ComplexArray ZM; ZM.assign( capa ,czero);

    for (unsigned k = 0; k < ((n3+2)*n1n2c); k++) {
  
      ZP[k] = A[k] + P[k];
      ZM[k] = A[k] - P[k];
  
    }

    for (unsigned l = 0; l < n3; l++) {

      kdx          = 0;

      for ( unsigned k = ((l+1) * n1n2c) ; k < ((l+2) * n1n2c); k++ ) {

        aevsz[l]   = (aevsz[l] + (dz*k2[kdx]           * ( ( A[k].real() * A[k].real())  + ( A[k].imag() *  A[k].imag()))  ));
        pevsz[l]   = (pevsz[l] + (dz*k2[kdx]           * ( ( P[k].real() * P[k].real())  + ( P[k].imag() *  P[k].imag()))  ));

        chvsz[l]   = (chvsz[l] + (dz*k2[kdx] * two     * ( ( P[k].real() * A[k].real())  + ( P[k].imag() *  A[k].imag()))  ));

        cevsz[l]   = (cevsz[l] + (dz*k2[kdx] * k2[kdx] * ( ( A[k].real() * A[k].real())  + ( A[k].imag() *  A[k].imag()))  ));
        oevsz[l]   = (oevsz[l] + (dz*k2[kdx] * k2[kdx] * ( ( P[k].real() * P[k].real())  + ( P[k].imag() *  P[k].imag()))  ));

        zpvsz[l]   = (zpvsz[l] + (dz*k2[kdx] * k2[kdx] * ( (ZP[k].real() * ZP[k].real()) + (ZP[k].imag() * ZP[k].imag())) ));
        zmvsz[l]   = (zmvsz[l] + (dz*k2[kdx] * k2[kdx] * ( (ZM[k].real() * ZM[k].real()) + (ZM[k].imag() * ZM[k].imag())) ));

        zepvsz[l]  = (zpvsz[l] + (dz*k2[kdx]           * ( (ZP[k].real() * ZP[k].real()) + (ZP[k].imag() * ZP[k].imag())) ));
        zemvsz[l]  = (zmvsz[l] + (dz*k2[kdx]           * ( (ZM[k].real() * ZM[k].real()) + (ZM[k].imag() * ZM[k].imag())) ));

        lplvsz[l]  = lplvsz[l]    + (sqrt(k2[kdx])     * ( (ZP[k].real() * ZP[k].real()) + (ZP[k].imag() * ZP[k].imag())) ); 
        lmivsz[l]  = lmivsz[l]    + (sqrt(k2[kdx])     * ( (ZM[k].real() * ZM[k].real()) + (ZM[k].imag() * ZM[k].imag())) ); 

        lplavsz[l] = lplvsz[l]    +                      ( (ZP[k].real() * ZP[k].real()) + (ZP[k].imag() * ZP[k].imag()))  ;
        lmiavsz[l] = lmivsz[l]    +                      ( (ZM[k].real() * ZM[k].real()) + (ZM[k].imag() * ZM[k].imag()))  ;

        ++kdx;

      }

      rzpvsz[l]    = sqrt(zepvsz[l]);
      rzmvsz[l]    = sqrt(zemvsz[l]);

      zpmvsz[l]    = rzpvsz[l]*rzmvsz[l];

      tevsz[l]     = pevsz[l] + aevsz[l];
      revsz[l]     = pevsz[l] - aevsz[l];

      if (l > 0) {

        consevsz[l] = ( tevsz[l] - QtyVsZ[(rank*n3) + l][i_tez] )                  \
                       + umean[l+1]*(  (tevsz[l] - tevsz[l-1])/ dz                 \
                                  - two*(EllA[l+1]*revsz[l] - Elln[l+1]*tevsz[l])  \
                                )                                                  \
                       - valfven[l+1] *(                                           \
                                          (chvsz[l] - chvsz[l-1]) / dz             \
                                         - two*Elln[l+1] * chvsz[l]                \
                                       )                                           \
                       +(nu * oevsz[l]) + (eta *cevsz[l] )                         \
                               ;
    }
    else {

        consevsz[l] = ( tevsz[l] - QtyVsZ[(rank*n3) + l][i_tez] )                  \
                       + umean[l+1]*(   tevsz[l] / dz                              \
                                  - two*(EllA[l+1]*revsz[l] - Elln[l+1]*tevsz[l])  \
                                )                                                  \
                       - valfven[l+1] *(                                           \
                                           chvsz[l]/ dz                            \
                                         - two*Elln[l+1] * chvsz[l]                \
                                       )                                           \
                       +(nu * oevsz[l]) + (eta *cevsz[l] )                         \
                                ;
    }

      if (std::abs(tevsz[l]) >=teensy) {
        nevsz[l] = revsz[l]/tevsz[l];
      }
      else { if      (revsz[l] > zero ) { nevsz[l] =  huge;}
             else if (revsz[l] < zero ) { nevsz[l] = -huge;}
             else                       { nevsz[l] =  zero;}
      }

      if (std::abs(zepvsz[l]) <= teensy ) { 
        lplvsz[l]  = huge;
        lplavsz[l] = huge;
      }
      else {
        lplvsz[l]  =      lplvsz[l] *dz / zepvsz[l];
        lplavsz[l] = sqrt(lplavsz[l]*dz / zepvsz[l]);
      }
      if (std::abs(zemvsz[l]) <= teensy ) {
        lmivsz[l]  = huge;
        lmiavsz[l] = huge;
      }
      else {
        lmivsz[l]  =      lmivsz[l] *dz / zemvsz[l];
        lmiavsz[l] = sqrt(lmiavsz[l]*dz / zemvsz[l]);
      }

/* ~ Normalized cross helicity ~ */

      if ( (std::abs(aevsz[l]) >= teensy) || (std::abs(pevsz[l])  >= teensy) ) {
        nchvsz[l] = chvsz[l] / (aevsz[l] + pevsz[l]) ;
      }
      else { nchvsz[l] = zero; }

/* ~ Current length scale ~ */

      if ( std::abs(cevsz[l]) >= teensy ) {
        lcvsz[l] = sqrt(std::abs(aevsz[l] / cevsz[l]));
        if ( std::abs(lcvsz[l]) > huge  ) { lcvsz[l] = huge; }
        if ( std::abs(lcvsz[l]) < teensy) { lcvsz[l] = zero; }
      }
      else{ lcvsz[l] = huge;}

/* ~ Vorticity length scale ~ */

      if ( std::abs(oevsz[l]) >= teensy) {
        lovsz[l]  = pevsz[l] / oevsz[l];
        if ( std::abs(lovsz[l]) > huge  ) { lovsz[l] = huge; }
        if ( std::abs(lovsz[l]) < teensy) { lovsz[l] = zero; }
      }
      else {lovsz[l] = huge ;}

/* ~ Elsasser z^+ length scale ~ */

      if (std::abs(zpvsz[l]) >= teensy) {
        lzpvsz[l] = sqrt(std::abs(zepvsz[l] / zpvsz[l]));
        if ( std::abs(lzpvsz[l]) > huge  ) { lzpvsz[l] = huge; }
        if ( std::abs(lzpvsz[l]) < teensy) { lzpvsz[l] = zero; }
      }
      else{lzpvsz[l] = huge;}

/* ~ Elsasser z^- length scale ~ */

      if (std::abs(zmvsz[l]) >= teensy) {
        lzmvsz[l] = zemvsz[l] / zmvsz[l];
        if ( std::abs(lzmvsz[l]) > huge  )  {lzmvsz[l] = huge; }
        if ( std::abs(lzmvsz[l]) < teensy)  {lzmvsz[l] = zero; }
      }
      else{ lzmvsz[l] = huge; }

      QtyVsZ[(rank*n3) + l][  i_z]   = z[l + 1];
      QtyVsZ[(rank*n3) + l][ i_ae]   = QtyVsZ[(rank*n3) + l][ i_ae]  + aevsz[ l];
      QtyVsZ[(rank*n3) + l][ i_pe]   = QtyVsZ[(rank*n3) + l][ i_pe]  + pevsz[ l];
      QtyVsZ[(rank*n3) + l][ i_ch]   = QtyVsZ[(rank*n3) + l][ i_ch]  + chvsz[ l];
      QtyVsZ[(rank*n3) + l][ i_ep]   = QtyVsZ[(rank*n3) + l][ i_ep]  + cevsz[ l];
      QtyVsZ[(rank*n3) + l][ i_em]   = QtyVsZ[(rank*n3) + l][ i_em]  + oevsz[ l];
      QtyVsZ[(rank*n3) + l][ i_ce]   = QtyVsZ[(rank*n3) + l][ i_ce]  + zpvsz[ l];
      QtyVsZ[(rank*n3) + l][ i_oe]   = QtyVsZ[(rank*n3) + l][ i_oe]  + zmvsz[ l];
      QtyVsZ[(rank*n3) + l][ i_zp]   = QtyVsZ[(rank*n3) + l][ i_zp]  + zepvsz[l];
      QtyVsZ[(rank*n3) + l][ i_zm]   = QtyVsZ[(rank*n3) + l][ i_zm]  + zemvsz[l];
      QtyVsZ[(rank*n3) + l][i_nch]   = QtyVsZ[(rank*n3) + l][i_nch]  + nchvsz[l];
      QtyVsZ[(rank*n3) + l][ i_lc]   = QtyVsZ[(rank*n3) + l][ i_lc]  + lcvsz[ l];
      QtyVsZ[(rank*n3) + l][ i_lo]   = QtyVsZ[(rank*n3) + l][ i_lo]  + lovsz[ l];
      QtyVsZ[(rank*n3) + l][i_lzp]   = QtyVsZ[(rank*n3) + l][i_lzp]  + lzpvsz[l];
      QtyVsZ[(rank*n3) + l][i_lzm]   = QtyVsZ[(rank*n3) + l][i_lzm]  + lzmvsz[l];
      QtyVsZ[(rank*n3) + l][i_tez]   = QtyVsZ[(rank*n3) + l][i_tez]  + tevsz[l];
      QtyVsZ[(rank*n3) + l][i_rez]   = QtyVsZ[(rank*n3) + l][i_rez]  + revsz[l];
      QtyVsZ[(rank*n3) + l][i_nez]   = QtyVsZ[(rank*n3) + l][i_nez]  + nevsz[l];
      QtyVsZ[(rank*n3) + l][i_rzp]   = QtyVsZ[(rank*n3) + l][i_rzp]  + rzpvsz[l];
      QtyVsZ[(rank*n3) + l][i_rzm]   = QtyVsZ[(rank*n3) + l][i_rzm]  + rzmvsz[l];
      QtyVsZ[(rank*n3) + l][i_zpm]   = QtyVsZ[(rank*n3) + l][i_zpm]  + zpmvsz[l];
      QtyVsZ[(rank*n3) + l][i_lpl]   = QtyVsZ[(rank*n3) + l][i_lpl]  + lplvsz[l];
      QtyVsZ[(rank*n3) + l][i_lmi]   = QtyVsZ[(rank*n3) + l][i_lmi]  + lmivsz[l];
      QtyVsZ[(rank*n3) + l][i_lpla]  = QtyVsZ[(rank*n3) + l][i_lpla] + lplavsz[l];
      QtyVsZ[(rank*n3) + l][i_lmia]  = QtyVsZ[(rank*n3) + l][i_lmia] + lmiavsz[l];
      QtyVsZ[(rank*n3) + l][i_cnse]  = consevsz[l];

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::reportEnergyQs ( stack& run ) {

  int calcqvz; run.palette.fetch("calcqvz", &calcqvz);

  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank  );

  if (calcqvz == 1) {
     if (rank == 0) {

       std::cout << "reportEnergyQs: reporting energies..." << std::endl;
       RealArray&  EnergyQs         = run.EnergyQs;

       std::string prefix;    run.palette.fetch(   "prefix",    &prefix   );
       std::string run_label; run.palette.fetch(   "run_label", &run_label);
       std::string res_str;   run.stack_data.fetch("res_str",   &res_str  );

       std::string energy_data_file = prefix + '_' + res_str + ".o" + run_label;
       const char *c_data_file      = energy_data_file.c_str();

       std::ofstream ofs;
       ofs.open( c_data_file, std::ios::out | std::ios::app );

       if (ofs.good()) {

         unsigned esize             = EnergyQs.size();
         for (unsigned k = 0; k < esize; k++) {

           ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << EnergyQs[k] << " ";

         }
         ofs   << std::endl;
         ofs.close();

       }
       else {std::cout << "reportEnergyQs: Warning - could not open file " << energy_data_file << std::endl;}
     }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::reportPowerSpectra ( stack& run ) {

  int calcsvz; run.palette.fetch("calcsvz", &calcsvz);
  if (calcsvz == 1) {

    int isp;physics_data.fetch(    "isp",     &isp);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int         nw; run.palette.fetch("nw",   &nw     );
    int         n3; run.palette.fetch("p3",   &n3     );
    int         np; run.palette.fetch("np",   &np     );
    RealVar  t_cur; physics_data.fetch("t_cur", &t_cur);
    std::string data_dir; run.palette.fetch("data_dir", &data_dir);

    RealArray& z = run.z;

    const char *c_spc_vs_z_out_file;
    std::ofstream ofs;

    RealVar nw_m1 = one / (RealVar (nw));

    std::string spout_pref; run.palette.fetch(    "spout_pref", &spout_pref     );
    std::string res_str;    run.stack_data.fetch( "res_str",    &res_str        );
    std::string run_label;  run.palette.fetch(    "run_label",  &run_label      );
    int srun;               run.palette.fetch(    "srun",       &srun           );

    std::string spc_vs_z_out_file_prefix = "./" + data_dir + "/" + spout_pref + "_" + res_str + ".";
    std::string spc_vs_z_out_file;

    for (unsigned l = 0; l < n3; l++) {

         std::string lyr_str   = static_cast<std::ostringstream*>( &(std::ostringstream() << l+1)  ) -> str();
         int lyr_len           = lyr_str.length();
         if (lyr_len <  3) { for (unsigned i = lyr_len; i < 3;++i){ lyr_str = '0' + lyr_str; }}
         std::string  rank_str = static_cast<std::ostringstream*>( &(std::ostringstream() << rank) ) -> str();
         int rank_len          = rank_str.length();
         if (rank_len < 3) { for (unsigned i = rank_len; i < 3;++i){ rank_str = '0' + rank_str; }}
         std::string srn_str   = static_cast<std::ostringstream*>( &(std::ostringstream() << srun) ) -> str();

         spc_vs_z_out_file     = spc_vs_z_out_file_prefix + rank_str +"_" + lyr_str + ".o" + run_label + srn_str;
         c_spc_vs_z_out_file   = spc_vs_z_out_file.c_str();
         ofs.open( c_spc_vs_z_out_file, std::ios::out );

         if (ofs.good() ) {
           ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << t_cur  << std::endl;
           ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << z[l+1] << std::endl;

           for (unsigned r=0; r < isp+1; ++r) {
             for (unsigned s=0; s < 11; ++s){
               if (s > 0) { SpcVsZ[r][l][s] = nw_m1 * SpcVsZ[r][l][s]; }
                 ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << SpcVsZ[r][l][s] << " ";
             }
             ofs << std::endl;
           }
         }
         ofs.close();
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::reportQtyVsZ ( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  MPI_Status status;

  int   nw; run.palette.fetch("nw",   &nw    );
  int   n3; run.palette.fetch("p3",   &n3    );
  int   np; run.palette.fetch("np",   &np    );
  RealVar  t_cur; physics_data.fetch("t_cur", &t_cur);

  const char *c_qty_vs_z_out_file;
  std::ofstream ofs;

  RealVar nw_m1 = one / (RealVar (nw));

 /*  ~ It should be possible to do the following with an MPI_Gather call, but for the life of me
  *  ~ I can't make it work.
  */

 for (unsigned l = 0; l < n3; l++) {

   if (rank != 0 ) { MPI_Send(&QtyVsZ[rank*n3 + l].front(), 25, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); }
   else { for (unsigned rnk_k = 1; rnk_k < np; rnk_k++) {
            MPI_Recv(&QtyVsZ[rnk_k*n3+l].front(), 25, MPI_DOUBLE, rnk_k, rnk_k, MPI_COMM_WORLD, &status);
          }
        }
 }

  if (rank == 0 ) {

    std::string qout_prefix;  run.palette.fetch(   "qout_pref", &qout_prefix);
    std::string res_str;      run.stack_data.fetch("res_str",   &res_str    );
    std::string run_label;    run.palette.fetch(   "run_label", &run_label  );
    int srun;                 run.palette.fetch(   "srun",      &srun       );
    std::string data_dir; run.palette.fetch("data_dir", &data_dir);
    std::string srn_str = static_cast<std::ostringstream*>( &(std::ostringstream() << srun) ) -> str();

    std::string qty_vs_z_out_file = "./" + data_dir + "/" + qout_prefix + "_" + res_str + ".o" + run_label + srn_str;
    c_qty_vs_z_out_file = qty_vs_z_out_file.c_str();

    std::cout << "reportQtyVsZ: qty_vs_z_out_file = " << qty_vs_z_out_file << std::endl;

    ofs.open( c_qty_vs_z_out_file, std::ios::out );

    if (ofs.good() ) {
      ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << t_cur << std::endl;
      for (unsigned l = 0; l < np*n3; l++) {
        for (unsigned k = 0; k < 26; k++) { 
          if (k > 0) { QtyVsZ[l][k] = nw_m1 * QtyVsZ[l][k]; }
          ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << QtyVsZ[l][k] << " ";
        }
        ofs << std::endl;
      }
    }
    ofs.close();
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalKineticEnergy ( stack& run ) {

  int        rank;  MPI_Comm_rank(MPI_COMM_WORLD,  &rank  );
  
  int        np;    run.palette.fetch(   "np",     &np    );

  int        n1n2c; run.stack_data.fetch("n1n2c",  &n1n2c );
  int        iu2;   run.stack_data.fetch("iu2"   , &iu2   );
  int        n3;    run.stack_data.fetch("n3"   ,  &n3    );
  double     dz;    run.stack_data.fetch("dz"   ,  &dz    );

  double     pe     = zero;
  double     pe_sum = zero;
  int        idx    = 0;
  int        kstart;
  int        kstop;

  RealArray& k2     = run.k2;

  if (rank  == 0) { kstart = 0;     }
  else            { kstart = n1n2c; }

  kstop             = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++){

    if (kdx % n1n2c == 0) { idx = 0; }

      pe            = pe +        k2[idx] * pow(std::abs(P[kdx]), 2);

      if (( rank    == np - 1 ) && ( kdx >= (four * n1n2c) )) {
        pe          = pe + half * k2[idx] * pow(std::abs(P[kdx]), 2);
      }
      if ((rank     == 0      ) && ( kdx <   n1n2c         )) {
        pe          = pe - half * k2[idx] * pow(std::abs(P[kdx]), 2);
      }

      ++idx;
  }

//  pe              = two_thirds * pe * dz;
  pe                = pe * dz;

  int i_red = MPI_Reduce(&pe, &pe_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return pe_sum;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalVorticitySqd( stack& run ) {

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  RealArray& k2 = run.k2;

  int np;    run.palette.fetch(   "np",    &np    );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;   run.stack_data.fetch("iu2"   , &iu2  );
  int n3;    run.stack_data.fetch("n3"   , &n3    );
  double dz; run.stack_data.fetch("dz"   , &dz    );

  double oe     = zero;
  double oe_sum = zero;

  int idx       = 0;

  int kstart;
  int kstop;

  if (rank  == 0) { kstart = 0;     }
  else            { kstart = n1n2c; }

  kstop         = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++){

    if (kdx % n1n2c  == 0) { idx = 0; }

      oe        = oe +        k2[idx] * k2[idx] * pow(std::abs(P[kdx]), 2);

      if (( rank == np - 1 ) && ( kdx >= (four * n1n2c) )) {
        oe      = oe + half * k2[idx] * k2[idx] * pow(std::abs(P[kdx]), 2);
      }
      if ((rank  == 0      ) && ( kdx <   n1n2c)           ) {
        oe      = oe - half * k2[idx] * k2[idx] * pow(std::abs(P[kdx]), 2);
      }

      ++idx;
  }

//  oe = two * two_thirds * oe * dz;  /* ~ NOTE!: When testing against reconnection scenario this ~ */
                                      /*          seems necessary but at moment I don't know why  ~ */

  oe            = two *  oe * dz;

  int i_red     = MPI_Reduce(&oe, &oe_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return oe_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalMagneticEnergy ( stack& run ) {

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  RealArray& k2       = run.k2;

  int np;    run.palette.fetch(   "np",     &np    );
  int n1n2c; run.stack_data.fetch("n1n2c",  &n1n2c );
  int iu2;   run.stack_data.fetch("iu2"   , &iu2   );

  int n3;    run.stack_data.fetch("n3"   ,  &n3    );
  double dz; run.stack_data.fetch("dz"   ,  &dz    );

  double me           = zero;
  double me_sum       = zero;

  int idx             = 0;

  int kstart          = n1n2c; 
  int kstop           = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    if (kdx % n1n2c  == 0) { idx = 0; }

//    me              = me + k2[idx] * pow(std::abs(A[kdx]), 2);
      me              = me + k2[idx] * ( (A[kdx].real()*A[kdx].real()) + (A[kdx].imag()*A[kdx].imag()) );
      ++idx;
  }

  me                  = me * dz;

  int i_red           = MPI_Reduce(&me, &me_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return me_sum;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalCurrentSqd( stack& run ) {

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  RealArray& k2       = run.k2;

  int    np;    run.palette.fetch(   "np",     &np    );
  int    n1n2c; run.stack_data.fetch("n1n2c",  &n1n2c );
  int    iu2;   run.stack_data.fetch("iu2",    &iu2   );
  int    n3;    run.stack_data.fetch("n3",     &n3    );
  double dz;    run.stack_data.fetch("dz",     &dz    );

  double ce           = zero;
  double ce_sum       = zero;

  int idx             = 0;
  int kstart          = n1n2c; 
  int kstop           = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    if (kdx % n1n2c  == 0) { idx = 0; }

//    ce              = ce + k2[idx] * k2[idx] * pow(std::abs(A[kdx]), 2);
      ce              = ce + k2[idx] * k2[idx] * ((A[kdx].real() * A[kdx].real()) + (A[kdx].imag()*A[kdx].imag()));
      ++idx;
  }

  ce                  = two * ce * dz;

  int i_red           = MPI_Reduce(&ce, &ce_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return ce_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalGradCurrentSqd( stack& run ) {

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  RealArray& k2       = run.k2;

  int np;    run.palette.fetch(   "np",    &np    );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;   run.stack_data.fetch("iu2",   &iu2   );
  int n3;    run.stack_data.fetch("n3",    &n3    );
  double dz; run.stack_data.fetch("dz",    &dz    );

  double cee          = zero;
  double cee_sum      = zero;

  cee                 = zero;

  int idx             = 0;
  int kstart          = n1n2c; 
  int kstop           = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    if (kdx % n1n2c  == 0) { idx = 0; }

//    cee             = cee + k2[idx] * k2[idx] * k2[idx] * pow(std::abs(A[kdx]), 2);
      cee             = cee + k2[idx] * k2[idx] * k2[idx] * ( (A[kdx].real()*A[kdx].real()) +(A[kdx].imag()*A[kdx].imag()) );
      ++idx;
  }

  cee                 = two * cee * dz;

  int i_red = MPI_Reduce(&cee, &cee_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return cee_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalFootPointKE( stack& run ) {

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  RealArray& k2       = run.k2;

  int    np; run.palette.fetch(   "np",    &np    );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int   iu2; run.stack_data.fetch("iu2",   &iu2   );
  int    n3; run.stack_data.fetch("n3",    &n3    );
  double dz; run.stack_data.fetch("dz",    &dz    );

  double fp           = zero;
  double fp_sum       = zero;

  assert(dz > 0.0);

  if (rank     == 0) {

    ComplexArray& O   = run.U0;
    int kstart        = 0;
    int kstop         = n1n2c;

    for (unsigned kdx = kstart; kdx < kstop; kdx++) {

//    fp              = fp + k2[kdx] * pow(std::abs(P[kdx]), 2);
      fp              = fp + k2[kdx] * ((P[kdx].real() * P[kdx].real()) + (P[kdx].imag() * P[kdx].imag()));
//    fp              = fp + pow(std::abs(O[kdx]), 2);
    }

// or
//  fp                = three * fp * dz;
    fp                = fp * dz;

//  fp                = two_thirds * fp;  /* ~ NOTE!:  do I need this? ~ */

  }

  return fp;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

RealVar redhallmhd::evalTotalPoyntingFlux ( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  RealArray& k2       = run.k2;

  int np;    run.palette.fetch(   "np",     &np    );
  int n1n2c; run.stack_data.fetch("n1n2c",  &n1n2c );
  int iu2;   run.stack_data.fetch("iu2",    &iu2   );
  int n3;    run.stack_data.fetch("n3",     &n3    );
  RealVar dz; run.stack_data.fetch("dz",    &dz    );
  RealVar n0; run.palette.fetch("n0",       &n0    ); /* ~ density at z_0 (i.e. z = zl / 2)     ~ */

  std::string model; run.palette.fetch("model",      &model);

  RealVar fe           = zero;
  RealVar fe_sum       = zero;

  RealVar Valf;
  RealVar valfmax; run.palette.fetch("valfmax",  &valfmax );
  RealVar Umean0;

  if (      model.compare("rmhd") == 0) { Valf = one;         Umean0 = zero;          }
  else if ( model.compare("inhm") == 0) {
    if (      rank == 0     ) {           Valf = valfmax;     Umean0 = one / nofz[0]; }
    else if ( rank == np -1 ) {           Valf = valfven[n3]; Umean0 = one / nofz[n3];}
  }

  int idx             = 0;

  int kstart;
  int kstop;

  if       (rank  ==  0)      { kstart = 0;                 }
  else if  (rank  ==  np - 1) { kstart = ( iu2 - 2 )*n1n2c; }

  kstop               = kstart + n1n2c;
  idx                 = 0;

  for (unsigned kdx = kstart; kdx < kstop; kdx++){

      if (( rank == np - 1 ) && ( kdx >= ( n3 * n1n2c) )) {

        fe            = fe + Valf   * k2[idx] * ( (P[kdx].real() * A[kdx].real()) + (P[kdx].imag() * A[kdx].imag()) )  \
                           + Umean0 * k2[idx] * ( (P[kdx].real() * P[kdx].real()) + (P[kdx].imag() * P[kdx].imag()) )  \
                           + Umean0 * k2[idx] * ( (A[kdx].real() * A[kdx].real()) + (A[kdx].imag() * A[kdx].imag()) );
      }

      if ((rank  == 0      ) && ( kdx <   n1n2c)           ) {

        fe            = fe - Valf   * k2[idx] * (   (P[kdx].real()        * A[kdx+n1n2c].real())    \
                                                  + (P[kdx].imag()        * A[kdx+n1n2c].imag()) )  \
                           + Umean0 * k2[idx] * (   (P[kdx].real()        * P[kdx].real())          \
                                                  + (P[kdx].imag()        * P[kdx].imag()) )        \
                           + Umean0 * k2[idx] * (   (A[kdx+n1n2c].real()  * A[kdx+n1n2c].real())    \
                                                  + (A[kdx+n1n2c].imag()  * A[kdx+n1n2c].imag()) );
      }
      ++idx;
  }

  fe                  = two * Valf * fe;

  int i_red           = MPI_Reduce(&fe, &fe_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return fe_sum;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalHelicalEnergy ( stack& run ) {

  static const int i_z    =  0;
  static const int i_ae   =  1;
  static const int i_pe   =  2;
  static const int i_ch   =  3;
  static const int i_ep   =  4;
  static const int i_em   =  5;
  static const int i_ce   =  6;
  static const int i_oe   =  7;
  static const int i_zp   =  8;
  static const int i_zm   =  9;
  static const int i_nch  = 10;
  static const int i_lc   = 11;
  static const int i_lo   = 12;
  static const int i_lzp  = 13;
  static const int i_lzm  = 14;
  static const int i_tez  = 15;
  static const int i_rez  = 16;
  static const int i_nez  = 17;
  static const int i_rzp  = 18;
  static const int i_rzm  = 19;
  static const int i_zpm  = 20;
  static const int i_lpl  = 21;
  static const int i_lmi  = 22;
  static const int i_lpla = 23;
  static const int i_lmia = 24;
  static const int i_cnse = 25;

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  RealArray& k2       = run.k2;

  int np;    run.palette.fetch(   "np",    &np    );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;   run.stack_data.fetch("iu2",   &iu2   );
  int n3;    run.stack_data.fetch("n3",    &n3    );
  double dz; run.stack_data.fetch("dz",    &dz    );

  double he           = zero;
  double he_sum       = zero;

  for ( unsigned l = 0; l < n3;++l) { 

    he = he + QtyVsZ[rank*n3+l][i_tez] * ( dudz[l+1] - two*umean[l+1]*Elln[l+1] ) \
            + QtyVsZ[rank*n3+l][i_rez] * two * umean[l+1]*EllA[l+1]               \
            - QtyVsZ[rank*n3+l][i_ch]  * (dvalfdz[l+1] + valfven[l+1]*Elln[l+1]);
  }

  int i_red           = MPI_Reduce(&he, &he_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return he_sum;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyBC( std::string str_step, stack& run ) {

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bdrys; run.palette.fetch("bdrys", &bdrys);
  int np;    run.palette.fetch("np"   , &np   );

  if ( rank == 0 || rank == np - 1) {

    if (bdrys > 0 ) { if (!str_step.compare("finalize") == 0) { applyFootPointDrivingBC( str_step, run ); }}
    else            {                                           applyLineTiedBC(         str_step, run );  }

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyFootPointDrivingBC( std::string str_step, stack& run ) {

  int rank;       MPI_Comm_rank(MPI_COMM_WORLD, &rank      );

  int bdrys;      run.palette.fetch(    "bdrys", &bdrys    );
  int np;         run.palette.fetch(    "np"   , &np       );
  RealVar dtau;   run.palette.fetch(    "dtau" , &dtau     );

  int n1n2;       run.stack_data.fetch( "n1n2",  &n1n2     );
  int n1n2c;      run.stack_data.fetch( "n1n2c", &n1n2c    );
  int n1;         run.stack_data.fetch( "n1"   , &n1       );
  int n2;         run.stack_data.fetch( "n2"   , &n2       );
  int iu3;        run.stack_data.fetch("iu3"   , &iu3      );
  int n_layers;   run.stack_data.fetch( "n3"   , &n_layers );

  RealVar t_cur;  physics_data.fetch(   "t_cur", &t_cur    );

  int i_ua = 0;

  int i_cpy;
  int oldnum, num;

  RealVar    ffp, kc;
  RealVar    dummy;
  RealVar    next_real, next_imag;
  ComplexVar tuple;

  RealVar    lowtau; 
  RealVar    bigtau; 

  RealVar    a, b;                                /* ~ i.e. Gilson's "a" and "b" ~ */

  ComplexArray::size_type nc = n1n2c;
  unsigned strt_idx, stop_idx, idx;


  if ( rank  == 0 || rank == (np - 1) ) {         /* ~ pevol "starts" here        ~ */

    RealArray&    k2         = run.k2;
    ComplexArray& U0         = run.U0;
    ComplexArray& Z          = run.U2; 

    if ( rank  == 0 ) {

      run.palette.fetch("oldnumlb", &oldnum );    /* ~ "numold"                   ~ */
      strt_idx               = 0;
      stop_idx               = strt_idx + n1n2c;

      idx = 0;
      for (unsigned k = strt_idx; k < stop_idx; k++) {

        if (k % n1n2c == 0){ idx = 0; }                    /* ~ reset idx when starting new layer           ~ */

        U0[k]                = czero;
        O[k]                 = czero;
        P[k]                 = czero;

        if (iu3 > 2){ 
        Z[k]                 = czero;
        }
        ++idx;
      }

      num                    = oldnum;
     
      while (  (num * dtau) <= t_cur ) { ++num; }

      num                    = num - 1; 
      lowtau                 = (num    * dtau) - t_cur;
      bigtau                 =((num+1) * dtau) - t_cur;
      a                      = cos((pi*lowtau)/( two * dtau)); /* ~ Gilson's "interp" ~ */
      b                      = cos((pi*bigtau)/( two * dtau));

      if (num == oldnum) {
    
        idx        = 0;

        for (unsigned k = strt_idx; k < stop_idx; k++) {
          
          if (k % n1n2c == 0){ idx  = 0; }                    /* ~ reset idx when starting new layer           ~ */
      
          U0[k]    = (a * roldlb[k]) + (b * rnewlb[k]); ++i_ua;
          O[k]     = U0[k];

          if (k2[idx] != 0) { P[k]  = O[k] / k2[idx]; }
          if (iu3 > 2)      { Z[k]  = czero;          }

          ++idx;
        }
      }   // num is oldnum
      else {

        int brcount; physics_data.fetch("brcount", &brcount);

        run.palette.fetch("ffp", &ffp);
        run.palette.fetch( "kc", &kc );

        for (unsigned k = 0; k < num - oldnum; k++) {

          for (unsigned l = 0; l < nc; l++) { roldlb[l] = rnewlb[l]; }

          i_ua = 0;
          for (unsigned l = 0; l < nc; l++ ) {

/* ~ fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  */

            if (l  < (n1n2c - (n2/2)-1)) {
              if ((l == 0) || (l %((n2/2)+1) !=0) ) {

                if ( sqrt(k2[l]) <= kc ) {

                    dummy        = ((double) rand() / RAND_MAX ); ++brcount; // why not just here?
                    dummy        = ((double) rand() / RAND_MAX ); ++brcount;

                    next_real    = ffp * (((double) rand() / RAND_MAX) * two - one);
                    ++brcount;
                    next_imag    = ffp * (((double) rand() / RAND_MAX) * two - one);
                    ++brcount;
                    tuple        = ComplexVar(next_real, next_imag);
                    rnewlb[l]    = tuple;

                }
                if ( l == 0) { rnewlb[l] = czero; }
                ++i_ua;
              } // not copying
              else {
                i_cpy    = (l/((n2/2)+1));
//              if (i_cpy < 0 || i_cpy > n2/2-1) {
//                std::cout << "applyFoot: WARNING - i_cpy = " << i_cpy << " for l = " << l << std::endl;
//              } // i_cpy out of range
//              else {
                  if ( i_cpy < n2/2 ) {
                    rnewlb[l] = rnewlb[i_cpy];
                  } // really copying
                  else {
                    if ( sqrt(k2[l]) <= kc ) {

                        dummy        = ((double) rand() / RAND_MAX ); ++brcount; // why not just here?
                        dummy        = ((double) rand() / RAND_MAX ); ++brcount;

                        next_real    = ffp * (((double) rand() / RAND_MAX) * two - one);
                        ++brcount;
                        next_imag    = ffp * (((double) rand() / RAND_MAX) * two - one);
                        ++brcount;
                        tuple        = ComplexVar(next_real, next_imag);
                        rnewlb[l]    = tuple;

                    }
                    ++i_ua;
                  } // not really copying
//              }   // i_cpy in range
              }     //copying
            }       // l is <  (n1n2c - (n2/2) - 1
            else {
              i_cpy = (l - n1n2c + n2/2+2);
              if (i_cpy < 0 || i_cpy > n2/2+1) {
                std::cout << "applyFoot: WARNING - i_cpy = "             << i_cpy << " for l = " << l <<  std::endl;
              }
              else {
//              std::cout << "applyFoot: l >= n1n2c - n2/2 <-> i_cpy = " << i_cpy << " for l = " << l << std::endl;
                if ( l != (n1n2c - 1)) {
                  rnewlb[l] = std::conj(rnewlb[i_cpy]);
                }
                else { rnewlb[l] = czero;}
              }
            }       // l is >= (n1n2c - (n2/2) - 1

/* ~ fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  fill rnewlb ~  */

          } // end of for loop in l

          if ( k == 0 ) { std::cout << "applyFoot: i_ua = " << i_ua << std::endl;}
        } // loop in k

        idx = 0;
        for (unsigned k = strt_idx; k < stop_idx; k++) {

          if (k % n1n2c == 0){ idx = 0; }                    /* ~ reset idx when starting new layer           ~ */

          U0[k]              = (a * roldlb[k]) + (b * rnewlb[k]);
          O[k]               = U0[k];

          if (k2[idx] != 0) { P[k]  = O[k] / k2[idx]; }
          if (iu3 > 2){ 
            Z[k]             = czero;
          }
          ++idx;
        } // end second loop in k

        physics_data.reset("brcount", brcount);
      } // num is not oldnum
      run.palette.reset("oldnumlb", num);
    }  //rank is zero
    else if ( rank == (np - 1) && (bdrys == 2)) {

      run.palette.fetch("oldnumub", &oldnum ); /* ~ "oldnum" ~ */
      strt_idx               = n_layers * nc;
      stop_idx               = strt_idx + n1n2c;

      idx = 0;
      for (unsigned k = strt_idx; k < stop_idx; k++) { 
        if (k % n1n2c == 0){ idx = 0; }                    /* ~ reset idx when starting new layer           ~ */
        U0[k]                = czero;
        O[k]                 = U0[k];
        P[k]                 = czero;

        if (iu3 > 2) {
          Z[k]               = czero;
        }
        ++idx;
      }

      num                    = oldnum;

      while (  (num * dtau) <= t_cur ) { ++num; }                 /* ~~~~~~~~~~~~~~~~~~~~~ */
                                                                  /*                       */
      num                    = num - 1;                           /* ~ Gilson's "getpas" ~ */
      lowtau                 =  (num    * dtau) - t_cur;          /*                       */
      bigtau                 = ((num+1) * dtau) - t_cur;          /* ~~~~~~~~~~~~~~~~~~~~~ */

      a                      = cos((pi*lowtau)/( two * dtau));    /* ~ Gilson's "interp" ~ */
      b                      = cos((pi*bigtau)/( two * dtau));

      if (num == oldnum) {
    
        int idx              = 0;
        for (unsigned k = strt_idx; k < stop_idx; k++) {
        if (k % n1n2c == 0){ idx = 0; }                    /* ~ reset idx when starting new layer           ~ */
           U0[k]             = (a * roldub[idx]) + (b * rnewub[idx]);
            O[k]             = U0[k];
            if (k2[idx] != 0) { P[k]  = O[k] / k2[idx]; }

          if (iu3 > 2) {
            Z[k]             = czero;
          }
          ++idx;
        }
      }

      else {

        int trcount; physics_data.fetch("trcount", &trcount);

        run.palette.fetch("ffp",      &ffp);
        run.palette.fetch("kc",       &kc);

        for (unsigned k = 0; k < num - oldnum; k++) {

          for (unsigned l = 0; l < nc; l++) { roldub[l] = rnewub[l]; }

          i_ua = 0;
          for (unsigned l = 0; l < nc; l++ ) {

// ~ fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  */

            if (l  < (n1n2c - (n2/2)-1)) {
              if ((l == 0) || (l %((n2/2)+1) !=0) ) {

                if ( sqrt(k2[l]) <= kc ) {

                    next_real    = ffp * (((double) rand() / RAND_MAX ) * two - one);
                    ++trcount;
                    next_imag    = ffp * (((double) rand() / RAND_MAX ) * two - one);
                    ++trcount;
                    tuple        = std::complex<double>(next_real, next_imag);
                    rnewub[l]    = tuple;

                }
                if ( l == 0) { rnewub[l] = czero; }
                ++i_ua;
              } // not copying
              else {

                i_cpy               = (l/((n2/2)+1));

                if ( i_cpy < n2/2 ) {
                  rnewub[l]         = rnewub[i_cpy];
                } // really copying
                else {
                  if ( sqrt(k2[l]) <= kc ) {

                      dummy         = ((double) rand() / RAND_MAX ); ++trcount; // why not just here?
                      dummy         = ((double) rand() / RAND_MAX ); ++trcount;

                      next_real     = ffp * (((double) rand() / RAND_MAX) * two - one);
                      ++trcount;
                      next_imag     = ffp * (((double) rand() / RAND_MAX) * two - one);
                      ++trcount;
                      tuple         = ComplexVar(next_real, next_imag);
                      rnewub[l]     = tuple;

                  }
                  ++i_ua;
                }   // not really copying
              }     //copying
            }       // l is <  (n1n2c - (n2/2) - 1
            else {
              i_cpy = (l - n1n2c + n2/2+2);
              if (i_cpy < 0 || i_cpy > n2/2+1) {
                std::cout << "applyFoot: WARNING - i_cpy = "             << i_cpy << " for l = " << l <<  std::endl;
              }
              else {
////            std::cout << "applyFoot: l >= n1n2c - n2/2 <-> i_cpy = " << i_cpy << " for l = " << l << std::endl;
                if (l != (n1n2c - 1)) {
                rnewub[l] = std::conj(rnewub[i_cpy]);
                }
                else { rnewub[l] = czero; }
              }
            }       // l is >= (n1n2c - (n2/2) - 1

// ~ fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  fill rnewub ~  */

          } // end of for loop in l

          if ( k == 0 ) { std::cout << "applyFoot: i_ua = " << i_ua << std::endl;}
        }   // loop in k

        int idx              = 0;
        for (unsigned k = strt_idx; k < stop_idx; k++) {

          if (k % n1n2c == 0){ idx = 0; }                    /* ~ reset idx when starting new layer           ~ */
          U0[k]               = (a * roldub[idx]) + (b * rnewub[idx]);
           O[k]               = U0[k];
          if (k2[idx] != 0) { P[k]  = O[k] / k2[idx]; }
          if (iu3 > 2) { Z[k] = czero; }
          ++idx;
        }
        physics_data.reset("trcount", trcount);
      }
        run.palette.reset("oldnumub", num);
    } // rank is np - 1 and bdrys is 2
  }   // rank is one of zero or  np - 1
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyLineTiedBC( std::string str_step, stack& run ) {

  int rank;          MPI_Comm_rank(MPI_COMM_WORLD, &rank       );

  int np;            run.palette.fetch(    "np",     &np       );
  int n1n2c;         run.stack_data.fetch( "n1n2c",  &n1n2c    );
  int n3;            run.palette.fetch(    "p3",     &n3       );
  std::string model; run.palette.fetch(    "model",  &model    );
  double tstart;     run.palette.fetch(    "tstart", &tstart   );

  unsigned strt_idx, stop_idx;

  if ((rank == 0) || rank == (np - 1)) {

    if (     rank    == 0       ) { strt_idx = 0;             }
    else if( rank    == (np - 1)) { strt_idx = ( n3 * n1n2c); }

    stop_idx          = strt_idx + n1n2c;

    ComplexArray& U0  = run.U0;
    ComplexArray& Z   = run.U2;

    ComplexArray& tU0 = run.tU0;
    ComplexArray& tZ  = run.tU2;

    if (str_step.compare("predict") == 0 || str_step.compare("finalize") == 0) {

      for (unsigned k = strt_idx; k < stop_idx; k++) {

        U0[k]         = czero;
        O[k]          = czero;
        P[k]          = czero;

        if (model.compare("hall") == 0 ) { 
          Z[k]        = czero; 
        }
      }
    }
    else if (str_step.compare("correct") == 0 ) {

      for (unsigned k = strt_idx; k < stop_idx; k++) {

        tU0[k]        = czero;
        O[k]          = czero;
        P[k]          = czero;

        if (model.compare("hall") == 0 ) {
          tZ[k]       = czero;
        }
      }
    }
    else { std::cout    << "applyLineTiedBC: WARNING - unknown str_step value " << std::endl; }
  } // rank has exterior boundary
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::updatePAOJ( std::string str_step, stack& run ) {

  int rank;     MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  int  np;      MPI_Comm_size(MPI_COMM_WORLD, &np);

  std::string model; run.palette.fetch("model", &model);
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c);    /* ~ number complex elements in layer            ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers); /* ~ number of layers in stack                   ~ */
  int n3;       run.palette.fetch(   "p3",    &n3);       /* ~ number of layers in stack                   ~ */

  ComplexArray&   U0 = run.U0;                            /* ~ for predictor case                          ~ */
  ComplexArray&    H = run.U1;

  ComplexArray&  tU0 = run.tU0;                           /* ~ for corrector case                          ~ */
  ComplexArray&   tH = run.tU1;

  RealArray&      k2 = run.k2;                            /* ~ square-magnitude of k-space vectors         ~ */
  RealArray&  inv_k2 = run.inv_k2;                        /* ~ inverse square magnitude of k-space vectors ~ */

  RealVar ssqd; run.palette.fetch(  "ssqd",  &ssqd);      /* ~ parameter sigma^2 relating A to H           ~ */

  unsigned kstart;                                        /* ~ lower and upper loop limits on k            ~ */
  unsigned kstop;                                         /* ~ note: pbot and atop boundary layers incl.'d ~ */

  if ( rank != 0 && rank != (np - 1)) {

    kstart    = 0;                                 /* ~ lower and upper loop limits on k            ~ */
    kstop     = n_layers * n1n2c;                  /* ~ note: pbot and atop boundary layers incl.'d ~ */

  }
  else {
    if (rank == 0 ) {

      kstart  = n1n2c;                             /* ~ lower and upper loop limits on k            ~ */
      kstop   = n_layers * n1n2c;                  /* ~ note: pbot and atop boundary layers incl.'d ~ */

    }
    else if (rank == (np - 1)) {

      kstart  = 0;                                 /* ~ lower and upper loop limits on k            ~ */
      kstop   = (n_layers - 1 ) * n1n2c;           /* ~ note: pbot and atop boundary layers incl.'d ~ */

    }
  }

  unsigned idx       = 0;                                 /* ~ index for k2 and inv_k2                     ~ */

  if (    str_step.compare("predict") == 0 ) {            /* ~ PRIOR to predictor step                     ~ */
    for (unsigned k = kstart; k < kstop; k++) { 

       if (k % n1n2c == 0){ idx = 0; }                    /* ~ reset idx when starting new layer           ~ */

         if (rank != (np - 1)) {

         if (k2[idx] != 0) { O[k]  = U0[k];               /* ~ last layer of highest process is an         ~ */
                             P[k]  = O[k] / k2[idx];      /* ~ boundary for P and O, but not A             ~ */
                           }                              /* ~ so we must be careful not to clobber A      ~ */
                                                          /* ~ while preserving boundary conditions        ~ */
         }
         else if ( rank == (np - 1) && k < (n3*n1n2c)){

         if (k2[idx] != 0) { O[k]  = U0[k];               /* ~ last layer of highest process is an         ~ */
                             P[k]  = O[k] / k2[idx];      /* ~ boundary for P and O, but not A             ~ */
                           }                              /* ~ so we must be careful not to clobber A      ~ */

         }

       if ( model.compare("hall") == 0 )  {           /* ~ initialize only if needed  ~ */
         A[k] = H[k] / ( one + ssqd*k2[idx] );
       }
       else { A[k] = H[k]; }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
//     J[k] = std::sqrt(two_thirds)*k2[idx] * A[k];               /* ~ J = -delperp^2 A                            ~ */
       J[k] = k2[idx] * A[k];                                     /* ~ J = -delperp^2 A                            ~ */
/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

       ++idx;

    }
  } // predict
  else if(str_step.compare("correct") == 0 ) {            /* ~ do as above using predictor results         ~ */
    for (unsigned k = kstart; k < kstop; k++) {

      if (k % n1n2c == 0){ idx = 0;}

         if (rank != (np - 1) ) {

         if (k2[idx] != 0) { O[k]  = tU0[k];
                             P[k]  = O[k] / k2[idx]; 
                           }
         }
         else if ( rank == (np - 1) && k < (n3*n1n2c)) {

         if (k2[idx] != 0) { O[k]  = tU0[k];
                             P[k]  = O[k] / k2[idx]; 
                           }

         }

      if ( model.compare("hall") == 0 )  {           /* ~ initialize only if needed  ~ */
      A[k] = tH[k] / ( one + ssqd*k2[idx] );
      }
      else { A[k] = tH[k]; }
/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
      J[k] = k2[idx] * A[k];
//    J[k] = std::sqrt(two_thirds)*k2[idx] * A[k];
/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

      ++idx;
    }
  } // correct
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::updateTimeInc( stack& run ) {

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank );

    RealArray glbMaxU;
    glbMaxU.assign(maxU.size(),zero);

    int n1n2;   run.stack_data.fetch("n1n2", &n1n2);
    RealVar dt; run.palette.fetch(   "dt",   &dt  );
    RealVar dz; run.stack_data.fetch("dz",   &dz  );
 
    RealVar q1; run.palette.fetch(   "q1",   &q1  );
    RealVar q2; run.palette.fetch(   "q2",   &q2  );
    RealVar qp; run.palette.fetch(   "qp",   &qp  );

    RealVar dtvb;
    RealVar dtr;

    int i_red   = MPI_Reduce(&maxU[0], &glbMaxU[0], maxU.size(), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
    int i_brd   = MPI_Bcast(&glbMaxU[0], maxU.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (unsigned i_f = 0; i_f < maxU.size(); i_f++) { maxU[i_f] = glbMaxU[i_f]; }

    dtvb        = zero;

    for (unsigned i_f = 0; i_f < maxU.size(); i_f++) { dtvb = dtvb + sqrt(maxU[i_f]); }

    dtvb        = dtvb * sqrt( (RealVar) (n1n2) );
    dtvb        = one / dtvb;
    dtr         = dtvb / dt;

    physics_data.reset("dtvb", dtvb);

    if (( dt < (q1 * dtvb) ) && ( dt <  (half * qp * dz) ) ) {

      dt        = two * dt;
      run.palette.reset("dt", dt);

    }
    else if ( dt > (q2 * dtvb) ) {

      dt        = half * dt;
      run.palette.reset("dt", dt);

    }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::checkState( int pair, stack &run, std::string roc) {

  int          rank  ; MPI_Comm_rank(MPI_COMM_WORLD ,   &rank  );
  int          srun  ; run.palette.fetch(    "srun" ,   &srun  );
  int          n1n2  ; run.stack_data.fetch( "n1n2" ,   &n1n2  );
  int          n1n2c ; run.stack_data.fetch( "n1n2c",   &n1n2c );
  int          iu2   ; run.stack_data.fetch( "iu2"  ,   &iu2   );

  unsigned     usize = iu2*n1n2c;
  ComplexArray Q0(usize,czero);
  ComplexArray Q1(usize,czero);

  RealArray RealQ0(iu2*n1n2,zero);
  RealArray RealQ1(iu2*n1n2,zero);

  RealArray& k2 = run.k2;

  switch(pair) {
    case(0) :
      for (unsigned k = 0; k < usize; ++k) { Q0[k] = P[k];      Q1[k] = A[k];      }
      break;
    case(1) :
      for (unsigned k = 0; k < usize; ++k) { Q0[k] = O[k];      Q1[k] = A[k];      }
      break;
    case(2) :
      for (unsigned k = 0; k < usize; ++k) { Q0[k] = run.U0[k]; Q1[k] = run.U1[k]; }
      break;
    case(3) :
      for (unsigned k = 0; k < usize; ++k) { Q0[k] = O[k];      Q1[k] = J[k];      }
      break;
    case(4) :
      for (unsigned k = 0; k < usize; ++k) { Q0[k] = A[k];      Q1[k] = J[k];      }
      break;
    case(5) :
      for (unsigned k = 0; k < usize; ++k) { Q0[k] = P[k];      Q1[k] = O[k];      }
      break;
  }

  std::string str_srun = static_cast<std::ostringstream*>( &(std::ostringstream() << srun)  ) -> str();
  std::string str_rank = static_cast<std::ostringstream*>( &(std::ostringstream() << rank)  ) -> str();
  std::string state_pa = "state_pa_rnk-" + str_rank + "_srun-" + str_srun;
  const char *c_state_pa;
  RealVar next;
  c_state_pa           = state_pa.c_str();
  std::ofstream ofs;

  if (rank == 0 ) { std::cout << "checkState: smallish = " << smallish << std::endl;}

  ofs.open( c_state_pa, std::ios::out );

  if (roc.compare("c") == 0) {
    int idx = 0;
    for (unsigned k = 0; k < iu2*n1n2c;++k) {

      if (k % n1n2c == 0){idx = 0;}

      if (std::abs(Q0[k].real()) > smallish ){
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << Q0[k].real() << " ";
      }
      else {
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << zero         << " ";
      }
      if (std::abs(Q0[k].imag()) > smallish ){
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << Q0[k].imag() << " ";
      }
      else {
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << zero         << " ";
      }
      if (std::abs(Q1[k].real()) > smallish ){
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << Q1[k].real() << " ";
      }
      else {
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << zero         << " ";
      }
      if (std::abs(Q1[k].imag()) > smallish ){
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << Q1[k].imag() << " ";
      }
      else {
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << zero         << " ";
      }
      ofs   << std::setw(26) << std::right << std::setprecision(18) << std::scientific <<  k2[idx] << std::endl;

      ++idx;
    }
  }
  else if (roc.compare("r") == 0) {

    fftw.fftwReverseRaw(run, Q0, RealQ0); 
    fftw.fftwReverseRaw(run, Q1, RealQ1); 

    int idx = 0;
    for (unsigned k = 0; k < iu2*n1n2;++k) {

      if (k % n1n2 == 0){idx = 0;}
      if (std::abs(RealQ0[k]) > smallish ){
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << RealQ0[k]    << " ";
      }
      else {
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << zero         << " ";
      }
      if (std::abs(RealQ1[k]) > smallish ){
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << RealQ1[k]    << " ";
      }
      else {
        ofs << std::setw(26) << std::right << std::setprecision(18) << std::scientific << zero         << " ";
      }
      ofs   << std::endl;

      ++idx;
    }
  }

  ofs.close();
}

#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::physicsFinalize ( stack& run) {

  finalizeBoundaries( run );

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~ Destructor ~ */

redhallmhd:: ~redhallmhd() {

//  fftw.fftwFinalize();

}
