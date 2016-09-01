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

  init_physics_data(     run   );   /* ~ physics - specific parameters               ~ */
  initU(                 run   );   /* ~ initialization of layers 1 - n3 of U        ~ */
                                    /* ~ boundary value application.                 ~ */

  initialize(            run   );   /* ~ not a good name, I'll probably revise this  ~ */
                                    /* ~ it might be to let initialize do everything ~ */
                                    /* ~ here or, alternatively to do away with it   ~ */

  initBoundaries(        run   );   /* ~ initialization of quantities needed for     ~ */

  evalValf(              run   );

  int srun; run.palette.fetch( "srun", &srun );

  if (srun == 1) {

    fftw.fftwReverseAll( run, J);
    run.writeUData  (          );   /* ~ initial conditions report                   ~ */

  }
  else { /* ~ NOTE! AUX should be read from data file here!                            ~ */ }

  initIRMHD(             run   );   /* ~ checks for inhomegeneous conditions         ~ */

#endif

}

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

    if (srun == 1) { run.zeroU();      }
    else           { readUData( run ); }

  }

  int n3;      run.palette.fetch("p3", &n3);
  int np;      run.palette.fetch("np", &np);

  int calcqvz; run.palette.fetch("calcqvz", &calcqvz);
  int calcsvz; run.palette.fetch("calcsvz", &calcsvz);

  if (calcqvz == 1) {

    RealArray Qlayer; Qlayer.assign(15,zero); QtyVsZ.assign(np*n3, Qlayer);

  }

  if (calcsvz == 1) {

    int n1;      run.stack_data.fetch("n1", &n1);
    int n2;      run.stack_data.fetch("n2", &n2);

    isp    = n2/2;
    dk     = (pi * n1) / isp;
    dk_m1  = one / dk;

    kb     = two * pi;
    kf     = kb + (isp * dk);

    ikb    = 1 + (int) (kb * dk_m1);
    ikf    = 1 + (int) (kf * dk_m1);

    nk  = ikf - ikb + 1;

//  RealArray     SpcLayer; SpcLayer.assign(n3,zero);
//  Real2DArray   SpcSet;   SpcSet.assign(7,SpcLayer); 
//  SpcVsZ.assign(isp,SpcSet);

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

  for (int i_f = 0; i_f < n_flds; ++i_f) {
    for (int i_l = 1; i_l < n_lyrs + 1; ++i_l) {
    
    for (int k = 0; k< n1n2; ++k) { U[k][i_l][i_f] = U[k][n3][i_f]; }

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
       real_part      =    -0.1L;
       imag_part      =     0.0L;
       tuple          = ComplexVar(real_part, imag_part);
       Cin[idx]       = tuple;

       idx            =  (n1-1) * (n2/2 + 1) + 1;
       real_part      =     0.1L;
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

  for (  int i_f = 0; i_f < n_flds; ++i_f) {

    for (int i_l = 1; i_l < n_lyrs ; ++i_l) {
    
    for (int k   = 0; k< n1n2; ++k) { U[k][i_l][i_f] = U[k][n3][i_f]; }

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::readUData( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int srun; run.palette.fetch("srun", &srun);

  std::string data_file                  = run.getLastDataFilename(srun-1);
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

      U[to_row_maj_idx][slab_index][0]   = next_p;
      U[to_row_maj_idx][slab_index][1]   = next_a;

      if(iu3 <= 2) {
        ifs >> next_o;
        ifs >> next_j;

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

    if (bdrys > 0) { finalizeFootPointDriving( run ); }
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

  MPI_Status * status   = 0;

  int tag_old           = 0;
  int tag_new           = 1;

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

  int srun;   run.palette.fetch(    "srun",   &srun  );
  int n1n2c;  run.stack_data.fetch( "n1n2c",  &n1n2c );
  RealVar kc; run.palette.fetch(    "kc",     &kc    );
  int np;     run.palette.fetch(    "np",     &np    );
  int bdrys;  run.palette.fetch( "bdrys",    &bdrys  );

  RealArray& k2         = run.k2;

  RealVar    next_real; 
  RealVar    next_imag;
  ComplexVar tuple;

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

        if ( sqrt(k2[l]) < kc ) {

            next_real   = ffp * ( ((double) rand() / RAND_MAX ) * two - one); ++brcount;
            dummy       =         ((double) rand() / RAND_MAX ) * two - one;  ++brcount;
            next_imag   = ffp * ( ((double) rand() / RAND_MAX ) * two - one); ++brcount;
            dummy       =         ((double) rand() / RAND_MAX ) * two - one;  ++brcount;

            tuple       = ComplexVar(next_real, next_imag);
            rnewlb[l]   = tuple;

        }
    }
  
    l_reset_success     = physics_data.reset("brcount", brcount);

    }
    else {
      
      std::string prefix;    run.palette.fetch(    "prefix",     &prefix    );
      std::string run_label; run.palette.fetch(    "run_label",  &run_label );
      std::string res_str;   run.stack_data.fetch( "res_str",    &res_str   );

      std::string srn_str              = static_cast<std::ostringstream*>( &(std::ostringstream() << ( srun - 1 ) ) ) -> str();
      std::string boundary_data_file   = prefix + '_' + res_str + "r" + srn_str;
      const char *c_boundary_data_file = boundary_data_file.c_str();

      std::ifstream ifs;
      ifs.open( c_boundary_data_file );

      ComplexVar cignore;
      if ( ifs.good() ) {
        int fld_count                  = 0;
        int line_count                 = 0;
        while ( !ifs.eof() ) {                             /* ~ order: roldlb rnewlb roldub rnewub pbot ~ */

          ifs >> next_real; ifs >> next_imag; tuple = ComplexVar(next_real, next_imag); roldlb[line_count] = tuple; ++fld_count;
          ifs >> next_real; ifs >> next_imag; tuple = ComplexVar(next_real, next_imag); rnewlb[line_count] = tuple; ++fld_count;
          ifs >> next_real; ifs >> next_imag; tuple = ComplexVar(next_real, next_imag); roldub[line_count] = tuple; ++fld_count;
          ifs >> next_real; ifs >> next_imag; tuple = ComplexVar(next_real, next_imag); rnewub[line_count] = tuple; ++fld_count;
          ifs >> next_real; ifs >> next_imag; tuple = ComplexVar(next_real, next_imag); cignore            = tuple; ++fld_count;

          if (fld_count == 5) { ++line_count; fld_count = 0;}

        }
        ifs.close();
        assert( line_count == n1n2c + 1);

        MPI_Send( &roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np - 1, tag_old, MPI_COMM_WORLD );
        MPI_Send( &rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np - 1, tag_new, MPI_COMM_WORLD );

      }
      else { std::cout << "initFootPointDriving: WARNING - could not open file " << boundary_data_file << std::endl; }

    }
  }

  if ( rank == np - 1 ) {
    if (bdrys  == 2) {

      int trcount; physics_data.fetch("trcount", &trcount);

      for (int l = 0; l < n1n2c; ++l) {

        if ( sqrt(k2[l]) < kc ) {

          if (srun == 1) {

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
          else {

            dummy       =        ((double) rand() / RAND_MAX );              ++trcount;
            dummy       =        ((double) rand() / RAND_MAX );              ++trcount;

          }
        }
      }

      l_reset_success   = physics_data.reset("trcount", trcount);
    }
    if ( srun > 1) {

       MPI_Recv(&roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_old, MPI_COMM_WORLD, status);
       MPI_Recv(&rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_new, MPI_COMM_WORLD, status);

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::finalizeFootPointDriving( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status * status = 0;

  int  rcount;  run.palette.fetch(  "rcount",  &rcount  );
  int brcount;  physics_data.fetch( "brcount", &brcount );
  int trcount;  physics_data.fetch( "trcount", &trcount );

  rcount = rcount + ( (brcount) > (trcount) ? (brcount) : (trcount) );

  run.palette.reset("rcount", rcount);

  int np;     run.palette.fetch(    "np",      &np      );
  int n1n2c;  run.stack_data.fetch( "n1n2c",   &n1n2c   );
  int srun;   run.palette.fetch(    "srun",    &srun    );

  int tag_old         = 0;
  int tag_new         = 1;

  if (rank            == 0)      {
     
    MPI_Recv(&roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np -1, tag_old, MPI_COMM_WORLD, status);
    MPI_Recv(&rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, np -1, tag_new, MPI_COMM_WORLD, status);

    std::string prefix;    run.palette.fetch(    "prefix",     &prefix    );
    std::string run_label; run.palette.fetch(    "run_label",  &run_label );
    std::string res_str;   run.stack_data.fetch( "res_str",    &res_str   );

    std::string srn_str              = static_cast<std::ostringstream*>( &(std::ostringstream() << ( srun ) ) ) -> str();
    std::string boundary_data_file   = prefix + '_' + res_str + "r" + srn_str;
    const char *c_boundary_data_file = boundary_data_file.c_str();

    std::ofstream ofs;
    ofs.open( c_boundary_data_file, std::ios::out | std::ios::trunc );

    if ( ofs.good() ) {
      for (unsigned k = 0; k < n1n2c; k++) {               /* ~ order: roldlb rnewlb roldub rnewub pbot ~ */

        
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << roldlb[k].real() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << roldlb[k].imag() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << rnewlb[k].real() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << rnewlb[k].imag() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << roldub[k].real() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << roldub[k].imag() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << rnewub[k].real() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << rnewub[k].imag() << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << P[k].real()      << " ";
        ofs << std::setw(30) << std::right << std::setprecision(16) << std::scientific << P[k].imag()      << " ";
        ofs << std::endl;

      }
      ofs.close();
    }
    else { std::cout << "finalizeFootPointDriving: WARNING - could not open file " << boundary_data_file << std::endl; }
  }
  else if (rank       == np - 1) {

    assert (roldub.size() == n1n2c);
    assert (rnewub.size() == n1n2c);

    MPI_Send( &roldub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_old, MPI_COMM_WORLD );
    MPI_Send( &rnewub.front(), n1n2c, MPI::DOUBLE_COMPLEX, 0, tag_new, MPI_COMM_WORLD );

  }
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

void redhallmhd::initIRMHD (stack& run ) {

  evalValf(  run );
  evalUmean( run );
  evalElls(  run );

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initialize (stack& run ) {

  std::string model;
  run.palette.fetch("model", &model);

#ifdef HAVE_CUDA_H
 
#else

  initTimeInc(         run );
  fftw.fftwForwardAll( run );

  OfromP(              run );
  HfromA(              run );

#endif

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::OfromP( stack& run )  {

    int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c    ); /* ~ number of complex elements per layer       ~ */
    int n_layers; run.stack_data.fetch("iu2",   &n_layers ); /* ~ number of layers in stack                  ~ */
  
    RealArray&    k2 = run.k2;                               /* ~ square magnitudes of k-space vectors       ~ */
    ComplexArray& U0 = run.U0;                               /* ~ holds phi (i.e. P ) at this point          ~ */ 
  
    ComplexArray::size_type usize;
    usize            = U0.capacity();                        /* ~ current capacity of U0 - should be known   ~ */
  
    assert(usize     == (n1n2c * n_layers));                 /* ~ test usize                                 ~ */
  
    ComplexArray O;                                          /* ~ temporary storage for vorticity            ~ */
    O.assign(usize,czero);
  
    P.assign(usize,czero);                                   /* ~ member P will be needed throughout run     ~ */
  
   for (unsigned k = 0; k <usize; k++) {P[k] = U0[k];}       /* ~ preserve stream funtion in P               ~ */
   
    unsigned  idx    = 0;                                    /* ~ index for k2                               ~ */
    for (unsigned k  = 0; k < usize; k++) {
  
      if (k % n1n2c  == 0 ) { idx = 0; }                     /* ~ reset idx when starting new layer          ~ */
      O[k] = k2[idx] * P[k];                                 /* ~ Omega = - delperp^2 P                      ~ */
  
      ++idx;
  
    }
  
   for (unsigned k = 0; k <usize; k++) {U0[k] = O[k];}       /* ~ U0 now holds Fourier transform of Vorticity ~ */
                                                             /* ~ and P is initialized                        ~ */
    O.resize(0);                                             /* ~ dispense with temporary storage             ~ */
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::HfromA( stack& run )  {

    int rank;     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   );  /* ~ number of complex elements per layer       ~ */
    int n_layers; run.stack_data.fetch("iu2",   &n_layers);  /* ~ number of layers in stack                  ~ */
  
    RealArray& k2    = run.k2;                               /* ~ square magnitude of k-space vectors        ~ */
    ComplexArray& U1 = run.U1;                               /* ~ holds A (flux function) at this point      ~ */
  
    ComplexArray::size_type usize;
    usize            = U1.capacity();                        /* ~ current capacity of U1 - should be known   ~ */
  
    assert(usize     == (n1n2c * n_layers));                 /* ~ test usize                                 ~ */
  
    ComplexArray H;                                          /* ~ temporary storage for H-function           ~ */
    H.reserve(usize);
  
    A.reserve(usize);                                        /* ~ members A and J (current density) needed   ~ */
    J.reserve(usize);                                        /* ~ throughout run                             ~ */
  
    RealVar ssqd;      run.palette.fetch(  "ssqd",  &ssqd ); /* ~ parameter sigma^2 needed for H             ~ */
    std::string model; run.palette.fetch(  "model", &model); 
  
    for (unsigned k = 0; k < usize; k++) {A[k] = U1[k];}     /* ~ preserve flux function in A                ~ */

    unsigned idx     = 0;                                    /* ~ index for k2                               ~ */
    for (unsigned k  = 0; k < usize; k++) {
  
      if (k % n1n2c  == 0 ) { idx = 0; }                     /* ~ reset idx when starting new layer          ~ */
  
      J[k]           = k2[idx] * A[k];                       /* ~ J = -delperp A                             ~ */
      if (model.compare("hall") == 0 ) {
        H[k]         = A[k] + ssqd * J[k];                   /* ~ H = A + sigma^2 J                          ~ */
      }
      else if (model.compare("rmhd") == 0 || model.compare("inhm") == 0) {
        H[k]         = A[k];                                 /* ~ H = A for reduced-MHD                      ~ */
      }
      ++idx;
  
    }
  
    for (unsigned k = 0; k < usize; k++) {U1[k] = H[k];}     /* ~ U1 now holds Fourier transform of H        ~ */
                                                             /* ~ and A and J are both initialized           ~ */
    H.resize(0);                                             /* ~ dispense with temporary storage            ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::PfromO( stack& run )  {
 
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c      );    /* ~ number of complex elements per layer        ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers);    /* ~ number of layers in stack                   ~ */

  RealArray&    inv_k2 = run.inv_k2;                         /* ~ inverse-square magnitude of k-space vectors ~ */
  ComplexArray& U0     = run.U0;                             /* ~ holds Omega (vorticity) at this point       ~ */

  ComplexArray::size_type usize;
  usize                = U0.capacity();                      /* ~ current capacity of U0 - should be known    ~ */

  assert(usize         == (n1n2c * n_layers));               /* ~ test usize                                  ~ */

  ComplexArray O;                                            /* ~ temporary storage for vorticity             ~ */
  O.reserve(usize);
                                                             /* ~ note: P is already known                    ~ */

  for (unsigned k = 0; k <usize; k++) {O[k] = U0[k];}        /* ~ not necessary. Could use U0 directly        ~ */

  unsigned idx         = 0;                                  /* ~ index for inv_k2                            ~ */

  for (unsigned k = 0; k < usize; k++) {

    if ( k % n1n2c == 0) { idx = 0; }                        /* ~ reset idx when starting new layer           ~ */
    P[k] = inv_k2[idx] * O[k];                               /* ~ Omega = - delperp^2 P                       ~ */
                                                             /* ~ NOTE: why not use known value?              ~ */
    ++idx;

  }

  for (unsigned k = 0; k <usize; k++) {U0[k] = P[k];}        /* ~ U0 now holds Fourier transform of P         ~ */

//  O.resize(0);                                             /* ~ vorticity is discarded here...              ~ */
                                                             /* ~ ... or perhaps not.                         ~ */
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

  ComplexArray H;                                         /* ~ temporary storage for H-function           ~ */
  H.reserve(usize);
                                                          /* ~ note: current values of A and J are known  ~ */

  RealVar ssqd;      run.palette.fetch("ssqd",   &ssqd ); /* ~ parameter sigma^2 needed for A             ~ */
  std::string model; run.palette.fetch( "model", &model);

  for (unsigned k = 0; k < usize; k++) {H[k] = U1[k];}

  unsigned idx     = 0;                                   /* ~ index for k2                               ~ */
  for (unsigned k  = 0; k < usize; k++) {

    if (k % n1n2c  == 0) { idx = 0; }                     /* ~ reset idx when starting new layer          ~ */

    if (model.compare("hall") == 0 ) {
      A[k]         = H[k] / (one + ssqd*k2[idx]);         /* ~ NOTE: why not use known values?            ~ */
    }
    else if(model.compare("rmhd") == 0 || model.compare("inhm") == 0 ) {
      A[k]         = H[k];                                /* ~ NOTE: why not use known values?            ~ */
    }

    ++idx;

  }

  for (unsigned k = 0; k < usize; k++) {U1[k] = A[k];}    /* ~  U1 now holds Fourier transform of A       ~ */

//  H.resize(0);                                          /* ~ H is discarded...                          ~ */
                                                          /* ~ ...or perhaps not.                         ~ */
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalElls(    stack& run ) {

  std::string nprofile; run.palette.fetch("nprofile", &nprofile);
  int p3;               run.palette.fetch("p3",       &p3      );

  Elln.assign(p3+1,zero);
  EllA.assign(p3+1,zero);
  EllB.assign(p3+1,zero);

   h11.assign(p3+1,zero);
   h12.assign(p3+1,zero);
   h21.assign(p3+1,zero);
   h22.assign(p3+1,zero);

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
  case(1) : for ( int k = 0; k <  p3+1;  ++k) {

              EllA[k] = abs( dvalfdz[k] / (two  * valfven[k]) );
              Elln[k] = abs( dndz[k]    / (four * nofz[k]   ) );

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

  for ( int k = 0; k <  p3+1;  ++k) {

    h11[k] =  umean[k]   * ( EllA[k] + EllB[k] - Elln[k] );
    h12[k] = -valfven[k] * ( EllA[k] + EllB[k] + Elln[k] );
    h21[k] =  h12[k];
    h22[k] =  umean[k]   * ( EllB[k] - Elln[k] - EllA[k] );

  }

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalValf( stack& run ) {
  
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

              std::cout << "evalValf: H0      = " << H0      << std::endl;    /* ~ maybe add this to physics_data? ~ */
//            std::cout << "evalValf: valfmin = " << valfmin << std::endl;    /* ~ maybe add this to physics_data? ~ */

                 nofz.assign(p3+1,one );
                 dndz.assign(p3+1,zero);
              valfven.assign(p3+1,one );
              dvalfdz.assign(p3+1,zero);

              for (unsigned k = 0; k < p3+1; ++k) {

                 nofz[k]    = ninf + ((n0 - ninf) * (exp(-half*pow((z[k] - (half*zl))/H0 ,two))));
                 dndz[k]    = -(n0 - ninf) * exp( -half * pow(((z[k]-(half*zl))/H0),2)) * (z[k] - (half*zl)) / ( pow(H0,2) );
                 valfven[k] = valfmin * sqrt( n0 / nofz[k] );
                 dvalfdz[k] = -half * valfven[k] * dndz[k] / nofz[k];

              }
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

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//  for (unsigned k = 0; k < p3+1; ++k) {
//    std::cout << "evalValf:  rank = " << rank << ": valfven[" << k << "] =  " << valfven[k] << std::endl;
//  }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalUmean( stack& run ) {
  
  std::string uprofile; run.palette.fetch("uprofile", &uprofile);
  int         p3;       run.palette.fetch("p3",       &p3      );

  int        i_profile;

  if      (uprofile.compare("noflow")   == 0) { i_profile =  0; }
  else if (uprofile.compare("uniform")  == 0) { i_profile =  1; }
  else if (uprofile.compare("whoknows") == 0) { i_profile =  2; }
  else                                        { i_profile = -1; }

  switch(i_profile) {

  case(0) : umean.assign(p3+1, zero); break;
  case(1) : umean.assign(p3+1, zero); break;
  case(2) : umean.assign(p3+1, zero); break;

  default : 

    std::cout << "evalValf: WARNING - the profile " << uprofile << " is not implemented. Assuming a flat profile." << std::endl;
    umean.assign(p3+1, one);

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::trackEnergies(double t_cur, stack& run ) {


  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  RealArray& EnergyQs    = run.EnergyQs;

  static const int i_tcr =  0; double tcr; /* ~ Current time                                    ~ */
  static const int i_pe  =  1; double pe ; /* ~ Total perpendicular kinetic energy              ~ */
  static const int i_ae  =  2; double ae ; /* ~ Total perpendicular perturbed magnetic energy   ~ */
  static const int i_mo  =  3; double mo ; /* ~ Maximum vorticity                               ~ */
  static const int i_imo =  4; double imo; /* ~ i index location of maximum vorticity           ~ */
  static const int i_jmo =  5; double jmo; /* ~ j index location of maximum vorticity           ~ */
  static const int i_oe  =  6; double oe ; /* ~ Total square magnitude of vorticity             ~ */
  static const int i_mj  =  7; double mj ; /* ~ Maximum current                                 ~ */
  static const int i_imj =  8; double imj; /* ~ i index location of maximum current             ~ */
  static const int i_jmj =  9; double jmj; /* ~ j index location of maximum current             ~ */
  static const int i_ce  = 10; double ce ; /* ~ Total square magnitude of current               ~ */
  static const int i_noe = 11; double noe; /* ~ viscous dissipation                             ~ */
  static const int i_ece = 12; double ece; /* ~ resistive dissipation                           ~ */
  static const int i_fe  = 13; double fe ; /* ~ Poynting Flux                                   ~ */
  static const int i_ftp = 14; double ftp; /* ~ Time-averaged Poynting Flux                     ~ */
  static const int i_eds = 15; double eds; /* ~ Time average of total dissipated energy         ~ */
  static const int i_dng = 16; double dng; /* ~ Average rate-of-change of total energy          ~ */
  static const int i_irc = 17; double irc; /* ~ Instantaneous rate-of-change of total energy    ~ */
  static const int i_gml = 18; double gml; /* ~ time-average of energy gains minus losses       ~ */
  static const int i_cns = 19; double cns; /* ~ Check on energy conservation                    ~ */
  static const int i_ttc = 20; double ttc; /* ~ Current time in units of correlation time       ~ */
  static const int i_dt  = 21; double dt ; /* ~ current value of time increment                 ~ */
  static const int i_dtv = 22; double dtv; /* ~ time increment adjustment parameter             ~ */
  static const int i_coc = 23; double coc; /* ~ root mean square J^2 / Del J^2                  ~ */
  static const int i_vkt = 24; double vkt; /* ~ average rate-of-change of j^2 / magnetic energy ~ */ 
  static const int i_avm = 25; double avm; /* ~ Time averaged magnetic field strength           ~ */
  static const int i_avp = 26; double avp; /* ~ Time-averaged footpoint velocity                ~ */
  static const int i_fp  = 27; double fp ; /* ~ Total footpoint kinetic energy at z = 0         ~ */                  
  static const int i_he  = 28; double he ; /* ~ Total energy in helicity due to inhomogeneity   ~ */

  double dtvb; physics_data.fetch("dtvb", &dtvb);
  double tauC; run.palette.fetch( "tauC", &tauC);
  double eta;  run.palette.fetch( "eta",  &eta );
  double nu;   run.palette.fetch( "nu",   &nu  );

  double cee;
  double t_old;
  double dt_old;
  double aeold;
  double peold;
  double heold;


  pe                     = evalTotalKineticEnergy(  run );
  ae                     = evalTotalMagneticEnergy( run );
  oe                     = evalTotalVorticitySqd(   run );
  ce                     = evalTotalCurrentSqd(     run );
  cee                    = evalTotalGradCurrentSqd( run );
  fp                     = evalTotalFootPointKE(    run );
  fe                     = evalTotalPoyntingFlux(   run );

  std::string model; run.palette.fetch("model", &model);

  if (     model.compare("rmhd") == 0) { he = zero;                           }
  else if (model.compare("inhm") == 0) { he = evalTotalHelicalEnergy(  run ); }

//  if (rank == 0) {
//    std::cout << "trackEnergies: he = " << he << std::endl;
//  }

  if ( EnergyQs.size()  >= i_he ) {

    t_old                = EnergyQs[i_tcr];
    dt_old               = EnergyQs[i_dt];
 
    aeold                = EnergyQs[i_ae];
    peold                = EnergyQs[i_pe];
    heold                = EnergyQs[i_he];

    tcr                  = t_cur;

/* ~ pe evaluated above ~ */
/* ~ ae evaluated above ~ */

    mo                   = zero;
    imo                  = zero;
    jmo                  = zero;

/* ~ oe evaluated above ~ */

    mj                   = zero;
    imj                  = zero;
    jmj                  = zero;

/* ~ ce evaluated above ~ */

    if (rank == 0 ) {

      RealVar AVEz;  run.palette.fetch("AVEz",  &AVEz );
      RealVar AVEpv; run.palette.fetch("AVEpv", &AVEpv);
      RealVar AVEpn; run.palette.fetch("AVEpn", &AVEpn);
      RealVar AVEpe; run.palette.fetch("AVEpe", &AVEpe);

      noe                = nu  * oe;
      ece                = eta * ce;

/* ~ fe evaluated above ~ */

      ftp                = ((EnergyQs[i_ftp]*t_old) + fe          * dt_old)        / t_cur;
      eds                = ((EnergyQs[i_eds]*t_old) + (noe + ece) * dt_old)        / t_cur;
      dng                = ( (EnergyQs[i_dng] * t_old) + (ae - aeold + pe - peold + he - heold) ) / t_cur;
      irc                = (ae - aeold + pe - peold + he - heold) / dt_old;
      gml                = ftp - eds;
      cns                = abs(( (fe - noe - ece) - ( (ae - aeold + pe - peold + he - heold ) / dt_old )) * dt_old );

      AVEz               = AVEz  + cns * dt_old;
      AVEpv              = AVEpv + noe * dt_old;
      AVEpn              = AVEpn + ece * dt_old;
      AVEpe              = AVEpe + fp  * dt_old;

      run.palette.reset("AVEz",  AVEz );
      run.palette.reset("AVEpv", AVEpv);
      run.palette.reset("AVEpn", AVEpn);
      run.palette.reset("AVEpe", AVEpe);

      ttc                = t_cur / tauC; 
      run.palette.fetch("dt",&dt);
      dtv                = dtvb;
      if (cee == zero){
        coc              = zero; 
      }
      else {
        coc              = sqrt(ce/cee);
      }
      if (ae != zero) {

        vkt              = EnergyQs[i_vkt]  * t_old;
        vkt              = (vkt +  (ce/ae)  * dt_old)     / t_cur;

      }
      else { vkt         = zero; }

      avm                = pow(EnergyQs[i_avm],2) * t_old;
      avm                = sqrt((avm + ae * dt_old)       / t_cur);
//    avp                = half * pow(EnergyQs[i_avp],2)  * t_old;
//    avp                = sqrt(two * (avp + fp * dt_old) / t_cur);
      avp                = sqrt(two * (AVEpe + fp * dt_old) / t_cur);

      EnergyQs[ i_tcr ]  = tcr;
      EnergyQs[ i_pe  ]  = pe ;
      EnergyQs[ i_ae  ]  = ae ;
      EnergyQs[ i_mo  ]  = mo ;
      EnergyQs[ i_imo ]  = imo;
      EnergyQs[ i_jmo ]  = jmo;
      EnergyQs[ i_oe  ]  = oe ;
      EnergyQs[ i_mj  ]  = mj ;
      EnergyQs[ i_imj ]  = imj;
      EnergyQs[ i_jmj ]  = jmj;
      EnergyQs[ i_ce  ]  = ce ;
      EnergyQs[ i_noe ]  = noe;      
      EnergyQs[ i_ece ]  = ece;
      EnergyQs[ i_fe  ]  = fe ;
      EnergyQs[ i_ftp ]  = ftp;
      EnergyQs[ i_eds ]  = eds;
      EnergyQs[ i_dng ]  = dng;
      EnergyQs[ i_irc ]  = irc;
      EnergyQs[ i_gml ]  = gml;
      EnergyQs[ i_cns ]  = cns;
      EnergyQs[ i_ttc ]  = ttc;
      EnergyQs[ i_dt  ]  = dt ;
      EnergyQs[ i_dtv ]  = dtv;
      EnergyQs[ i_coc ]  = coc;
      EnergyQs[ i_vkt ]  = vkt;
      EnergyQs[ i_avm ]  = avm;
      EnergyQs[ i_avp ]  = avp;
      EnergyQs[ i_fp  ]  =  fp;
      EnergyQs[ i_he  ]  =  he;

    }
  }

  else {

    int srun; run.palette.fetch("srun", &srun);

    if (srun == 1) {

      run.palette.fetch("dt", &dt_old);
      t_old                = t_cur;

      aeold                = zero;
      peold                = zero;
      heold                = zero;

      tcr                  = t_cur;

/* ~   pe evaluated above ~ */
/* ~   ae evaluated above ~ */

      mo                   = zero;
      imo                  = zero;
      jmo                  = zero;

/* ~   oe evaluated above ~ */

      mj                   = zero;
      imj                  = zero;
      jmj                  = zero;

/* ~   ce evaluated above ~ */

      if (rank == 0 ){

        noe                = nu  * oe;
        ece                = eta * ce;

/* ~   fe evaluated above ~ */

        ftp                = zero;
        eds                = zero;
        dng                = zero;
        irc                = (ae - aeold + pe - peold) / dt_old;
        gml                = (ae - aeold + pe - peold )/ dt_old;
        cns                = abs(  ( ( fe - noe - ece ) - ((ae - aeold + pe - peold) / dt_old)) * dt_old );
        ttc                = t_cur / tauC; 
        dt                 = dt_old;
        dtv                = dtvb;
        if (cee == zero){
          coc              = zero; 
        }
        else {
          coc              = sqrt(ce/cee);
        }
        vkt                = zero;
        avm                = zero;
        avp                = zero;

        EnergyQs.push_back(tcr);
        EnergyQs.push_back(pe );
        EnergyQs.push_back(ae );
        EnergyQs.push_back(mo );
        EnergyQs.push_back(imo);
        EnergyQs.push_back(jmo);
        EnergyQs.push_back(oe );
        EnergyQs.push_back(mj );
        EnergyQs.push_back(imj);
        EnergyQs.push_back(jmj);
        EnergyQs.push_back(ce );
        EnergyQs.push_back(noe);
        EnergyQs.push_back(ece);
        EnergyQs.push_back(fe );
        EnergyQs.push_back(ftp);
        EnergyQs.push_back(eds);
        EnergyQs.push_back(dng);
        EnergyQs.push_back(irc);
        EnergyQs.push_back(gml);
        EnergyQs.push_back(cns);
        EnergyQs.push_back(ttc);
        EnergyQs.push_back(dt );
        EnergyQs.push_back(dtv);
        EnergyQs.push_back(coc);
        EnergyQs.push_back(vkt);
        EnergyQs.push_back(avm);
        EnergyQs.push_back(avp);
        EnergyQs.push_back(fp );
        EnergyQs.push_back(he );

      }
    }
    else {
      if (rank == 0) {

        std::string prefix;    run.palette.fetch(   "prefix",    &prefix   );
        std::string run_label; run.palette.fetch(   "run_label", &run_label);
        std::string res_str;   run.stack_data.fetch("res_str",   &res_str  );

        std::string energy_data_file = prefix + '_' + res_str + ".o" + run_label;
        const char *c_data_file      = energy_data_file.c_str();

        std::ifstream ifs;
        ifs.open( c_data_file, std::ios::in );

        if (ifs.good()) {

          unsigned esize                = 29;

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

        } // ifs is good

        else { std::cout << "trackEnergyQs: Warning - could not open file " << energy_data_file << std::endl; }

      } // rank is zero
    } // srun is not 1
  } // i_fp is zero
} // end function

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::trackPowerSpectra(RealVar t_cur, stack& run ) {

  int calcsvz; run.palette.fetch("calcsvz", &calcsvz);
  if (calcsvz == 1) {

     int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     int   n3; run.palette.fetch("p3",   &n3    );

     static const int i_k     = 0;
     static const int i_spe   = 1;
     static const int i_sae   = 2;
     static const int i_ts    = 3;
     static const int i_szp   = 4;
     static const int i_szm   = 5;
     static const int i_tz    = 6;

     int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c);

     RealArray& k2            = run.k2;
     int ksize                = k2.capacity();

     RealArray spe;
     RealArray sae;
     RealArray szp;
     RealArray szm;


     assert(n1n2c == ksize);

     
     for (unsigned l = 1; l < n3 + 1; ++l ) {

       spe.assign(isp+1,zero);
       sae.assign(isp+1,zero);
       szp.assign(isp+1,zero);
       szm.assign(isp+1,zero);

       int m;
       int idx;
       for (unsigned k = 0; k < n1n2c; ++k) {

         m       = 1 + sqrt(k2[k]) * dk_m1;
         idx     = (l*n1n2c) + k;

         spe[m]  = spe[m] + k2[k] * ( pow(std::norm(P[idx]),         2) );
         sae[m]  = sae[m] + k2[k] * ( pow(std::norm(A[idx]),         2) );
         szp[m]  = szp[m] + k2[k] * ( pow(std::norm(P[idx] + A[idx]),2) );
         szm[m]  = szm[m] + k2[k] * ( pow(std::norm(P[idx] - A[idx]),2) );

       }

       for (unsigned j = 0; j < isp+1; ++j) {

         spe[j] = two * spe[j] * dk_m1;
         sae[j] = two * sae[j] * dk_m1;
         szp[j] = two * szp[j] * dk_m1;
         szm[j] = two * szm[j] * dk_m1;

       }

//     for (unsigned j = 0; j < isp+1; ++j) {
//
//       SpcVsZ[j][0][0] = j * dk;
//       if (rank == 0) { std::cout << "j = " << j << "jmax = " << isp+1 << std::endl; }
//         SpcVsZ[j][i_k  ][l] = j * dk;
//         SpcVsZ[j][i_spe][l] = SpcVsZ[j][i_spe][l] + spe[j];
//         SpcVsZ[j][i_sae][l] = SpcVsZ[j][i_sae][l] + sae[j];
//         SpcVsZ[j][i_ts ][l] = SpcVsZ[j][i_ts ][l] + spe[j] + sae[j];
//         SpcVsZ[j][i_szp][l] = SpcVsZ[j][i_szp][l] + szp[j];
//         SpcVsZ[j][i_szm][l] = SpcVsZ[j][i_szm][l] + szm[j];
//         SpcVsZ[j][i_tz ][l] = SpcVsZ[j][i_tz ][l] + szp[j] + szm[j];
//
//     }

     }

     // need dk and dk_m1
     //  isp -> defaults to n2/2    
     //  kb  -> defaults to 2 pi
     //  dk  -> defaults (pi *n1) / isp
     //  kf  -> defaults to kb + isp * dk 
     //  nk  -> given by krange
     //  ke  -> given by krange
     //  ikb -> given by krange also (?) defaults to kb / dk
     //  ikf -> given by krange

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::trackQtyVsZ(RealVar t_cur, stack& run ) {

  int calcqvz; run.palette.fetch("calcqvz", &calcqvz);

  if (calcqvz == 1) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int   n3; run.palette.fetch("p3",   &n3    );

    static const int i_z   =  0;
    static const int i_ae  =  1;
    static const int i_pe  =  2;
    static const int i_ch  =  3;
    static const int i_ep  =  4;
    static const int i_em  =  5;
    static const int i_ce  =  6;
    static const int i_oe  =  7;
    static const int i_zp  =  8;
    static const int i_zm  =  9;
    static const int i_nch = 10;
    static const int i_km  = 11;
    static const int i_kp  = 12;
    static const int i_kzp = 13;
    static const int i_kzm = 14;

    RealVar dz; run.stack_data.fetch("dz",&dz    );
    RealArray& z           = run.z;

    RealArray aevsz;
    RealArray pevsz;
    RealArray chvsz;
    RealArray cevsz;
    RealArray oevsz;
    RealArray zpvsz;
    RealArray zmvsz;
    RealArray zepvsz;
    RealArray zemvsz;

    RealArray nchvsz;
    RealArray kcvsz;
    RealArray kovsz;
    RealArray kzpvsz;
    RealArray kzmvsz;

    aevsz.assign( n3,zero);
    pevsz.assign( n3,zero);
    chvsz.assign( n3,zero);
    cevsz.assign( n3,zero);
    oevsz.assign( n3,zero);
    zpvsz.assign( n3,zero);
    zmvsz.assign( n3,zero);
    zepvsz.assign(n3,zero);
    zemvsz.assign(n3,zero);
    nchvsz.assign(n3,zero);
    kcvsz.assign( n3,zero);
    kovsz.assign( n3,zero);
    kzpvsz.assign(n3,zero);
    kzmvsz.assign(n3,zero);

    int kdx;
    RealArray& k2 = run.k2;
    int n1n2c     = k2.capacity();

    for (unsigned l = 0; l < n3; l++) {

      kdx         = 0;

      for ( unsigned k = ((l+1) * n1n2c) ; k < ((l+2) * n1n2c); k++ ) {

        aevsz[l]  = dz*(aevsz[l] + k2[kdx] * pow(std::norm(A[k]),2));             /* ~ Magnetic Energy layer l  ~ */
        pevsz[l]  = dz*(pevsz[l] + k2[kdx] * pow(std::norm(P[k]),2));             /* ~ Kinetic  Energy layer l  ~ */

        chvsz[l]  = dz*(chvsz[l] + k2[kdx] * (P[k].real() * A[k].real()) + (P[k].imag() * A[k].imag())) * two; 
        cevsz[l]  = dz*(cevsz[l] + k2[kdx] * k2[kdx] * pow(std::norm(A[k]),2) );
        oevsz[l]  = dz*(oevsz[l] + k2[kdx] * k2[kdx] * pow(std::norm(P[k]),2) );
        zpvsz[l]  = dz*(zpvsz[l] + k2[kdx] * pow(std::norm((A[k] + P[k])),2 ) );
        zmvsz[l]  = dz*(zmvsz[l] + k2[kdx] * pow(std::norm((A[k] - P[k])),2 ) );

        ++kdx;

      }

      zepvsz[l] = aevsz[l] + pevsz[l] + chvsz[l];
      zemvsz[l] = aevsz[l] + pevsz[l] - chvsz[l];

      if ( (abs(aevsz[l]) >= teensy) || (abs(pevsz[l])  >= teensy) ) {
        nchvsz[l] = chvsz[l] / (aevsz[l] + pevsz[l]) ;
      }
      else { nchvsz[l] = zero; }
      if ( abs(aevsz[l]) >= teensy ) {
         kcvsz[l] = cevsz[l] / aevsz[l];
         if ( (abs(kcvsz[l]) < teensy) || (abs(kcvsz[l]) > huge)) { kcvsz[l] = zero; }
      }
      else{ kcvsz[l] = zero;}
      if ( abs(pevsz[l]) >= teensy) {
        kovsz[l]  = oevsz[l] / pevsz[l];
        if ((abs(kovsz[l]) < teensy) || (abs(kovsz[l]) > huge)) { kovsz[l] = zero; }
      }
      else {kovsz[l] = zero ;}
      if (abs(zepvsz[l]) >= teensy) {
        kzpvsz[l] = zpvsz[l] / zepvsz[l];
        if ((abs(kzpvsz[l]) < teensy) || (abs(kzpvsz[l]) > huge) ) {kzpvsz[l] = zero; }
      }
      else{kzpvsz[l] =zero;}
      if (abs(zemvsz[l]) >= teensy) {
        kzmvsz[l] = zmvsz[l] / zemvsz[l];
        if ((abs(kzmvsz[l]) < teensy) || (abs(kzmvsz[l]) > huge) ) {kzmvsz[l] = zero; }
      }
      else{ kzmvsz[l] =zero; }

      QtyVsZ[(rank*n3) + l][ 0]   = z[l + 1];
      QtyVsZ[(rank*n3) + l][ 1]   = QtyVsZ[(rank*n3) + l][ 1] + aevsz[ l];
      QtyVsZ[(rank*n3) + l][ 2]   = QtyVsZ[(rank*n3) + l][ 2] + pevsz[ l];
      QtyVsZ[(rank*n3) + l][ 3]   = QtyVsZ[(rank*n3) + l][ 3] + chvsz[ l];
      QtyVsZ[(rank*n3) + l][ 4]   = QtyVsZ[(rank*n3) + l][ 4] + cevsz[ l];
      QtyVsZ[(rank*n3) + l][ 5]   = QtyVsZ[(rank*n3) + l][ 5] + oevsz[ l];
      QtyVsZ[(rank*n3) + l][ 6]   = QtyVsZ[(rank*n3) + l][ 6] + zpvsz[ l];
      QtyVsZ[(rank*n3) + l][ 7]   = QtyVsZ[(rank*n3) + l][ 7] + zmvsz[ l];
      QtyVsZ[(rank*n3) + l][ 8]   = QtyVsZ[(rank*n3) + l][ 8] + zepvsz[l];
      QtyVsZ[(rank*n3) + l][ 9]   = QtyVsZ[(rank*n3) + l][ 9] + zemvsz[l];
      QtyVsZ[(rank*n3) + l][10]   = QtyVsZ[(rank*n3) + l][10] + nchvsz[l];
      QtyVsZ[(rank*n3) + l][11]   = QtyVsZ[(rank*n3) + l][11] + kcvsz[ l];
      QtyVsZ[(rank*n3) + l][12]   = QtyVsZ[(rank*n3) + l][12] + kovsz[ l];
      QtyVsZ[(rank*n3) + l][13]   = QtyVsZ[(rank*n3) + l][13] + kzpvsz[l];
      QtyVsZ[(rank*n3) + l][14]   = QtyVsZ[(rank*n3) + l][14] + kzmvsz[l];

    }
      //if (rank == 0) {
      //  for (int m = 0; m < n3; ++m) {
      //    std::cout << "trackQ: z[" << rank*n3 + m  << "] = " << QtyVsZ[(rank*n3)+m][0] << std::endl;
      //  }
      //}
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::reportEnergyQs ( stack& run ) {

  int calcqvz; run.palette.fetch("calcqvz", &calcqvz);

  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank  );

  if (calcqvz == 1) {
     if (rank == 0) {

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

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::reportQtyVsZ ( RealVar t_cur, stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  MPI_Status * status = 0;

  int   nw; run.palette.fetch("nw",   &nw    );
  int   n3; run.palette.fetch("p3",   &n3    );
  int   np; run.palette.fetch("np",   &np    );

  const char *c_qty_vs_z_out_file;
  std::ofstream ofs;

  RealVar nw_m1 = one / (RealVar (nw));

 /*  ~ It should be possible to do the following with an MPI_Gather call, but for the life of me
  *  ~ I can't make it work.
  */

 for (unsigned l = 0; l < n3; l++) {

   if (rank != 0 ) { MPI_Send(&QtyVsZ[rank*n3 + l].front(), 15, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD); }
   else { for (unsigned rnk_k = 1; rnk_k < np; rnk_k++) {
            MPI_Recv(&QtyVsZ[rnk_k*n3+l].front(), 15, MPI_DOUBLE, rnk_k, rnk_k, MPI_COMM_WORLD, status);
          }
        }
 }

  if (rank == 0 ) {

    std::string qout_prefix;  run.palette.fetch(   "qout_pref", &qout_prefix);
    std::string res_str;      run.stack_data.fetch("res_str",   &res_str    );
    std::string run_label;    run.palette.fetch(   "run_label", &run_label  );
    int srun;                 run.palette.fetch(   "srun",      &srun       );
    std::string srn_str = static_cast<std::ostringstream*>( &(std::ostringstream() << srun) ) -> str();

    std::string qty_vs_z_out_file = qout_prefix + "_" + res_str + ".o" + run_label + srn_str;
    c_qty_vs_z_out_file = qty_vs_z_out_file.c_str();

    std::cout << "reportQtyVsZ: qty_vs_z_out_file = " << qty_vs_z_out_file << std::endl;

    ofs.open( c_qty_vs_z_out_file, std::ios::out );

    if (ofs.good() ) {
      ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << t_cur << std::endl;
      for (unsigned l = 0; l < np*n3; l++) {
        for (unsigned k = 0; k < 15; k++) { 
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

      pe            = pe +        k2[idx] * pow(abs(P[kdx]), 2);

      if (( rank    == np - 1 ) && ( kdx >= (four * n1n2c) )) {
        pe          = pe + half * k2[idx] * pow(abs(P[kdx]), 2);
      }
      if ((rank     == 0      ) && ( kdx <   n1n2c         )) {
        pe          = pe - half * k2[idx] * pow(abs(P[kdx]), 2);
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

      oe        = oe +        k2[idx] * k2[idx] * pow(abs(P[kdx]), 2);

      if (( rank == np - 1 ) && ( kdx >= (four * n1n2c) )) {
        oe      = oe + half * k2[idx] * k2[idx] * pow(abs(P[kdx]), 2);
      }
      if ((rank  == 0      ) && ( kdx <   n1n2c)           ) {
        oe      = oe - half * k2[idx] * k2[idx] * pow(abs(P[kdx]), 2);
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

      me              = me + k2[idx] * pow(abs(A[kdx]), 2);
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

      ce              = ce + k2[idx] * k2[idx] * pow(abs(A[kdx]), 2);
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

      cee             = cee + k2[idx] * k2[idx] * k2[idx] * pow(abs(A[kdx]), 2);
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

//    fp              = fp + k2[kdx] * pow(abs(P[kdx]), 2);
      fp              = fp + pow(abs(O[kdx]), 2);
    }

// or
//  fp                = three * fp * dz;
    fp                = fp * dz;

//  fp                = two_thirds * fp;  /* ~ NOTE!:  do I need this? ~ */

  }

  return fp;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalPoyntingFlux ( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  RealArray& k2       = run.k2;

  int np;    run.palette.fetch(   "np",    &np    );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;   run.stack_data.fetch("iu2",   &iu2   );
  int n3;    run.stack_data.fetch("n3",    &n3    );
  double dz; run.stack_data.fetch("dz",    &dz    );
  double n0; run.palette.fetch("n0",       &n0    ); /* ~ density at z_0 (i.e. z = zl / 2)     ~ */

  std::string model; run.palette.fetch("model",      &model);

  double fe           = zero;
  double dfe_r        = zero;
  double dfe_i        = zero;
  double fe_sum       = zero;

  double Valf;
  double valfmax;       run.palette.fetch("valfmax",  &valfmax );

  if (      model.compare("rmhd") == 0) { Valf = one;              }
  else if ( model.compare("inhm") == 0) {
    if (      rank == 0     ) {           Valf = valfmax;          }
    else if ( rank == np -1 ) {           Valf = valfven[ n3 ];    }
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

        dfe_r         =       k2[idx] * ( (A[kdx].imag() * P[kdx].real()) + (A[kdx].real() * P[kdx].imag()) );
        dfe_i         =       k2[idx] * ( (A[kdx].imag() * P[kdx].imag()) - (A[kdx].real() * P[kdx].real()) );

        fe            = fe + sqrt(pow(dfe_r,2) + pow(dfe_i,2));

      }
      if ((rank  == 0      ) && ( kdx <   n1n2c)           ) {

        dfe_r         =       k2[idx] * ( (A[idx + n1n2c].imag() * P[idx].real()) + (A[idx + n1n2c].real() * P[kdx].imag()) );
        dfe_i         =       k2[idx] * ( (A[idx + n1n2c].imag() * P[idx].imag()) - (A[idx + n1n2c].real() * P[kdx].real()) );

        fe            = fe - sqrt(pow(dfe_r,2) + pow(dfe_i,2));

      }
      ++idx;
  }

  fe                  = two * Valf * fe;

  int i_red           = MPI_Reduce(&fe, &fe_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

// fe_sum = two*two_thirds * fe_sum;  /* ~ NOTE!: why the two_thirds. Seems necessary but at moment I don't know why ~ */
// fe_sum = fe_sum;  /* ~ NOTE!: why the two_thirds. Seems necessary but at moment I don't know why ~ */

  return fe_sum;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalHelicalEnergy ( stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  RealArray& k2       = run.k2;

  int np;    run.palette.fetch(   "np",    &np    );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;   run.stack_data.fetch("iu2",   &iu2   );
  int n3;    run.stack_data.fetch("n3",    &n3    );
  double dz; run.stack_data.fetch("dz",    &dz    );

  double he           = zero;
  double dhe_r        = zero;
  double dhe_i        = zero;
  double he_sum       = zero;

  int idx             = 0;
  int ldx             = 0;

  int kstart          = n1n2c;
  int kstop           = (n3+1)*n1n2c;
  RealVar lcoff;

  for (unsigned kdx = kstart; kdx < kstop; kdx++){

        if ( idx % n1n2c == 0) {
          ++ldx;
          std::cout << "" << std::endl;
          std::cout << "evalTotalH: for rank = " << rank << ": Elln["    << ldx  <<   "]    = " << Elln[ldx]    << std::endl;
          std::cout << "evalTotalH: for rank = " << rank << ": EllB["    << ldx  <<   "]    = " << EllB[ldx]    << std::endl;
          std::cout << "evalTotalH: for rank = " << rank << ": EllA["    << ldx  <<   "]    = " << EllA[ldx]    << std::endl;
          std::cout << "evalTotalH: for rank = " << rank << ": valfven[" << ldx  <<   "] = "    << valfven[ldx] << std::endl;
          std::cout << "evalTotalH: for rank = " << rank << ": dvalfdz[" << ldx  <<   "] = "    << dvalfdz[ldx] << std::endl;

          idx         = 0 ; 

          lcoff       = (dvalfdz[ldx] + two * valfven[ldx] * (Elln[ldx] + EllB[ldx] + EllA[ldx] ));

        }

        dhe_r         =       k2[idx] * ( (A[kdx].imag() * P[kdx].real()) + (A[kdx].real() * P[kdx].imag()) );
        dhe_i         =       k2[idx] * ( (A[kdx].imag() * P[kdx].imag()) - (A[kdx].real() * P[kdx].real()) );

        he            = he + lcoff * sqrt(pow(dhe_r,2) + pow(dhe_i, 2));

      ++idx;
  }
  
  he                  = he * dz;

  int i_red           = MPI_Reduce(&he, &he_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return he_sum;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyBC( std::string str_step, stack& run ) {

  int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bdrys; run.palette.fetch("bdrys", &bdrys);
  int np;    run.palette.fetch("np"   , &np   );

  if ( rank == 0 || rank == np - 1) {

    if (bdrys > 0 ) {
      if (!str_step.compare("finalize") == 0) { applyFootPointDrivingBC( str_step, run ); }
    }
    else            { applyLineTiedBC( str_step, run );                           }

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyFootPointDrivingBC( std::string str_step, stack& run ) {

  int rank;       MPI_Comm_rank(MPI_COMM_WORLD, &rank      );

  int bdrys;      run.palette.fetch(    "bdrys", &bdrys    );
  int np;         run.palette.fetch(    "np"   , &np       );
  int n1n2c;      run.stack_data.fetch( "n1n2c", &n1n2c    );
  int n_layers;   run.stack_data.fetch( "n3"   , &n_layers );
  RealVar dtau;   run.palette.fetch(    "dtau" , &dtau     );
  RealVar t_cur;  physics_data.fetch(   "t_cur", &t_cur    );

  int oldnum, num;

  RealVar    ffp, kc;
  RealVar    dummy;
  RealVar    next_real, next_imag;
  ComplexVar tuple;

  RealVar    lowtau; 
  RealVar    bigtau; 

  RealVar    a, b;                                /* ~ i.e. Gilson's "a" and "b" ~ */

  ComplexArray::size_type nc = n1n2c;

  int iu3;
  run.stack_data.fetch("iu3", &iu3);

  unsigned strt_idx, stop_idx;

  if ( rank  == 0 || rank == (np - 1) ) {         /* ~ pevol "starts" here        ~ */

    RealArray&    k2         = run.k2;

    ComplexArray& O          = run.U0;
    ComplexArray& Z          = run.U2; 

    if ( rank  == 0 ) {

      run.palette.fetch("oldnumlb", &oldnum );    /* ~ "numold"                   ~ */
      strt_idx               = 0;
      stop_idx               = strt_idx + n1n2c;

      for (unsigned k = strt_idx; k < stop_idx; k++) {
        O[k]                 = czero;
        if (iu3 > 2){ 
        Z[k]                 = czero;
        }
      }
      num                    = oldnum;
     
      while (  (num * dtau) <= t_cur ) { ++num; }

      num                    = num - 1; 
      lowtau                 = (num    * dtau) - t_cur;
      bigtau                 =((num+1) * dtau) - t_cur;
      a                      = cos((pi*lowtau)/( two * dtau)); /* ~ Gilson's "interp" ~ */
      b                      = cos((pi*bigtau)/( two * dtau));

      if (num == oldnum) {
    
        for (unsigned k = strt_idx; k < stop_idx; k++) {
          O[k]               = (a * roldlb[k]) + (b * rnewlb[k]);
          if (iu3 > 2){ 
            Z[k]             = czero;
          }
        }
      }
      else {

        int brcount; physics_data.fetch("brcount", &brcount);

        run.palette.fetch("ffp", &ffp);
        run.palette.fetch( "kc", &kc );

        for (unsigned k = 0; k < num - oldnum; k++) {

          for (unsigned l = 0; l < nc; l++) { roldlb[l] = rnewlb[l]; }

          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) <= kc ) {

                dummy        = ((double) rand() / RAND_MAX );
                ++brcount;
                dummy        = ((double) rand() / RAND_MAX );
                ++brcount;
            }
          }

          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) <= kc ) {

                next_real    = ffp * (((double) rand() / RAND_MAX) * two - one);
                ++brcount;
                next_imag    = ffp * (((double) rand() / RAND_MAX) * two - one);
                ++brcount;
                tuple        = ComplexVar(next_real, next_imag);
                rnewlb[l]    = tuple;

            }
          }
        }
        for (unsigned k = strt_idx; k < stop_idx; k++) {
          O[k]               = (a * roldlb[k]) + (b * rnewlb[k]);
          if (iu3 > 2){ 
            Z[k]             = czero;
          }
        }
        physics_data.reset("brcount", brcount);
      }
      run.palette.reset("oldnumlb", num);
    }
    else if ( rank == np - 1 && bdrys == 2) {

      run.palette.fetch("oldnumub", &oldnum ); /* ~ "oldnum" ~ */
      strt_idx               = n_layers * nc;
      stop_idx               = strt_idx + n1n2c;

      for (unsigned k = strt_idx; k < stop_idx; k++) { 
        O[k]                 = czero;
        if (iu3 > 2) {
          Z[k]               = czero;
        }
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
    
        int kdk              = 0;
        for (unsigned k = strt_idx; k < stop_idx; k++) {
          O[k]               = (a * roldub[kdk]) + (b * rnewub[kdk]);
          if (iu3 > 2) {
            Z[k]             = czero;
          }
          ++kdk;
        }
      }
      else {

        int trcount;
        physics_data.fetch("trcount", &trcount);
        run.palette.fetch("ffp",      &ffp);
        run.palette.fetch("kc",       &kc);

        for (unsigned k = 0; k < num - oldnum; k++) {
          for (unsigned l = 0; l < nc; l++) { roldub[l] = rnewub[l]; }
          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) <= kc ) {

                dummy        = ((double) rand() / RAND_MAX);
                ++trcount;
                dummy        = ((double) rand() / RAND_MAX);
                ++trcount;
            }
          }

          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) <= kc ) {

                next_real    = ffp * (((double) rand() / RAND_MAX ) * two - one);
                ++trcount;
                next_imag    = ffp * (((double) rand() / RAND_MAX ) * two - one);
                ++trcount;
                tuple        = std::complex<double>(next_real, next_imag);
                rnewub[l]    = tuple;

            }
          }
        }
        int kdk              = 0;
        for (unsigned k = strt_idx; k < stop_idx; k++) {
          O[k]               = (a * roldub[kdk]) + (b * rnewub[kdk]);
          if (iu3 > 2) {
            Z[k]             = czero;
          }
          ++kdk;
        }
        physics_data.reset("trcount", trcount);
      }
        run.palette.reset("oldnumub", num);
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyLineTiedBC( std::string str_step, stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank          );

  int np;            run.palette.fetch(    "np",     &np       );
  int n1n2c;         run.stack_data.fetch( "n1n2c",  &n1n2c    );
  int n_layers;      run.stack_data.fetch( "n3",     &n_layers );
  std::string model; run.palette.fetch("model",      &model    );

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx, stop_idx;

  if (     rank   == 0     ) { strt_idx = 0;                 }
  else if( rank   == np - 1) { strt_idx = ( n_layers  * nc); }

  stop_idx         = strt_idx + n1n2c;

  ComplexArray& O  = run.U0;
  ComplexArray& Z  = run.U2;

  ComplexArray& tO = run.tU0;
  ComplexArray& tZ = run.tU2;

  for (unsigned k  = strt_idx; k < stop_idx; k++) {

    if (str_step.compare("predict") == 0 || str_step.compare("finalize") == 0) {
      O[k]         = czero;
    }
    else if (str_step.compare("correct") == 0 ) {
      tO[k]        = czero;
    }
    else {
      std::cout   << "applyLineTiedBC: WARNING - unknown str_step value " << std::endl;
    }

    if (model.compare("hall") == 0 ) {
      if (str_step.compare("predict") == 0 || str_step.compare("finalize") == 0) {
        Z[k]       = czero;
      }
      else if (str_step.compare("correct") == 0 ) {
        tZ[k]      = czero;
      }
      else {
        std::cout << "applyLineTiedBC: WARNING - unknown str_step value " << std::endl;
      }
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::updatePAJ( std::string str_step, stack& run ) {

  int rank;     MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  int  np;      MPI_Comm_size(MPI_COMM_WORLD, &np);

  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c);    /* ~ number complex elements in layer            ~ */
  int n_layers; run.stack_data.fetch("iu2",   &n_layers); /* ~ number of layers in stack                   ~ */

  ComplexArray&    O = run.U0;                            /* ~ for predictor case                          ~ */
  ComplexArray&    H = run.U1;

  ComplexArray&   tO = run.tU0;                           /* ~ for corrector case                          ~ */
  ComplexArray&   tH = run.tU1;

  RealArray&      k2 = run.k2;                            /* ~ square-magnitude of k-space vectors         ~ */
  RealArray&  inv_k2 = run.inv_k2;                        /* ~ inverse square magnitude of k-space vectors ~ */

  RealVar ssqd; run.palette.fetch(  "ssqd",  &ssqd);      /* ~ parameter sigma^2 relating A to H           ~ */

  unsigned kstart    = 0;                                 /* ~ lower and upper loop limits on k            ~ */
  unsigned kstop     = n_layers * n1n2c;                  /* ~ note: pbot and atop boundary layers incl.'d ~ */

  unsigned idx       = 0;                                 /* ~ index for k2 and inv_k2                     ~ */

  if (    str_step.compare("predict") == 0 ) {            /* ~ PRIOR to predictor step                     ~ */
    for (unsigned k = kstart; k < kstop; k++) { 

       if (k % n1n2c == 0){ idx = 0; }                    /* ~ reset idx when starting new layer           ~ */

       P[k] = inv_k2[idx] * O[k];                         /* ~ O = -delperp^2 P                            ~ */
       A[k] = H[k] / ( one + ssqd*k2[idx] );
       J[k] = k2[idx] * A[k];                             /* ~ J = -delperp^2 A                            ~ */

       ++idx;

    }
  }
  else if(str_step.compare("correct") == 0 ) {            /* ~ do as above using predictor results         ~ */
    for (unsigned k = kstart; k < kstop; k++) {

      if (k % n1n2c == 0){ idx = 0;}

      if (rank != 0 && rank != np-1) {
      P[k] = inv_k2[idx] * tO[k];                         /* ~ P, A, and J are now ready for use in        ~ */
      }
      else {
        if (rank == 0 && k < n1n2c ) {
          P[k] = inv_k2[idx] * tO[k];                     /* ~ P, A, and J are now ready for use in        ~ */
        }
        else if (rank == np-1 && k < 4*n1n2c ) {
          P[k] = inv_k2[idx] * tO[k];                     /* ~ P, A, and J are now ready for use in        ~ */
        }
      }
      A[k] = tH[k] / ( one + ssqd*k2[idx] );
      J[k] = k2[idx] * A[k];

      ++idx;

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::updateTimeInc( stack& run ) {

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank );

    RealArray glbMaxU;
    glbMaxU.reserve(maxU.size());

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
