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
 *        FILE: Implemenation of class "redhallmhd"
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

  init_physics_data( run         );   /* ~ physics - specific parameters               ~ */
  initU(             run         );   /* ~ initialization of layers 1 - n3 of U        ~ */
  initBoundaries(    run         );   /* ~ initialization of quantities needed for     ~ */
                                      /* ~ boundary value application.                 ~ */

  initialize(        run         );   /* ~ not a good name, I'll probably revise this  ~ */
                                      /* ~ it might be to let initialize do everything ~ */
                                      /* ~ here or, alternatively to do away with it   ~ */
  int srun;
  run.palette.fetch("srun", &srun);
  if (srun == 1) {

    fftw.fftwReverseAll( run,   J);
    run.writeUData  (            );   /* ~ initial conditions report                   ~ */

  }
  else {

    /* ~ NOTE! AUX should be read from data file here!                                 ~ */
   
  }

  initIRMHD(         run         );   /* ~ checks for inhomegeneous conditions         ~ */

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

  std::string init;

  int srun;

  fftw.fftwInitialize( run );

  MPI_Barrier(MPI_COMM_WORLD);

  run.palette.fetch("srun", &srun);

  if (srun == 1) {

    run.palette.fetch("initMode", &init);

    if (init.compare("fourierspace") == 0) computeFourierU( run );
    if (init.compare("realspace")    == 0) computeRealU(    run );
    if (init.compare("from_data" )   == 0) readUData(       run );

    int ilnr;
    run.palette.fetch("ilnr", &ilnr);

    if (ilnr != 0 && init.compare("from_data") !=0 ) pLinzEnv( run );

  }
  else { std::cout << "reading data for subrun " << srun << std::endl;                                  

                                           readUData(       run );

       }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::init_physics_data( stack& run ) {

  std::string model;
  run.palette.fetch("model",   &model  ); 

  int      bdrys;
  run.palette.fetch("bdrys",   &bdrys  );

  RealVar   tstart;
  run.palette.fetch("tstart" , &tstart );
  RealVar nu;  
  run.palette.fetch("nu"     , &nu     );
  RealVar eta; 
  run.palette.fetch("eta"    , &eta    );

  RealVar qs0    = nu;
  RealVar qs1    = eta;

  /* ~ hall - related ~ */

  RealVar delta;
  RealVar epratio;
  RealVar beta; 
  RealVar kappa;
  RealVar qs2; 
  RealVar qs3;
  RealVar ssqd;
  RealVar rho;

  if ( model.compare("hall") == 0 )  {           /* ~ initialize only if needed  ~ */

   run.palette.fetch("delta"  , &delta  );
   run.palette.fetch("epratio", &epratio);
   run.palette.fetch("beta"   , &beta   );
   run.palette.fetch("kappa"  , &kappa  );

   qs2           = kappa + (half * beta * eta);
   qs3           = nu;
   ssqd          = two * delta * sqrt(epratio);
   rho           = sqrt(beta) * delta;

  }

  std::string pname;
  std::string padjust;

  padjust.assign("adj" );
  pname.assign( "t_cur");
  physics_data.emplace(pname, tstart, padjust);
  pname.assign( "dtvb");
  physics_data.emplace(pname, zero, padjust);


  if (bdrys > 0) {

    int brcount  = 0;
    int trcount  = 0;

     padjust.assign("adj");


     pname.assign("brcount");
     physics_data.emplace(pname, brcount, padjust);
     pname.assign("trcount");
     physics_data.emplace(pname, trcount, padjust);

  }

  padjust.assign("rfx");

  pname.assign( "qs0" );
  physics_data.emplace(pname, qs0,   padjust);
  pname.assign( "qs1" );
  physics_data.emplace(pname, qs1,   padjust);

  if ( model.compare("hall") == 0 )  {

    pname.assign( "qs2" );
    physics_data.emplace(pname, qs2,   padjust);
    pname.assign( "qs3" );
    physics_data.emplace(pname, qs3,   padjust);

    pname.assign( "ssqd");
    physics_data.emplace(pname, ssqd,  padjust);
    pname.assign( "rho");
    physics_data.emplace(pname, rho,  padjust);

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::computeRealU( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int n1; 
  run.stack_data.fetch( "n1",    &n1    );
  int n2; 
  run.stack_data.fetch( "n2",    &n2    );
  int n3; 
  run.stack_data.fetch( "n3",    &n3    );
  int n1n2;
  run.stack_data.fetch( "n1n2",  &n1n2  );
  int n_flds;
  run.stack_data.fetch("iu3",    &n_flds);

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

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1; 
  run.stack_data.fetch( "n1",    &n1    );
  int n2; 
  run.stack_data.fetch( "n2",    &n2    );
  int n3; 
  run.stack_data.fetch( "n3",    &n3    );
  int n1n2; 
  run.stack_data.fetch( "n1n2",  &n1n2  );
  int n1n2c; 
  run.stack_data.fetch( "n1n2c", &n1n2c );
  int iu1; 
  run.stack_data.fetch("iu1",    &iu1   );
  int iu2; 
  run.stack_data.fetch("iu2",    &iu2   );
  int iu3;
  run.stack_data.fetch("iu3",    &iu3   );

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

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int srun;
  run.palette.fetch("srun", &srun);

  std::string data_file;
  data_file                              = run.getLastDataFilename(srun-1);
  const char *c_data_file;
  c_data_file                            = data_file.c_str();

  std::ifstream ifs;
  ifs.open( c_data_file, std::ios::in );

  if ( ifs.good() ) {

    InputOutputArray& U                  = run.U;
    
    int iu3; 
    run.stack_data.fetch("iu3",  &iu3);
    int n1; 
    run.stack_data.fetch("n1",   &n1);
    int n2; 
    run.stack_data.fetch("n2",   &n2);
    int n3; 
    run.stack_data.fetch("n3",   &n3);
    int n1n2;
    run.stack_data.fetch("n1n2", &n1n2);

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
    
    while ( !ifs.eof() ) {

      if (slab_index > n3) break; 

      ifs >> next_p;
      ++point_count;

      ifs >> next_a;
      ++point_count;

      U[to_row_maj_idx][slab_index][0]   = next_p;
      U[to_row_maj_idx][slab_index][1]   = next_a;

      if(iu3 > 2) {

        ifs >> next_bz;
        ++point_count;
        ifs >> next_vz;
        ++point_count;

        U[to_row_maj_idx][slab_index][2] = next_bz;
        U[to_row_maj_idx][slab_index][3] = next_vz;

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

  int    n1n2;
  run.stack_data.fetch("n1n2", &n1n2);
  int    n3;
  run.stack_data.fetch("n3",   &n3  );
  RealVar zl;
  run.palette.fetch(   "zl",   &zl  );

  InputOutputArray&  U  = run.U;
  RealArray&         z  = run.z;

  int i, j;

  for (i = 0; i < n1n2; ++i) {
    for (j = 0; j < n3+1; ++j) {

      U[i][j][0]        = U[i][j][0] * (1.0 - (std::abs(z[j] - (0.5*zl))/(0.5*zl)));

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initBoundaries( stack& run) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  

  int np;
  run.palette.fetch("np", &np);

  int bdrys;

  if ( rank == 0 || rank == np - 1 ) {
    
  /* ~ do exterior boundary stuff ~ */  
 
    run.palette.fetch("bdrys", &bdrys);

    if (bdrys > 0 ) { initFootPointDriving( run ); }
    else            { initNoDrive(          run ); }

  }

  /* ~ probably want to put an MPI_Barrier here so the stacks 
   *   with exterior layers don't get out of sync with the others ~ */

  /* ~ do everything else  ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::countModes( stack& run ) {

  RealVar kc;
  run.palette.fetch("kc", &kc);


  bool   l_reset_success = false;

  int m                  = 0;
  while ( ((RealVar) m) <= kc ) { ++m; } 

  RealVar arg;
  int nf                 = 0;
  for (int ix = -m; ix < (m + 1);  ++ix) {
    for (int iy = -m; iy < (m + 1);  ++iy) {

    arg                  = sqrt( pow((two_pi * ((RealVar) ix)),2) + pow((two_pi * ((RealVar) iy)),2) );
    if (arg != zero && arg <= kc) { ++nf; }
    
    }
  }
 
  l_reset_success        = run.palette.reset("nf", nf);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initFootPointDriving( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  countModes( run );

  int    nf;
  run.palette.fetch("nf"  , &nf  );
  RealVar tauC;
  run.palette.fetch("tauC", &tauC);
  RealVar tauE;
  run.palette.fetch("tauE", &tauE);

  RealVar dtau          = ((RealVar) (1.383)) * tauC;
  RealVar qfp           = zero;
  RealVar ffp           = (two_pi / tauE) * (one / sqrt( (RealVar) nf));

  bool l_reset_success  = false;

  l_reset_success       = run.palette.reset("dtau", dtau);
  l_reset_success       = run.palette.reset("qfp" , qfp);
  l_reset_success       = run.palette.reset("ffp" , ffp);

  const int    seed     = 1234567;
  srand(seed);

  int rcount;
  run.palette.fetch(   "rcount",  &rcount);

  RealVar dummy;
  if (rcount > 0) { for (int l = 0; l < rcount; ++l) { dummy = rand(); } }

  int srun;
  run.palette.fetch(   "srun",     &srun );

  int n1n2c;
  run.stack_data.fetch("n1n2c",   &n1n2c);

  RealVar kc;
  run.palette.fetch(   "kc", &kc);

  RealArray& k2         = run.k2;

  RealVar    next_real; 
  RealVar    next_imag;
  ComplexVar tuple;

  if (rank == 0) {

    int brcount;

    physics_data.fetch("brcount", &brcount);

    roldlb.reserve(n1n2c);
    rnewlb.reserve(n1n2c);

    roldub.reserve(n1n2c);
    rnewub.reserve(n1n2c);

    if (srun == 1) {

      for (int l = 0; l < n1n2c; ++l) {

        roldlb.push_back(czero);

        if ( sqrt(k2[l]) < kc ) {

            next_real   = ffp * (rand() * two - one);
            ++brcount;
            next_imag   = ffp * (rand() * two - one);
            ++brcount;
            tuple       = ComplexVar(next_real, next_imag);
            rnewlb.push_back(tuple);

/* ~        next_real = ffp * (rand() * two - one);
            ++brcount;
            next_imag = ffp * (rand() * two - one);
            ++brcount;
            tuple = std::complex<double>(next_real, next_imag);
            rnewub.push_back(tuple);
 ~ */
        }
      }
  
      l_reset_success = physics_data.reset("brcount", brcount);

    }
    else {

           std::cout << "initFootPointDriving: read rmct2r.in not yet implemented" << std::endl;

    }
  }

  int bdrys;
  run.palette.fetch("bdrys", &bdrys);

  if (bdrys  == 2) {

    int np;
    run.palette.fetch("np", &np);

    if (rank == np - 1 ) {

      int trcount;
      physics_data.fetch("trcount", &trcount);

      for (int l = 0; l < n1n2c; ++l) {

        roldub.push_back(czero);

        if ( sqrt(k2[l]) < kc ) {

          if (srun == 1) {

            dummy       = rand();
            ++trcount;
            dummy       = rand();
            ++trcount;

            next_real   = ffp * (rand() * two - one);
            ++trcount;
            next_imag   = ffp * (rand() * two - one);
            ++trcount;
            tuple       = ComplexVar(next_real, next_imag);
            rnewub.push_back(tuple);

          }
          else {

            dummy       = rand();
            ++trcount;
            dummy       = rand();
            ++trcount;

           std::cout << "initFootPointDriving: read rmct2r.in not yet implemented" << std::endl;

          }
        }
      }
      l_reset_success   = physics_data.reset("trcount", trcount);
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void redhallmhd::finalizeFootPointDriving( stack& run, lcsolve& solve ) {
//
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  
//    int bdrys;
//    run.palette.fetch("bdrys", &bdrys);
//  
//    if (bdrys > 0) {
//      if (rank == 0) {
//  
//         roldlb.resize(0);
//         roldub.resize(0);
//         rnewlb.resize(0);
//         rnewub.resize(0);
//   
//       }
//     }
//  }

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initNoDrive( stack& run) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  int srun;
  run.palette.fetch(   "srun" ,  &srun );
  int np;
  run.palette.fetch(    "np"  ,  &np   );      /* ~ number of processes ~ */
  int n1n2; 
  run.stack_data.fetch( "n1n2",  &n1n2 );
  int n3; 
  run.stack_data.fetch( "n3"  ,  &n3   );      /* ~ number of layers    ~ */
  int iu2; 
  run.stack_data.fetch("iu2",    &iu2   );     /* ~ n3 + 2              ~ */
  int iu3;
  run.stack_data.fetch("iu3"  ,  &iu3   );     /* ~ number of fields    ~ */

  InputOutputArray& U = run.U;

  if ( rank == 0 || rank == np - 1) {

//    std::cout << "initNoDrive: initializing p for rank " << rank << std::endl;
//    std::cout << "initNoDrive: n3  = "                   << n3   << std::endl;
//    std::cout << "initNoDrive: iu2 = "                   << iu2  << std::endl;

    if ( rank == 0 ) {
      for (int k   = 0; k< n1n2; ++k) { U[k][0][0]      = zero; }
      for (int k   = 0; k< n1n2; ++k) { U[k][n3+1][1]   = zero; }  /* ~ atop is zero on rank 0                   ~ */
                                                                   /* ~ NOTE: should only be done for first run! ~ */
      if ( iu3  > 2) {
        for (int k   = 0; k< n1n2; ++k) { U[k][0][2]    = zero; }
        for (int k   = 0; k< n1n2; ++k) { U[k][n3+1][3] = zero; }  /* ~ vztop is zero on rank 0                  ~ */
                                                                   /* ~ NOTE: should only be done for first run! ~ */
      }
    }

    if ( rank == np - 1 ) {
      for (int k   = 0; k< n1n2; ++k) { U[k][n3][0]     = zero; }
      for (int k   = 0; k< n1n2; ++k) { U[k][n3+1][1]   = zero; }  /* ~ atop is zero on rank 0                   ~ */
                                                                   /* ~ NOTE: should only be done for first run! ~ */
      if ( iu3  > 2) {
        for (int k   = 0; k< n1n2; ++k) { U[k][n3][2]   = zero; }
        for (int k   = 0; k< n1n2; ++k) { U[k][n3+1][3] = zero; }  /* ~ vztop is zero on rank 0                  ~ */
                                                                   /* ~ NOTE: should only be done for first run! ~ */
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
 

    int n1n2c;
    run.stack_data.fetch("n1n2c", &n1n2c   );          /* ~ number of complex elements per layer       ~ */
    int n_layers;                           
    run.stack_data.fetch("iu2",   &n_layers);          /* ~ number of layers in stack                  ~ */
  
    RealArray&    k2 = run.k2;                         /* ~ square magnitudes of k-space vectors       ~ */
    ComplexArray& U0 = run.U0;                         /* ~ holds phi (i.e. P ) at this point          ~ */ 
  
    ComplexArray::size_type usize;
    usize            = U0.capacity();                  /* ~ current capacity of U0 - should be known   ~ */
  
    assert(usize     == (n1n2c * n_layers));           /* ~ test usize                                 ~ */
  
    ComplexArray O;                                    /* ~ temporary storage for vorticity            ~ */
    O.reserve(usize);
  
    P.reserve(usize);                                  /* ~ member P will be needed throughout run     ~ */
  
   for (unsigned k = 0; k <usize; k++) {P[k] = U0[k];} /* ~ preserve stream funtion in P               ~ */
   
    unsigned  idx    = 0;                              /* ~ index for k2                               ~ */
    for (unsigned k  = 0; k < usize; k++) {
  
      if (k % n1n2c  == 0 ) { idx = 0; }               /* ~ reset idx when starting new layer          ~ */
      O[k] = k2[idx] * P[k];                           /* ~ Omega = - delperp^2 P                      ~ */
  
      ++idx;
  
    }
  
   for (unsigned k = 0; k <usize; k++) {U0[k] = O[k];} /* ~ U0 now holds Fourier transform of Vorticity ~ */
                                                       /* ~ and P is initialized                        ~ */
    O.resize(0);                                       /* ~ dispense with temporary storage             ~ */
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::HfromA( stack& run )  {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n1n2c; 
    run.stack_data.fetch("n1n2c", &n1n2c);               /* ~ number of complex elements per layer       ~ */
    int n_layers;
    run.stack_data.fetch("iu2",   &n_layers);            /* ~ number of layers in stack                  ~ */
  
    RealArray& k2    = run.k2;                           /* ~ square magnitude of k-space vectors        ~ */
    ComplexArray& U1 = run.U1;                           /* ~ holds A (flux function) at this point      ~ */
  
    ComplexArray::size_type usize;
    usize            = U1.capacity();                    /* ~ current capacity of U1 - should be known   ~ */
  
    assert(usize     == (n1n2c * n_layers));             /* ~ test usize                                 ~ */
  
    ComplexArray H;                                      /* ~ temporary storage for H-function           ~ */
    H.reserve(usize);
  
    A.reserve(usize);                                    /* ~ members A and J (current density) needed   ~ */
    J.reserve(usize);                                    /* ~ throughout run                             ~ */
  
    RealVar ssqd;                                        /* ~ parameter sigma^2 needed for H             ~ */
    physics_data.fetch( "ssqd", &ssqd);
  
    std::string model;
    run.palette.fetch("model", &model);
  
    for (unsigned k = 0; k < usize; k++) {A[k] = U1[k];} /* ~ preserve flux function in A                ~ */

    unsigned idx     = 0;                                /* ~ index for k2                               ~ */
    for (unsigned k  = 0; k < usize; k++) {
  
      if (k % n1n2c  == 0 ) { idx = 0; }                 /* ~ reset idx when starting new layer          ~ */
  
      J[k]           = k2[idx] * A[k];                   /* ~ J = -delperp A                             ~ */
      if (model.compare("hall") == 0 ) {
        H[k]         = A[k] + ssqd * J[k];               /* ~ H = A + sigma^2 J                          ~ */
      }
      else if (model.compare("rmhd") == 0 || model.compare("inhm") == 0) {
        H[k]         = A[k];                             /* ~ H = A for reduced-MHD                      ~ */
      }
      ++idx;
  
    }
  
    for (unsigned k = 0; k < usize; k++) {U1[k] = H[k];} /* ~ U1 now holds Fourier transform of H        ~ */
                                                         /* ~ and A and J are both initialized           ~ */
    H.resize(0);                                         /* ~ dispense with temporary storage            ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::PfromO( stack& run )  {
 
  int n1n2c;
  run.stack_data.fetch("n1n2c", &n1n2c      );        /* ~ number of complex elements per layer        ~ */
  int n_layers;
  run.stack_data.fetch("iu2",   &n_layers);           /* ~ number of layers in stack                   ~ */

  RealArray&    inv_k2 = run.inv_k2;                  /* ~ inverse-square magnitude of k-space vectors ~ */
  ComplexArray& U0     = run.U0;                      /* ~ holds Omega (vorticity) at this point       ~ */

  ComplexArray::size_type usize;
  usize                = U0.capacity();               /* ~ current capacity of U0 - should be known    ~ */

  assert(usize         == (n1n2c * n_layers));        /* ~ test usize                                  ~ */

  ComplexArray O;                                     /* ~ temporary storage for vorticity             ~ */
  O.reserve(usize);
                                                      /* ~ note: P is already known                    ~ */

  for (unsigned k = 0; k <usize; k++) {O[k] = U0[k];} /* ~ not necessary. Could use U0 directly        ~ */

  unsigned idx         = 0;                           /* ~ index for inv_k2                            ~ */

  for (unsigned k = 0; k < usize; k++) {

    if ( k % n1n2c == 0) { idx = 0; }                 /* ~ reset idx when starting new layer           ~ */
    P[k] = inv_k2[idx] * O[k];                        /* ~ Omega = - delperp^2 P                       ~ */
                                                      /* ~ NOTE: why not use known value?              ~ */
    ++idx;

  }

  for (unsigned k = 0; k <usize; k++) {U0[k] = P[k];} /* ~ U0 now holds Fourier transform of P         ~ */

//  O.resize(0);                                      /* ~ vorticity is discarded here                 ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::AfromH( stack& run )  {

  int n1n2c;
  run.stack_data.fetch("n1n2c", &n1n2c);               /* ~ number of complex elements per layer       ~ */
  int n_layers;
  run.stack_data.fetch("iu2",   &n_layers);            /* ~ number of layers in stack                  ~ */

  RealArray& k2    = run.k2;                           /* ~ square magnitude of k-space vectors        ~ */
  ComplexArray& U1 = run.U1;                           /* ~ holds H function at this point             ~ */

  ComplexArray::size_type usize;
  usize            = U1.capacity();                    /* ~ current capacity of U1 - should be known   ~ */

  assert(usize     == (n1n2c * n_layers));             /* ~ test usize                                 ~ */

  ComplexArray H;                                      /* ~ temporary storage for H-function           ~ */
  H.reserve(usize);
                                                       /* ~ note: current values of A and J are known  ~ */

  RealVar ssqd;
  physics_data.fetch("ssqd", &ssqd);                   /* ~ parameter sigma^2 needed for A             ~ */

  std::string model;
  run.palette.fetch("model", &model);

  for (unsigned k = 0; k < usize; k++) {H[k] = U1[k];}

  unsigned idx     = 0;                                /* ~ index for k2                               ~ */
  for (unsigned k  = 0; k < usize; k++) {

    if (k % n1n2c  == 0) { idx = 0; }                  /* ~ reset idx when starting new layer          ~ */

    if (model.compare("hall") == 0 ) {
      A[k]         = H[k] / (one + ssqd*k2[idx]);      /* ~ NOTE: why not use known values?            ~ */
    }
    else if(model.compare("rmhd") == 0 || model.compare("inhm") == 0 ) {
      A[k]         = H[k];                             /* ~ NOTE: why not use known values?            ~ */
    }

    ++idx;

  }

  for (unsigned k = 0; k < usize; k++) {U1[k] = A[k];} /* ~  U1 now holds Fourier transform of A       ~ */

//  H.resize(0);                                         /* ~ H is discarded                             ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalElls(    stack& run ) {

  std::string nprofile;
  run.palette.fetch("nprofile", &nprofile);
  int p3;
  run.palette.fetch("p3", &p3);

  Elln.assign(p3,zero);
  EllA.assign(p3,zero);
  EllB.assign(p3,zero);

   h11.assign(p3,zero);
   h12.assign(p3,zero);
   h21.assign(p3,zero);
   h22.assign(p3,zero);

  int i_profile;

  if      (nprofile.compare("flat")   == 0) { i_profile =  0; }
  else if (nprofile.compare("torus")  == 0) { i_profile =  1; }
  else if (nprofile.compare("cloop")  == 0) { i_profile =  2; }
  else                                      { i_profile = -1; }

  switch(i_profile) {

  case(0) : Elln.assign(p3, zero); 
            EllA.assign(p3, zero);
            EllB.assign(p3, zero);
            break;
  case(1) : for ( int k = 0; k <  p3;  ++k) {

              EllA[k] = abs( dvalfdz[k] / (two  * valfven[k]) );
              Elln[k] = abs( dndz[k]    / (four * nofz[k]   ) );
//            Elln[k] = EllA[k];

            }
            EllB.assign(p3, zero);
            break;
  case(2) : Elln.assign(p3, zero);
            EllA.assign(p3, zero);
            EllB.assign(p3, zero);
            break;

  default :

    std::cout << "evalN: WARNING - the profile " << nprofile << " is not implemented. Assuming a flat profile." << std::endl;

    Elln.assign(p3, zero);
    EllA.assign(p3, zero);
    EllB.assign(p3, zero);

  }

  for ( int k = 0; k <  p3;  ++k) {

    h11[k] =  umean[k]   * ( EllA[k] + EllB[k] - Elln[k] );
    h12[k] = -valfven[k] * ( EllA[k] + EllB[k] + Elln[k] );
    h21[k] =  h12[k];
    h22[k] =  umean[k]   * ( EllB[k] - Elln[k] - EllA[k] );

  }

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalValf( stack& run ) {
  
  std::string nprofile;
  run.palette.fetch("nprofile", &nprofile);                          /* ~ specification for density profile    ~ */
  int p3;
  run.palette.fetch("p3", &p3);                                      /* ~ layers per stack not incl. ghosts      ~ */
  double zl;                   
  run.palette.fetch("zl", &zl);                                      /* ~ size of domain along z               ~ */
  double ninf;
  run.palette.fetch("ninf", &ninf);                                  /* ~ density at z = +/- infinity          ~ */
  double n0;
  run.palette.fetch("n0", &n0);                                      /* ~ density at z_0 (i.e. z = zl / 2)     ~ */

  double H0;                                                         /* ~ density scale-height constant        ~ */

  RealArray& z = run.z;

//  RealArray nofz;
//  RealArray dndz;

  int i_profile;
  if      (nprofile.compare("flat")   == 0) { i_profile =  0; }      /* ~ switch statements don't like strings ~ */
  else if (nprofile.compare("torus")  == 0) { i_profile =  1; }
  else if (nprofile.compare("cloop")  == 0) { i_profile =  2; }
  else                                      { i_profile = -1; }

  switch(i_profile) {

  case(0) : valfven.assign(p3, one );
            dvalfdz.assign(p3, zero);
               nofz.assign(p3, one );
               dndz.assign(p3, zero);
            break;
  case(1) : H0 = zl / sqrt( - eight * log( (one - ninf) / (n0 - ninf) ));

            std::cout << "evalValf: H0 = " << H0 << std::endl;

               nofz.assign(p3,one );
               dndz.assign(p3,zero);
            valfven.assign(p3,one );
            dvalfdz.assign(p3,zero);

            for (unsigned k = 0; k < p3; ++k) {

               nofz[k]    = ninf + ((n0 - ninf) * (exp(-half*pow((z[k] - (half*zl))/H0 ,two))));
               dndz[k]    = -(n0 - ninf) * exp( -half * pow(((z[k]-(half*zl))/H0),2)) * (z[k] - (half*zl)) / ( pow(H0,2) );
               valfven[k] = sqrt( n0 / nofz[k] );
               dvalfdz[k] = -half * valfven[k] * dndz[k] / nofz[k];

//               valfven[k] = sqrt(n0 / (ninf + (n0-ninf) * (exp(-half*pow((z[k] - (half*zl))/H0 ,two)))));

            }
            break;
  case(2) : valfven.assign(p3, one ); 
            dvalfdz.assign(p3, zero);
               nofz.assign(p3, one );
               dndz.assign(p3, zero);
            break;

  default : 

    std::cout << "evalValf: WARNING - the profile " << nprofile << " is not implemented. Assuming a flat profile." << std::endl;
    valfven.assign(p3, one );
    dvalfdz.assign(p3, zero);

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::evalUmean( stack& run ) {
  
  std::string uprofile;
  run.palette.fetch("uprofile", &uprofile);
  int p3;
  run.palette.fetch("p3", &p3);

  int i_profile;

  if      (uprofile.compare("noflow")   == 0) { i_profile =  0; }
  else if (uprofile.compare("uniform")  == 0) { i_profile =  1; }
  else if (uprofile.compare("whoknows") == 0) { i_profile =  2; }
  else                                        { i_profile = -1; }

  switch(i_profile) {

  case(0) : umean.assign(p3, zero); break;
  case(1) : umean.assign(p3, zero); break;
  case(2) : umean.assign(p3, zero); break;

  default : 

    std::cout << "evalValf: WARNING - the profile " << uprofile << " is not implemented. Assuming a flat profile." << std::endl;
    umean.assign(p3, one);

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::trackEnergies(double t_cur, stack& run ) {

  RealArray& EnergyQs = run.EnergyQs;

  static const int i_tcr =  0;
  static const int i_pe  =  1;
  static const int i_ae  =  2;
  static const int i_mo  =  3;
  static const int i_imo =  4;
  static const int i_jmo =  5;
  static const int i_oe  =  6;
  static const int i_mj  =  7;
  static const int i_imj =  8;
  static const int i_jmj =  9;
  static const int i_ce  = 10;
  static const int i_noe = 11;
  static const int i_ece = 12;
  static const int i_fe  = 13;
  static const int i_ftp = 14;
  static const int i_eds = 15;
  static const int i_dng = 16;
  static const int i_irc = 17;
  static const int i_gml = 18;
  static const int i_cns = 19;
  static const int i_ttc = 20;
  static const int i_dt  = 21;
  static const int i_dtv = 22;
  static const int i_coc = 23;
  static const int i_vkt = 24;
  static const int i_avm = 25;
  static const int i_avp = 26;
  static const int i_fp  = 27;

  double tcr;
  double pe ;
  double ae ;
  double mo ;
  double imo;
  double jmo;
  double oe ;
  double mj ;
  double imj;
  double jmj;
  double ce ;
  double cee;
  double noe;
  double ece;
  double fe ;
  double ftp;
  double eds;
  double dng;
  double irc;
  double gml;
  double cns;
  double ttc;
  double dt ;
  double dtv;
  double coc;
  double vkt;
  double avm;
  double avp;
  double fp ;

  pe  = evalTotalKineticEnergy(  run );
  ae  = evalTotalMagneticEnergy( run );
  oe  = evalTotalVorticitySqd(   run );
  ce  = evalTotalCurrentSqd(     run );
  cee = evalTotalGradCurrentSqd( run );
  fp  = evalTotalFootPointKE(    run );
  fe  = evalTotalPoyntingFlux(   run );

  double t_old;
  double dt_old;
  double aeold;
  double peold;
  double dtvb;
  double tauC;
  double eta;
  double nu;
  run.palette.fetch("tauC", &tauC);
  run.palette.fetch("eta", &eta);
  run.palette.fetch("nu", &nu);

  physics_data.fetch("dtvb", &dtvb);

  if (EnergyQs.size() >= i_fp ) {

   
    t_old  = EnergyQs[i_tcr];
    dt_old = EnergyQs[i_dt];
 
    aeold  = EnergyQs[i_ae];
    peold  = EnergyQs[i_pe];
    tcr    = t_cur;
// pe
// ae
    mo     = zero;
    imo    = zero;
    jmo    = zero;
// oe
    mj     = zero;
    imj    = zero;
    jmj    = zero;
//ce
    noe    = nu  * oe;
    ece    = eta * ce;
//fe
    ftp    = (EnergyQs[i_ftp] + fe * dt_old) / t_cur;
    eds    = (EnergyQs[i_eds] + (noe + ece) * dt_old) / t_cur;
    dng    = (EnergyQs[i_dng] + (ae - aeold + pe - peold)) /t_cur;
    irc    = (ae - aeold + pe - peold) / dt_old;
    gml    = ftp - eds;
    cns    = abs(  ( (fe - noe - ece) - ((ae - aeold + pe - peold) / dt_old)) * dt_old );
    ttc    = t_cur / tauC; 
    run.palette.fetch("dt",&dt);
    dtv    = dtvb;
    if (cee == zero){
      coc  = zero; 
    }
    else {
      coc  = sqrt(ce/cee);
    }
    if (ae != zero) {

    vkt    = EnergyQs[i_vkt] * t_old;
    vkt    = (vkt +  (ce/ae)  * dt_old) / t_cur;

    }
    else { vkt = zero; }
    avm    = pow(EnergyQs[i_avm],2) * t_old;
    avm    = sqrt((avm + ae * dt_old) / t_cur);
    avp    = half * pow(EnergyQs[i_avp],2) * t_old;
    avp    = sqrt(two * (avp + fp * dt_old) / t_cur);

    EnergyQs[ i_tcr ] = tcr;
    EnergyQs[ i_pe  ] = pe ;
    EnergyQs[ i_ae  ] = ae ;
    EnergyQs[ i_mo  ] = mo ;
    EnergyQs[ i_imo ] = imo;
    EnergyQs[ i_jmo ] = jmo;
    EnergyQs[ i_oe  ] = oe ;
    EnergyQs[ i_mj  ] = mj ;
    EnergyQs[ i_imj ] = imj;
    EnergyQs[ i_jmj ] = jmj;
    EnergyQs[ i_ce  ] = ce ;
    EnergyQs[ i_noe ] = noe;      
    EnergyQs[ i_ece ] = ece;
    EnergyQs[ i_fe  ] = fe ;
    EnergyQs[ i_ftp ] = ftp;
    EnergyQs[ i_eds ] = eds;
    EnergyQs[ i_dng ] = dng;
    EnergyQs[ i_irc ] = irc;
    EnergyQs[ i_gml ] = gml;
    EnergyQs[ i_cns ] = cns;
    EnergyQs[ i_ttc ] = ttc;
    EnergyQs[ i_dt  ] = dt ;
    EnergyQs[ i_dtv ] = dtv;
    EnergyQs[ i_coc ] = coc;
    EnergyQs[ i_vkt ] = vkt;
    EnergyQs[ i_avm ] = avm;
    EnergyQs[ i_avp ] = avp;
    EnergyQs[ i_fp  ] =  fp;
  }
  else {

    run.palette.fetch("dt", &dt_old);
    t_old  = t_cur;

    aeold  = zero;
    peold  = zero;

    tcr    = t_cur;
// pe
// ae
    mo     = zero;
    imo    = zero;
    jmo    = zero;
// oe
    mj     = zero;
    imj    = zero;
    jmj    = zero;
//ce
    noe    = nu  * oe;
    ece    = eta * ce;

/* ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ */
    noe    = zero;
    ece    = zero;
/* ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ */

//fe
    ftp    = zero;
    eds    = zero;
    cns    = abs(  ( (fe - noe - ece) - ((ae - aeold + pe - peold) / dt_old)) * dt_old );

    dng    = zero;
    gml    = (ae - aeold + pe - peold )/ dt;
    dt     = dt_old;
    dtv    = dtvb;
    if (cee == zero){
      coc  = zero; 
    }
    else {
      coc  = sqrt(ce/cee);
    }
    vkt    = zero;
    avm    = zero;
    avp    = zero;

    EnergyQs.push_back(tcr);
//    EnergyQs.push_back(pe );
//    EnergyQs.push_back(ae );
/* ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ */
    EnergyQs.push_back(peold );
    EnergyQs.push_back(aeold );
/* ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ TEST  ~ */

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

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::reportEnergyQs ( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank  );

  if (rank == 0) {

    RealArray& EnergyQs = run.EnergyQs;

    std::string prefix;
    run.palette.fetch("prefix",     &prefix );

    std::string run_label;
    run.palette.fetch("run_label",  &run_label);

    std::string res_str;
    run.stack_data.fetch("res_str", &res_str);

    std::string energy_data_file = prefix + '_' + res_str + ".o" + run_label;

    const char *c_data_file;
    c_data_file              = energy_data_file.c_str();

    std::ofstream ofs;
    ofs.open( c_data_file, std::ios::out | std::ios::app );

    if (ofs.good()) {

      unsigned esize  = EnergyQs.size();
      for (unsigned k = 0; k < esize; k++) {

        ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << EnergyQs[k] << " ";

      }
      ofs << std::endl;
      ofs.close();

    }
    else {std::cout << "reportEnergyQs: Warning - could not open file " << energy_data_file << std::endl;}

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalKineticEnergy ( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int wsize;
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  int np;
  run.palette.fetch("np", &np);

  assert (wsize == np);
  
  RealArray& k2       = run.k2;

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;
  run.stack_data.fetch("iu2"   , &iu2  );

  int n3;
  run.stack_data.fetch("n3"   , &n3    );
  double dz;
  run.stack_data.fetch("dz"   , &dz    );

  double pe;
  double pe_sum;

  pe        = zero;

  int idx   = 0;

  int kstart;
  int kstop;

  if (rank  == 0) { kstart = 0;     }
  else            { kstart = n1n2c; }

  kstop  = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++){

    if (kdx % n1n2c  == 0) { idx = 0; }

      pe    = pe + k2[idx] * pow(abs(P[kdx]), 2);

      if (( rank == np - 1 ) && ( kdx >= (four * n1n2c) )) {
        pe = pe - half * k2[idx] * pow(abs(P[kdx]),2);
      }
      if ((rank  == 0      ) && ( kdx <   n1n2c)           ) {
        pe = pe - half * k2[idx] * pow(abs(P[kdx]),2);
      }

      ++idx;
  }

  pe = pe * dz;

  int i_red;

  i_red = MPI_Reduce(&pe, &pe_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {

    pe_sum = two_thirds * pe_sum;  /* ~ NOTE!: why the two_thirds. Seems necessary but at moment I don't know why ~ */

  }

  return pe_sum;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalVorticitySqd( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int wsize;
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  int np;
  run.palette.fetch("np", &np);
  
  RealArray& k2       = run.k2;

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;
  run.stack_data.fetch("iu2"   , &iu2  );

  int n3;
  run.stack_data.fetch("n3"   , &n3    );
  double dz;
  run.stack_data.fetch("dz"   , &dz    );

  double oe;
  double oe_sum;

  oe        = zero;

  int idx   = 0;

  int kstart;
  int kstop;

  if (rank  == 0) { kstart = 0;     }
  else            { kstart = n1n2c; }

// kstart =  n1n2c; 
  kstop  = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++){

    if (kdx % n1n2c  == 0) { idx = 0; }

      oe    = oe + k2[idx]*k2[idx] * pow(abs(P[kdx]), 2);

      if (( rank == np - 1 ) && ( kdx >= (four * n1n2c) )) {
        oe = oe - half * k2[idx]*k2[idx] * pow(abs(P[kdx]),2);
      }
      if ((rank  == 0      ) && ( kdx <   n1n2c)           ) {
        oe = oe - half * k2[idx]*k2[idx] * pow(abs(P[kdx]),2);
      }

      ++idx;
  }

  oe = oe * dz;

  int i_red;

  i_red = MPI_Reduce(&oe, &oe_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {

    oe_sum = two * two_thirds * oe_sum;  /* ~ NOTE!: why the two_thirds. Seems necessary but at moment I don't know why ~ */

//    if (EnergyQs.size() >= (i_oe + 1)) {
//      EnergyQs[i_oe] = oe_sum;
//    }
//    else {
//      EnergyQs.push_back(oe_sum);
//    }

  }
  return oe_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalMagneticEnergy ( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int wsize;
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  int np;
  run.palette.fetch("np", &np);

  assert (wsize == np);
  
  RealArray& k2       = run.k2;

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;
  run.stack_data.fetch("iu2"   , &iu2  );

  int n3;
  run.stack_data.fetch("n3"   , &n3    );
  double dz;
  run.stack_data.fetch("dz"   , &dz    );

  double me;
  double me_sum;

  me        = zero;

  int idx   = 0;

  int kstart;
  int kstop;

  kstart = n1n2c; 
  kstop  = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    if (kdx % n1n2c  == 0) { idx = 0; }

      me    = me + k2[idx] * pow(abs(A[kdx]), 2);
      ++idx;
  }

  me = me * dz;

  int i_red;

  i_red = MPI_Reduce(&me, &me_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//  if (rank == 0) {
//
//    if (EnergyQs.size() >= (i_me + 1)) {
//      EnergyQs[i_me] = me_sum;
//    }
//    else {
//      EnergyQs.push_back(me_sum);
//    }
//  }
    return me_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalCurrentSqd( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int wsize;
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  int np;
  run.palette.fetch("np", &np);

  assert (wsize == np);
  
  RealArray& k2       = run.k2;

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;
  run.stack_data.fetch("iu2"   , &iu2  );

  int n3;
  run.stack_data.fetch("n3"   , &n3    );
  double dz;
  run.stack_data.fetch("dz"   , &dz    );

  double ce;
  double ce_sum;

  ce        = zero;

  int idx   = 0;

  int kstart;
  int kstop;

  kstart = n1n2c; 
  kstop  = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    if (kdx % n1n2c  == 0) { idx = 0; }

      ce    = ce + k2[idx] * k2[idx] * pow(abs(A[kdx]), 2);
      ++idx;
  }

  ce = two * ce * dz;

  int i_red;

  i_red = MPI_Reduce(&ce, &ce_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//  if (rank == 0) {
//
//    if (EnergyQs.size() >= (i_ce + 1)) {
//      EnergyQs[i_ce] = ce_sum;
//    }
//    else {
//      EnergyQs.push_back(ce_sum);
//    }
//  }
    return ce_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalGradCurrentSqd( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int wsize;
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  int np;
  run.palette.fetch("np", &np);

  assert (wsize == np);
  
  RealArray& k2       = run.k2;

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;
  run.stack_data.fetch("iu2"   , &iu2  );

  int n3;
  run.stack_data.fetch("n3"   , &n3    );
  double dz;
  run.stack_data.fetch("dz"   , &dz    );

  double cee;
  double cee_sum;

  cee        = zero;

  int idx   = 0;

  int kstart;
  int kstop;

  kstart = n1n2c; 
  kstop  = n1n2c * (iu2 - 1);

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    if (kdx % n1n2c  == 0) { idx = 0; }

      cee    = cee + k2[idx] * k2[idx] * k2[idx] * pow(abs(A[kdx]), 2);
      ++idx;
  }

  cee = two * cee * dz;

  int i_red;

  i_red = MPI_Reduce(&cee, &cee_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//  if (rank == 0) {
//
//    if (EnergyQs.size() >= (i_ce + 1)) {
//      EnergyQs[i_ce] = ce_sum;
//    }
//    else {
//      EnergyQs.push_back(ce_sum);
//    }
//  }
    return cee_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalFootPointKE( stack& run ) {

  /* ~ Note: this is the equivalent to vbot from the old code. I'm not convinced this is correct in
   *         either of the codes. But I'm getting in implemented for now. ~ */

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int wsize;
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  int np;
  run.palette.fetch("np", &np);


  assert (wsize == np);
  RealArray& k2       = run.k2;

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;
  run.stack_data.fetch("iu2"   , &iu2  );

  int n3;
  run.stack_data.fetch("n3"   , &n3    );
  double dz;
  run.stack_data.fetch("dz"   , &dz    );

  double fp;
  double fp_sum;

  fp           = zero;

  if (rank     == 0) {

    int idx    = 0;

    int kstart = 0;
    int kstop  = n1n2c;

    for (unsigned kdx = kstart; kdx < kstop; kdx++) {

      fp       = fp + k2[kdx] * pow(abs(P[kdx]), 2);
    }

    fp         = fp * dz;

//    int i_red;
//   i_red = MPI_Reduce(&fp, &fp_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    fp         = two_thirds * fp;  /* ~ NOTE!:  do I need this? ~ */

  }
    return fp;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

double redhallmhd::evalTotalPoyntingFlux ( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int wsize;
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  int np;
  run.palette.fetch("np", &np);

  assert (wsize == np);
  
  RealArray& k2       = run.k2;

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c );
  int iu2;
  run.stack_data.fetch("iu2"   , &iu2  );

  int n3;
  run.stack_data.fetch("n3"   , &n3    );
  double dz;
  run.stack_data.fetch("dz"   , &dz    );

  double fe;
  double fe_sum;

  fe        = zero;

  int idx   = 0;

  int kstart;
  int kstop;

  if       (rank  ==  0)      { kstart = 0;                 }
  else if  (rank  ==  np - 1) { kstart = ( iu2 - 2 )*n1n2c; }

  kstop  = kstart + n1n2c;

  for (unsigned kdx = kstart; kdx < kstop; kdx++){

    if (kdx % n1n2c  == 0) { idx = 0; }

      if (( rank == np - 1 ) && ( kdx >= ( n3 * n1n2c) )) {
        fe = fe + k2[idx] * ( A[kdx].real() * P[kdx].real() + A[kdx].imag() * P[kdx].imag() );
      }
      if ((rank  == 0      ) && ( kdx <   n1n2c)           ) {
        fe = fe - k2[idx] * ( A[kdx + n1n2c].real() * P[kdx].real() + A[kdx + n1n2c].imag() * P[kdx].imag() );
      }

      ++idx;
  }
  fe = two * fe;

  int i_red;
  i_red = MPI_Reduce(&fe, &fe_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {

    fe_sum = two_thirds * fe_sum;  /* ~ NOTE!: why the two_thirds. Seems necessary but at moment I don't know why ~ */

  }
  return fe_sum;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyBC( std::string str_step, stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bdrys;
  run.palette.fetch("bdrys", &bdrys);
  int np;
  run.palette.fetch("np"   , &np   );

  if ( rank == 0 || rank == np - 1) {

    if (bdrys > 0 ) { applyFootPointDrivingBC(   run ); }
    else            { applyLineTiedBC( str_step, run ); }

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyFootPointDrivingBC( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank     );

  int bdrys;
  run.palette.fetch("bdrys"   , &bdrys    );

  int np;  
  run.palette.fetch("np"      , &np       );
  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c    );
  int n_layers;
  run.stack_data.fetch("n3"   , &n_layers );

  int oldnum, num;

  RealVar    ffp, kc;
  RealVar    dummy;
  RealVar    next_real, next_imag;
  ComplexVar tuple;

  RealVar dtau;
  run.palette.fetch("dtau"    , &dtau     );
  RealVar t_cur; 
  physics_data.fetch("t_cur"  , &t_cur    );

  RealVar lowtau; 
  RealVar bigtau; 

  RealVar a, b;                                   /* ~ i.e. Gilson's "a" and "b" ~ */

  ComplexArray::size_type nc = n1n2c;

  int iu3;
  run.stack_data.fetch("iu3", &iu3);

  unsigned strt_idx, stop_idx;

  if ( rank  == 0 || rank == (np - 1) ) {    /* ~ pevol "starts" here ~ */

    RealArray&    k2       = run.k2;
    ComplexArray& O        = run.U0;
    ComplexArray& Z        = run.U2; 

    if ( rank  == 0 ) {

      run.palette.fetch("oldnumlb", &oldnum ); /* ~ "numold" ~ */
      strt_idx             = 0;
      stop_idx             = strt_idx + n1n2c;

      ComplexArray& oldr   = roldlb;
      ComplexArray& newr   = rnewlb;
      for (unsigned k = strt_idx; k < stop_idx; k++) { 
        O[k]               = czero;
        if (iu3 > 2){ 
        Z[k]               = czero;
        }
      }
      num                  = oldnum;
     
      while (  (num * dtau) <= t_cur ) { ++num; }

      bigtau               = (num * dtau) - t_cur;
      lowtau               = bigtau       - dtau;
      a                    = cos((pi*lowtau)/( two * dtau)); /* ~ Gilson's "interp" ~ */
      b                    = cos((pi*bigtau)/( two * dtau));
      if (num == oldnum) {
    
        for (unsigned k = strt_idx; k < stop_idx; k++) {

          O[k]             = (a * oldr[k]) + (b * newr[k]);
          if (iu3 > 2){ 
            Z[k]           = czero;
          }

        }
      }
      else {

        int brcount;
        physics_data.fetch("brcount", &brcount);
        physics_data.fetch("ffp",     &ffp);
        run.palette.fetch( "kc",      &kc);

        for (unsigned k = 0; k < num - oldnum; k++) {

          for (unsigned l = 0; l < nc; l++) { oldr[l] = newr[l]; }

          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) < kc ) {

                dummy      = rand();
                ++brcount;
                dummy      = rand();
                ++brcount;
            }
          }

          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) < kc ) {

                next_real  = ffp * (rand() * two - one);
                ++brcount;
                next_imag  = ffp * (rand() * two - one);
                ++brcount;
                tuple      = ComplexVar(next_real, next_imag);
                newr[l]    = tuple;

            }
          }
        }
        physics_data.reset("brcount", brcount);

      }
    }
    else if ( rank == np - 1 && bdrys == 2) {

      run.palette.fetch("oldnumub", &oldnum ); /* ~ "oldnum" ~ */
      strt_idx             = (n_layers - 1) * nc;
      stop_idx             = strt_idx + n1n2c;

      ComplexArray& oldr   = roldub;
      ComplexArray& newr   = rnewub;
      for (unsigned k = strt_idx; k < stop_idx; k++) { 
        O[k]               = czero;
        if (iu3 > 2) {
          Z[k]             = czero;
        }
      }
      num                  = oldnum;
      while (  (num * dtau) <= t_cur ) { ++num; }
      bigtau               = (num * dtau) - t_cur;
      lowtau               = bigtau       - dtau;
      a                    = cos((pi*lowtau)/( two * dtau));    /* ~ Gilson's "interp" ~ */
      b                    = cos((pi*bigtau)/( two * dtau));
      if (num == oldnum) {
    
        for (unsigned k = strt_idx; k < stop_idx; k++) {
          O[k]             = (a * oldr[k]) + (b * newr[k]);
          if (iu3 > 2) {
            Z[k]           = czero;
          }
        }
      }
      else {

        int trcount;
        physics_data.fetch("trcount", &trcount);
        physics_data.fetch("ffp",     &ffp);
        run.palette.fetch("kc",       &kc);

        for (unsigned k = 0; k < num - oldnum; k++) {
//        for (unsigned l = 0; l < nc; l++) { oldr[l] = newr[l]; }
//        for (unsigned l = 0; l < nc; l++ ) {
//          if ( sqrt(k2[l]) < kc ) {

//              dummy      = rand();
//              ++trcount;
//              dummy      = rand();
//              ++trcount;
//          }
//        }

//        for (unsigned l = 0; l < nc; l++ ) {
//          if ( sqrt(k2[l]) < kc ) {

//              next_real  = ffp * (rand() * two - one);
//              ++trcount;
//              next_imag  = ffp * (rand() * two - one);
//              ++trcount;
//              tuple      = std::complex<double>(next_real, next_imag);
//              newr[l]    = tuple;

//          }
//        }
        }
//      physics_data.reset("trcount", trcount);
      }
    }
    
//    if (num == oldnum) {
//
//      for (unsigned k = strt_idx; k < stop_idx; k++) {
//    
//        O[k]             = (a * oldr[k]) + (b * newr[k]);
//        Z[k]               = czero;
//    
//      }
//
//    }
//    else {
//
//      for (unsigned k = 0; k < num - numold; k++) { 
//
//        for (unsigned l = 0; l < nc; l++) { oldr[k] = newr[k]; }
//
//      }
//
//      /* ~ reset newr ~ */
//      /* ~ increment either trcount or brcount and reset ~ */
//
//    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyLineTiedBC( std::string str_step, stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank    );  

  int np;  
  run.palette.fetch(   "np"   , &np      );
  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n_layers;
  run.stack_data.fetch("n3"   , &n_layers);

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx, stop_idx;

  if (     rank   == 0     ) { strt_idx = 0;                 }
  else if( rank   == np - 1) { strt_idx = ( n_layers  * nc); }

  stop_idx         = strt_idx + n1n2c;

  std::string model;
  run.palette.fetch("model", &model);

  ComplexArray& O  = run.U0;
  ComplexArray& Z  = run.U2;

  ComplexArray& tO = run.tU0;
  ComplexArray& tZ = run.tU2;

  for (unsigned k  = strt_idx; k < stop_idx; k++) {

      if (str_step.compare("predict") == 0 ) {
        O[k]           = czero;
      }
      else if (str_step.compare("correct") == 0 ) {
        tO[k]          = czero;
      }
      else {
        std::cout << "applyLineTiedBC: WARNING - unknown str_step value " << std::endl;
      }

      if (model.compare("hall") == 0 ) {
        if (str_step.compare("predict") == 0 ) {
          Z[k]         = czero;
        }
        else if (str_step.compare("correct") == 0 ) {
          tZ[k]        = czero;
        }
        else {
          std::cout << "applyLineTiedBC: WARNING - unknown str_step value " << std::endl;
        }
      }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::updatePAJ( std::string str_step, stack& run ) {

  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c);         /* ~ number complex elements in layer            ~ */

  int n_layers;
  run.stack_data.fetch("iu2",   &n_layers);      /* ~ number of layers in stack                   ~ */

  ComplexArray&    O = run.U0;                   /* ~ for predictor case                          ~ */
  ComplexArray&    H = run.U1;

  ComplexArray&   tO = run.tU0;                  /* ~ for corrector case                          ~ */
  ComplexArray&   tH = run.tU1;

  RealArray&      k2 = run.k2;                   /* ~ square-magnitude of k-space vectors         ~ */
  RealArray&  inv_k2 = run.inv_k2;               /* ~ inverse square magnitude of k-space vectors ~ */

  RealVar ssqd;                                  /* ~ parameter sigma^2 relating A to H           ~ */
  physics_data.fetch(  "ssqd",  &ssqd);

  unsigned kstart    = 0;                        /* ~ lower and upper loop limits on k            ~ */
  unsigned kstop     = n_layers * n1n2c;         /* ~ note: pbot and atop boundary layers incl.'d ~ */

  unsigned idx;                                  /* ~ index for k2 and inv_k2                     ~ */
  if (    str_step.compare("predict") == 0 ) {   /* ~ PRIOR to predictor step                     ~ */
    for (unsigned k = kstart; k < kstop; k++) { 

       if (k % n1n2c == 0){ idx = 0; }           /* ~ reset idx when starting new layer           ~ */

       P[k] = inv_k2[idx] * O[k];                /* ~ O = -delperp^2 P                            ~ */
       A[k] = H[k] / ( one + ssqd*k2[idx] );
       J[k] = k2[idx] * A[k];                    /* ~ J = -delperp^2 A                            ~ */

       ++idx;

    }
  }
  else if(str_step.compare("correct") == 0 ) {   /* ~ do as above using predictor results         ~ */
    for (unsigned k = kstart; k < kstop; k++) {

      if (k % n1n2c == 0){ idx = 0;}

      P[k] = inv_k2[idx] * tO[k];                /* ~ P, A, and J are now ready for use in        ~ */
      A[k] = tH[k] / ( one + ssqd*k2[idx] );     /* ~ corrector step                              ~ */
      J[k] = k2[idx] * A[k];

      ++idx;

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::updateTimeInc( stack& run ) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    RealArray glbMaxU;
    glbMaxU.reserve(maxU.size());

    int n1n2;
    run.stack_data.fetch("n1n2", &n1n2);
    RealVar dt;
    run.palette.fetch("dt", &dt);
    RealVar dz;
    run.stack_data.fetch("dz", &dz);
 
    RealVar q1, q2, qp;

    run.palette.fetch("q1", &q1);
    run.palette.fetch("q2", &q2);
    run.palette.fetch("qp", &qp);

    RealVar dtvb;
    RealVar dtr;

    int i_red;
    i_red       = MPI_Reduce(&maxU[0], &glbMaxU[0], maxU.size(), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
    int i_brd;
    i_brd       = MPI_Bcast(&glbMaxU[0], maxU.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

/* ~ currently a stub, but eventually intended to replace the stuff at end of Loop in lcsolve ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~ Destructor ~ */

redhallmhd:: ~redhallmhd() {

//  fftw.fftwFinalize();

}
