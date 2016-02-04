/* class redhallmhd (implementation)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * For defining and implementing the reduced mhd Hall physics for coronos.
 * this class is responsible for "filling" lcsolve's data structures with the 
 * appropriate values - based on its "knowledge" of the physical model of
 * the plasma. These values are needed by lcsolve so that lcsolve can update 
 * its stack.
 *
 */

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
    init_physics_data( run       );    /* ~ physics - specific parameters              ~ */
    initU(             run       );    /* ~ initialization of layers 1 - n3 of U       ~ */

  int srun;
  run.palette.fetch("srun", &srun);

  initBoundaries(    run         );   /* ~ initialization of quantities needed for     ~ */
                                      /* ~ boundary value application.                 ~ */
  if (srun == 1) {

  run.writeUData  (              );   /* ~ initial conditions report                   ~ */

  }

  initialize(        run         );   /* ~ not a good name, I'll probably revise this  ~ */
                                      /* ~ it might be to let initialize do everything ~ */
                                      /* ~ here or, alternatively to do away with it   ~ */
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

    std::cout << "initNoDrive: initializing p for rank " << rank << std::endl;
    std::cout << "initNoDrive: n3  = "                   << n3   << std::endl;
    std::cout << "initNoDrive: iu2 = "                   << iu2  << std::endl;

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

void redhallmhd::initialize (stack& run ) {

    std::string model;
    run.palette.fetch("model", &model);

#ifdef HAVE_CUDA_H
 
#else

            initTimeInc(         run );
 
            std::string init;
            run.palette.fetch("initMode", &init);
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
      else if(model.compare("rmhd") == 0 ) {
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
    else if(model.compare("rmhd") == 0) {
      A[k]         = H[k];                             /* ~ NOTE: why not use known values?            ~ */
    }

    ++idx;

  }

  for (unsigned k = 0; k < usize; k++) {U1[k] = A[k];} /* ~  U1 now holds Fourier transform of A       ~ */

//  H.resize(0);                                         /* ~ H is discarded                             ~ */

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
