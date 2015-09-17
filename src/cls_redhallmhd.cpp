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

  init_physics_data( run       );     /* ~ physics - specific parameters               ~ */
  initU(             run       );     /* ~ initialization of layers 1 - n3 of U        ~ */

  int srun;
  run.palette.fetch("srun", &srun);

  if (srun == 1) {

  run.writeUData(              );     /* ~ initial conditions report                   ~ */

  }

  initBoundaries(    run       );     /* ~ initialization of quantities needed for     ~ */
                                      /* ~ boundary value application.                 ~ */

  initialize(        run       );     /* ~ not a good name, I'll probably revise this  ~ */
                                      /* ~ it might be to let initialize do everything ~ */
                                      /* ~ here or, alternatively to do away with it   ~ */
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

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

  double   tstart;
  run.palette.fetch("tstart" , &tstart );
  double nu;  
  run.palette.fetch("nu"     , &nu     );
  double eta; 
  run.palette.fetch("eta"    , &eta    );

  double qs0  = nu;
  double qs1  = eta;

  /* ~ hall - related ~ */

  double delta;
  double epratio;
  double beta; 
  double kappa;
  double qs2; 
  double qs3;
  double ssqd;
  double rho;

  if ( model.compare("hall") == 0 )  {           /* ~ initialize only if needed  ~ */

   run.palette.fetch("delta"  , &delta  );
   run.palette.fetch("epratio", &epratio);
   run.palette.fetch("beta"   , &beta   );
   run.palette.fetch("kappa"  , &kappa  );

   qs2        = kappa + (half * beta * eta);
   qs3        = nu;
   ssqd       = two * delta * sqrt(epratio);
   rho        = sqrt(beta) * delta;

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
//  run.palette.emplace( pname, qs0,   padjust);
  pname.assign( "qs1" );
  physics_data.emplace(pname, qs1,   padjust);
//  run.palette.emplace( pname, qs1,   padjust);

  if ( model.compare("hall") == 0 )  {

    pname.assign( "qs2" );
    physics_data.emplace(pname, qs2,   padjust);
//    run.palette.emplace( pname, qs2,   padjust);
    pname.assign( "qs3" );
    physics_data.emplace(pname, qs3,   padjust);
//    run.palette.emplace( pname, qs3,   padjust);

    pname.assign( "ssqd");
    physics_data.emplace(pname, ssqd,  padjust);
//    run.palette.emplace(pname, ssqd,  padjust);
    pname.assign( "rho");
    physics_data.emplace(pname, rho,  padjust);
//    run.palette.emplace(pname, rho,  padjust);

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

  int n_lyrs          = n3;

  InputOutputArray& U = run.U;
  RealArray&        x = run.x;
  RealArray&        y = run.y;

  int idx             = 0;

  for (int i_f = 0; i_f < n_flds; ++i_f) {

    switch(i_f) {

    case(0) :
      for (int i_x=0;i_x < n1; ++i_x) {
        for (int j_y=0;j_y < n1; ++j_y) {

          idx             = (i_x * n1) + j_y;

          U[idx][n3][i_f] = -half * 0.001L * ( cos(two_pi*x[i_x]) - cos(two_pi * y[j_y]) );

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

//  fft fftw;

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

  double               real_part;
  double               imag_part;

  std::complex<double> tuple;

  unsigned idx        = 0;

  for (int i_f = 0; i_f < n_flds; ++i_f) {

     for (unsigned k  = 0; k < n1n2c; ++k) { Cin[k]    = czero; }
     for (unsigned k  = 0; k < n1n2;  ++k) { Rout[k]   =  zero; }
     
     switch(i_f) {

     case(0) :
       idx            =        0 * (n2/2 + 1) + 1;
       real_part      =  2.5e-04 * n1n2;
       imag_part      =  0.0;
       tuple          = std::complex<double>(real_part, imag_part);
       Cin[idx]       = tuple;

       idx            =        1 * (n2/2 + 1);
       real_part      = -2.5e-04 * n1n2;
       imag_part      =  0.0;
       tuple          = std::complex<double>(real_part, imag_part);
       Cin[idx]       = tuple;

       idx            = (n1-1) * (n2/2 + 1);
       Cin[idx]       = tuple;

       break;
     case(1) :
       idx            =       1 * (n2/2 + 1) + 1;
       real_part      =    -0.1 * n1n2;
       imag_part      =     0.0;
       tuple          = std::complex<double>(real_part, imag_part);
       Cin[idx]       = tuple;

       idx            =  (n1-1) * (n2/2 + 1) + 1;
       real_part      =     0.1 * n1n2;
       imag_part      =     0.0;
       tuple          = std::complex<double>(real_part, imag_part);
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

  for (int i_f = 0; i_f < n_flds; ++i_f) {
    for (int i_l = 1; i_l < n_lyrs + 1; ++i_l) {
    
    for (int k = 0; k< n1n2; ++k) { U[k][i_l][i_f] = U[k][n3][i_f]; }

    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~ STOPPED HERE ~ */

void redhallmhd::readUData( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int srun;
  run.palette.fetch("srun", &srun);

  std::cout << "readUData << srun = " << srun << std::endl;

  std::string      data_file;
//  if (srun == 1) { srun = 0; }
  data_file = run.getLastDataFilename(srun-1);
  std::cout << "readUData: opening file - " << data_file << " for reading..." << std::endl;
  const char    *c_data_file;
  c_data_file = data_file.c_str();

  std::ifstream    ifs;
  ifs.open( c_data_file, std::ios::in );

  if ( ifs.good() ) {

    InputOutputArray& U  = run.U;
    
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

    int n_slab_points    = n1n2 * iu3;
    int point_count      = 0;
    int slab_index       = 1;
    int from_col_maj_idx = 0;
    int to_row_maj_idx   = 0;
    int i                = 0;
    int j                = 0;

    double next_p; 
    double next_a; 

    double next_bz; 
    double next_vz;
    
    while ( !ifs.eof() ) {

      if (slab_index > n3) break; 

      ifs >> next_p;
      ++point_count;

      ifs >> next_a;
      ++point_count;

      U[to_row_maj_idx][slab_index][0] = next_p;
      U[to_row_maj_idx][slab_index][1] = next_a;

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
      else to_row_maj_idx  = 0;

      if (from_col_maj_idx == n1n2) {

        from_col_maj_idx   = 0;
        i                  = 0;
        j                  = 0;
      }

      if(point_count == n_slab_points) {

        point_count        = 0;
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
  double zl;
  run.palette.fetch(   "zl",   &zl  );

  InputOutputArray&  U  = run.U;
  RealArray&         z  = run.z;

  int i, j;

  for (i = 0; i < n1n2; ++i) {
    for (j = 0; j < n3; ++j) {

      U[i][j][0] = U[i][j][0] * (1.0 - (std::abs(z[j] - (0.5*zl))/(0.5*zl)));

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

  double kc;
  run.palette.fetch("kc", &kc);


  bool   l_reset_success = false;

  int m   = 0;
  while ( ((double) m) <= kc ) { ++m; } 

  double arg;
  int nf  = 0;
  for (int ix = -m; ix < (m + 1);  ++ix) {
    for (int iy = -m; iy < (m + 1);  ++iy) {

    arg   = sqrt( pow((two_pi * ((double) ix)),2) + pow((two_pi * ((double) iy)),2) );
    if (arg != zero && arg <= kc) { ++nf; }
    
    }
  }
 
  l_reset_success = run.palette.reset("nf", nf);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initFootPointDriving( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  countModes( run );

  int    nf;
  run.palette.fetch("nf"  , &nf  );
  double tauC;
  run.palette.fetch("tauC", &tauC);
  double tauE;
  run.palette.fetch("tauE", &tauE);

  double dtau          = ((double) (1.383)) * tauC;
  double qfp           = zero;
  double ffp           = (two_pi / tauE) * (one / sqrt( (double) nf));

  bool l_reset_success = false;

  l_reset_success      = run.palette.reset("dtau", dtau);
  l_reset_success      = run.palette.reset("qfp" , qfp);
  l_reset_success      = run.palette.reset("ffp" , ffp);

  const int    seed    = 1234567;
  srand(seed);

  int rcount;
  run.palette.fetch(   "rcount",  &rcount);

  double dummy;
  if (rcount > 0) { for (int l = 0; l < rcount; ++l) { dummy = rand(); } }

  int srun;
  run.palette.fetch(   "srun",     &srun );

  int n1n2c;
  run.stack_data.fetch("n1n2c",   &n1n2c);

  double kc;
  run.palette.fetch(   "kc", &kc);

  RealArray& k2        = run.k2;

  double               next_real; 
  double               next_imag;
  std::complex<double> tuple;

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
//      roldub.push_back(czero);

        if ( sqrt(k2[l]) < kc ) {

            next_real = ffp * (rand() * two - one);
            ++brcount;
            next_imag = ffp * (rand() * two - one);
            ++brcount;
            tuple = std::complex<double>(next_real, next_imag);
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

            dummy     = rand();
            ++trcount;
            dummy     = rand();
            ++trcount;

            next_real = ffp * (rand() * two - one);
            ++trcount;
            next_imag = ffp * (rand() * two - one);
            ++trcount;
            tuple     = std::complex<double>(next_real, next_imag);
            rnewub.push_back(tuple);

          }
          else {

            dummy     = rand();
            ++trcount;
            dummy     = rand();
            ++trcount;

           std::cout << "initFootPointDriving: read rmct2r.in not yet implemented" << std::endl;

          }
        }
      }
      l_reset_success = physics_data.reset("trcount", trcount);
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void redhallmhd::finalizeFootPointDriving( stack& run, lcsolve& solve ) {

//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  int bdrys;
//  run.palette.fetch("bdrys", &bdrys);
//
//  if (bdrys > 0) {
//    if (rank == 0) {
//
//       roldlb.resize(0);
//       roldub.resize(0);
//       rnewlb.resize(0);
//       rnewub.resize(0);
// 
//     }
//   }
//}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::initNoDrive( stack& run) {

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void redhallmhd::writeUData( stack& run ) {
//
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD,&rank  );
//
//  int    srun;
//  run.palette.fetch("srun",    &srun  );
//
//  std::string data_file;
//  data_file   = run.getLastDataFilename(srun-1);
//  const char  *c_data_file;
//  c_data_file = data_file.c_str();
//
//  std::ofstream ofs;
//  ofs.open( c_data_file, std::ios::out );
//
//  if ( ofs.good() ) {
//
//    InputOutputArray& U = run.U;
//    
//    int iu3; 
//    run.stack_data.fetch("iu3",  &iu3);
//    int n1; 
//    run.stack_data.fetch("n1",   &n1);
//    int n2; 
//    run.stack_data.fetch("n2",   &n2);
//    int n3; 
//    run.stack_data.fetch("n3",   &n3);
//    int n1n2;
//    run.stack_data.fetch("n1n2", &n1n2);
//
//    int n_slab_points    = n1n2 * iu3;
//    int point_count      = 0;
//    int slab_index       = 1;
//    int from_col_maj_idx = 0;
//    int to_row_maj_idx   = 0;
//    int i                = 0;
//    int j                = 0;
//
//    double next_p; 
//    double next_a; 
//    double next_bz; 
//    double next_vz;
//    
//    while ( slab_index < n3 + 1 ) {
//
////    U[to_row_maj_idx][slab_index][0] = next_p;
////    U[to_row_maj_idx][slab_index][1] = next_a;
//
//      ++point_count;
//      next_p = U[to_row_maj_idx][slab_index][0];
//      ++point_count;
//      next_a = U[to_row_maj_idx][slab_index][1];
//
//      if(iu3 > 2) {
//
////        U[to_row_maj_idx][slab_index][2] = next_bz;
////        U[to_row_maj_idx][slab_index][3] = next_vz;
//
//        ++point_count;
//        next_bz = U[to_row_maj_idx][slab_index][2];
//        ++point_count;
//        next_vz = U[to_row_maj_idx][slab_index][3];
//
//      }
//
////      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_p << " ";
////      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_a << " ";
////
//      ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_p << " ";
//      ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_a << " ";
//
//      if (iu3 > 2)  {
//
////      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_bz << " ";
////      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_vz << " ";
//
//        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_bz << " ";
//        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_vz << " ";
//
//      }
//
//      ofs << std::endl;
//     
//      if (from_col_maj_idx < n1n2) {
//
//        ++from_col_maj_idx;
//        if (from_col_maj_idx % n2 != 0) ++j;
//        else {
//          j = 0;
//          ++i;
//        }
//      }
//
//      if (to_row_maj_idx < n1n2 - 1) to_row_maj_idx = i + (j*n1);
//      else to_row_maj_idx  = 0;
//
//      if (from_col_maj_idx == n1n2) {
//
//        from_col_maj_idx   = 0;
//        i                  = 0;
//        j                  = 0;
//      }
//
//      if(point_count == n_slab_points) {
//
//        point_count        = 0;
//        ++slab_index;
//      }
//    }
//
//    ofs.close();
//
//  }
//}

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

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

//  int rank;
//  int np,  n1n2c, n_layers;
//  unsigned strt_idx, stop_idx;
//
//  MPI_Comm_rank(MPI_COMM_WORLD,   &rank  );
//  run.stack_data.fetch("n1n2c", &n1n2c   );
//  run.stack_data.fetch("n3"   , &n_layers);
//  run.palette.fetch(   "np"   , &np      );
//
//  ComplexArray& P = run.U0;
////  ComplexArray& A = run.U1;

//  ComplexArray::size_type nc = n1n2c;
//
//  if (rank == 1) {
//
//    strt_idx  = nc;
//    stop_idx  = strt_idx + nc;
//
//       for (unsigned k = strt_idx; k < stop_idx; k++) {
//
//          if (abs(A[k].real()) >= 1.0e-12) {
//          std::cout << std::setw(12) << std::right << std::setprecision(4) << std::scientific << A[k].real() << std::endl;
//          }
//          else {
//          std::cout << std::setw(12) << std::right << std::setprecision(4) << std::scientific << zero        << std::endl;
//          }
//
//       }
//
//    }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */


/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

//   int rank;
//   int np,  n1n2c, n_layers;
//   unsigned strt_idx, stop_idx;
// 
//   MPI_Comm_rank(MPI_COMM_WORLD,   &rank);
// 
//   run.stack_data.fetch("n1n2c", &n1n2c   );
//   run.stack_data.fetch("n3"   , &n_layers);
//   run.palette.fetch(   "np"   , &np      );
// 
//   ComplexArray& O            = solve.U1;
//   ComplexArray::size_type nc = n1n2c;
// 
//   if (rank == 1) {
// 
//     strt_idx = nc;
// 
//   }
//   else if (rank == 0) {
// 
//     strt_idx = (n_layers + 1) * nc;
// 
//   }
//   stop_idx   = strt_idx + n1n2c;
// 
//   if (rank == 1) {
// 
//     for (unsigned k = strt_idx; k < stop_idx; k++) {
// 
//       O[k] = cone;
   
   //      std::cout << std::setw(24) << std::right << std::setprecision(16) << std::scientific << O[k].real() << std::endl;
   
//       }
//     }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

#endif

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::OfromP( stack& run )  {
 

    int n1n2c;
    run.stack_data.fetch("n1n2c", &n1n2c   );      /* ~ number of complex elements per layer       ~ */
    int n_layers;                           
    run.stack_data.fetch("iu2",   &n_layers);      /* ~ number of layers in stack                  ~ */
  
    RealArray&    k2 = run.k2;                     /* ~ square magnitudes of k-space vectors       ~ */
    ComplexArray& U0 = run.U0;                   /* ~ holds phi (i.e. P ) at this point          ~ */ 
  
    ComplexArray::size_type usize;
    usize            = U0.capacity();              /* ~ current capacity of U0 - should be known   ~ */
  
    assert(usize     == (n1n2c * n_layers));       /* ~ test usize                                 ~ */
  
    ComplexArray O;                                /* ~ temporary storage for vorticity            ~ */
    O.reserve(usize);
  
    P.reserve(usize);                              /* ~ member P will be needed throughout run     ~ */
  
//    P                = U0;                         /* ~ preserve stream funtion in P               ~ */

   for (unsigned k = 0; k <usize; k++) {P[k] = U0[k];}
   
    unsigned  idx    = 0;                          /* ~ index for k2                               ~ */
    for (unsigned k  = 0; k < usize; k++) {
  
      if (k % n1n2c  == 0 ) { idx = 0; }           /* ~ reset idx when starting new layer          ~ */
      O[k] = k2[idx] * P[k];                       /* ~ Omega = - delperp^2 P                      ~ */
  
      ++idx;
  
    }
  
//    U0               = O;                          /* ~ U0 now holds Fourier transform of Vorticity ~ */
   for (unsigned k = 0; k <usize; k++) {U0[k] = O[k];}
                                                   /* ~ and P is initialized                        ~ */
    O.resize(0);                                   /* ~ dispense with temporary storage             ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::HfromA( stack& run )  {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n1n2c; 
    run.stack_data.fetch("n1n2c", &n1n2c);         /* ~ number of complex elements per layer       ~ */
    int n_layers;
    run.stack_data.fetch("iu2",   &n_layers);      /* ~ number of layers in stack                  ~ */
  
    RealArray& k2    = run.k2;                     /* ~ square magnitude of k-space vectors        ~ */
    ComplexArray& U1 = run.U1;                   /* ~ holds A (flux function) at this point      ~ */
  
    ComplexArray::size_type usize;
    usize            = U1.capacity();              /* ~ current capacity of U1 - should be known   ~ */
  
    assert(usize     == (n1n2c * n_layers));       /* ~ test usize                                 ~ */
  
    ComplexArray H;                                /* ~ temporary storage for H-function           ~ */
    H.reserve(usize);
  
    A.reserve(usize);                              /* ~ members A and J (current density) needed   ~ */
    J.reserve(usize);                              /* ~ throughout run                             ~ */
  
    double ssqd;                                   /* ~ parameter sigma^2 needed for H             ~ */
    physics_data.fetch( "ssqd", &ssqd);
  
//    std::cout << "HfromA: ssqd = " << ssqd << std::endl;
  
    std::string model;
    run.palette.fetch("model", &model);
  
//    A                = U1;                         /* ~ preserve flux function in A                ~ */

    for (unsigned k = 0; k < usize; k++) {A[k] = U1[k];}

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

// if (rank == 0) {
//   for( unsigned k = n1n2c; k < 2*n1n2c; ++k ) {
//    if (abs(A[k].real()) >= 1.0e-10) {
//    std::cout << std::setw(12) << std::right << std::setprecision(4) << std::scientific << A[k].real() << std::endl;
//    }
//    else {
//    std::cout << std::setw(12) << std::right << std::setprecision(4) << std::scientific << zero  << std::endl;
//    }
//   }
// }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */



    unsigned idx     = 0;                          /* ~ index for k2                               ~ */
    for (unsigned k  = 0; k < usize; k++) {
  
      if (k % n1n2c  == 0 ) { idx = 0; }           /* ~ reset idx when starting new layer          ~ */
  
      J[k]           = k2[idx] * A[k];             /* ~ J = -delperp A                             ~ */
      if (model.compare("hall") == 0 ) {
        H[k]         = A[k] + ssqd * J[k];         /* ~ H = A + sigma^2 J                          ~ */
      }
      else if(model.compare("rmhd") == 0 ) {
        H[k]         = A[k];                       /* ~ H = A for reduced-MHD                      ~ */
      }
      ++idx;
  
    }
  
    for (unsigned k = 0; k < usize; k++) {U1[k] = H[k];}

//    U1               = H;                          /* ~ U1 now holds Fourier transform of H        ~ */
                                                     /* ~ and A and J are both initialized           ~ */
    H.resize(0);                                     /* ~ dispense with temporary storage            ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::PfromO( stack& run )  {
 
  int n1n2c;
  run.stack_data.fetch("n1n2c", &n1n2c      );   /* ~ number of complex elements per layer        ~ */
  int n_layers;
  run.stack_data.fetch("iu2",   &n_layers);      /* ~ number of layers in stack                   ~ */

  RealArray&    inv_k2 = run.inv_k2;             /* ~ inverse-square magnitude of k-space vectors ~ */
  ComplexArray& U0     = run.U0;                 /* ~ holds Omega (vorticity) at this point       ~ */

  ComplexArray::size_type usize;
  usize                = U0.capacity();          /* ~ current capacity of U0 - should be known    ~ */

  assert(usize         == (n1n2c * n_layers));   /* ~ test usize                                  ~ */

  ComplexArray O;                                /* ~ temporary storage for vorticity             ~ */
  O.reserve(usize);
                                                 /* ~ note: P is already known                    ~ */

//  O                    = U0;                     /* ~ not necessary. Could use U0 directly        ~ */
  for (unsigned k = 0; k <usize; k++) {O[k] = U0[k];}

  unsigned idx         = 0;                      /* ~ index for inv_k2                            ~ */
  for (unsigned k = 0; k < usize; k++) { 

    if ( k % n1n2c == 0) { idx = 0; }            /* ~ reset idx when starting new layer           ~ */
    P[k] = inv_k2[idx] * O[k];                   /* ~ Omega = - delperp^2 P                       ~ */
                                                 /* ~ NOTE: why not use known value?              ~ */
    ++idx;

  }

//  U0                   = P;                     /* ~ U0 now holds Fourier transform of P          ~ */

  for (unsigned k = 0; k <usize; k++) {U0[k] = P[k];}

  O.resize(0);                                  /* ~ vorticity is discarded here                  ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::AfromH( stack& run )  {

  int n1n2c;
  run.stack_data.fetch("n1n2c", &n1n2c);         /* ~ number of complex elements per layer       ~ */
  int n_layers;
  run.stack_data.fetch("iu2",   &n_layers);      /* ~ number of layers in stack                  ~ */

  RealArray& k2    = run.k2;                     /* ~ square magnitude of k-space vectors        ~ */
  ComplexArray& U1 = run.U1;                     /* ~ holds H function at this point             ~ */

  ComplexArray::size_type usize;
  usize            = U1.capacity();              /* ~ current capacity of U1 - should be known   ~ */

  assert(usize     == (n1n2c * n_layers));       /* ~ test usize                                 ~ */

  ComplexArray H;                                /* ~ temporary storage for H-function           ~ */
  H.reserve(usize);
                                                 /* ~ note: current values of A and J are known  ~ */

  double ssqd;
  physics_data.fetch("ssqd", &ssqd);             /* ~ parameter sigma^2 needed for A             ~ */

  std::string model;
  run.palette.fetch("model", &model);

//  H                = U1;
  for (unsigned k = 0; k < usize; k++) {H[k] = U1[k];}

  unsigned idx     = 0;                          /* ~ index for k2                               ~ */
  for (unsigned k  = 0; k < usize; k++) {

    if (k % n1n2c  == 0) { idx = 0; }            /* ~ reset idx when starting new layer          ~ */

    if (model.compare("hall") == 0 ) {
      A[k] = H[k] / (one + ssqd*k2[idx]);        /* ~ NOTE: why not use known values?            ~ */
    }
    else if(model.compare("rmhd") == 0) {
      A[k] = H[k];                               /* ~ NOTE: why not use known values?            ~ */
    }

    ++idx;

  }

//  U1               = A;                          /* ~  U1 now holds Fourier transform of A       ~ */
  for (unsigned k = 0; k < usize; k++) {U1[k] = A[k];}

  H.resize(0);                                   /* ~ H is discarded                             ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void redhallmhd::applyBC(    stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bdrys;
  run.palette.fetch("bdrys", &bdrys);
  int np;
  run.palette.fetch("np"   , &np   );

  if ( rank == 0 || rank == np - 1) {

    if (bdrys > 0 ) { applyFootPointDrivingBC( run ); }
    else            { applyLineTiedBC(         run ); }

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

  double ffp, kc;
  double dummy;
  double next_real, next_imag;
  std::complex<double> tuple;

  double dtau;
  run.palette.fetch("dtau"    , &dtau     );
  double t_cur; 
  physics_data.fetch("t_cur"  , &t_cur    );

  double lowtau; 
  double bigtau; 

  double a, b;                                   /* ~ i.e. Gilson's "a" and "b" ~ */

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
        Z[k]             = czero;
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
                tuple      = std::complex<double>(next_real, next_imag);
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
          for (unsigned l = 0; l < nc; l++) { oldr[l] = newr[l]; }
          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) < kc ) {

                dummy      = rand();
                ++trcount;
                dummy      = rand();
                ++trcount;
            }
          }

          for (unsigned l = 0; l < nc; l++ ) {
            if ( sqrt(k2[l]) < kc ) {

                next_real  = ffp * (rand() * two - one);
                ++trcount;
                next_imag  = ffp * (rand() * two - one);
                ++trcount;
                tuple      = std::complex<double>(next_real, next_imag);
                newr[l]    = tuple;

            }
          }
        }
        physics_data.reset("trcount", trcount);
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

void redhallmhd::applyLineTiedBC( stack& run ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank    );  

  std::cout << "applyLineTiedBC: WARNING! - I have not been tested! " << std::endl;

  int np;  
  run.palette.fetch(   "np"   , &np      );
  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n_layers;
  run.stack_data.fetch("n3"   , &n_layers);

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx, stop_idx;

  if (     rank   == 0     ) { strt_idx = 0;                   }
  else if( rank   == np - 1) { strt_idx = (n_layers - 1) * nc; }

  stop_idx          = strt_idx + n1n2c;

  std::string model;
  run.palette.fetch("model", &model);

  ComplexArray& O   = run.U0;
  ComplexArray& Z   = run.U2;

  for (unsigned k   = strt_idx; k < stop_idx; k++) {

    O[k]            = czero;
    if (model.compare("hall") == 0 ) {
      Z[k]          = czero;
    }

  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void redhallmhd::setS( std::string str_step, stack& run, lcsolve& solve) {

//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  std::string model;
//  run.palette.fetch("model", &model);
//
//  double  s0,  s1,  s2,   s3;                    /* ~ implicit fraction s - parameters                 ~ */
//  run.palette.fetch(  "s0",   &s0 );
//  run.palette.fetch(  "s1",   &s1 );
//
//  if (model.compare("hall") == 0) {
//
//    run.palette.fetch("s2",   &s2 );
//    run.palette.fetch("s3",   &s3 );
//
//  }
//
//  double qs0, qs1, qs2,  qs3;                    /* ~ dissipation 'q' - parameters see documentation   ~ */
//
//  physics_data.fetch("qs0", &qs0);
//  physics_data.fetch("qs1", &qs1);
//
//  if (model.compare("hall") == 0) {
//
//    physics_data.fetch("qs2", &qs2);
//    physics_data.fetch("qs3", &qs3);
//
//  }
//
//  double ge0, ge1, ge2,  ge3;
//  double gi0, gi1, gi2,  gi3;
//
//  double  dt;                                    /* ~ the current time increment                       ~ */
//  run.palette.fetch("dt",   &dt );
//
//  double  pfrac;                                 /* ~ fraction of dt to use in predictor step          ~ */
//  run.palette.fetch("pfrac", &pfrac);
//
//  if ( str_step.compare("predict") == 0 ) {      /* ~ use partial time-step for predictor case         ~ */
//    dt           = dt * pfrac;                   /* ~ dt is local to setS so no harm done here         ~ */
//  }
//
//  ge0            = (one - s0) * qs0 * dt;        /* ~ th g^(ex)'s  see documentation                   ~ */
//  ge1            = (one - s1) * qs1 * dt;
//
//  if (model.compare("hall") == 0) {
//
//    ge2          = (one - s2) * qs2 * dt;
//    ge3          = (one - s3) * qs3 * dt;
//
//  }
//
//  gi0            =        s0  * qs0 * dt;        /* ~ th g^(im)'s  see documentation                   ~ */
//  gi1            =        s1  * qs1 * dt;
//
//  if (model.compare("hall") == 0) {
//
//    gi2          =        s2  * qs2 * dt;
//    gi3          =        s3  * qs3 * dt;
//
//  }
//
//  RealArray& SE0 = solve.SE0;                    /* ~ the S^(ex)'s - see documentation                 ~ */
//  RealArray& SE1 = solve.SE1;
//
//  RealArray& SE2 = solve.SE2;
//  RealArray& SE3 = solve.SE3;
//
//  RealArray& SI0 = solve.SI0;                    /* ~ the S^(im)'s - see documentation                 ~ */
//  RealArray& SI1 = solve.SI1;
//
//  RealArray& SI2 = solve.SI2;
//  RealArray& SI3 = solve.SI3;
//
//  RealArray& k2  = run.k2;                       /* ~ square magnitude of k-space vectors              ~ */
//
//  int n1n2c;                                     /* ~ number of complex elements in layer              ~ */
//  run.stack_data.fetch( "n1n2c", &n1n2c );
//
//  RealArray::size_type nc = SE0.capacity();      /* ~ S's should already by sized. This is a check     ~ */
//  assert (nc     == n1n2c);
//
//  for (unsigned k = 0; k < n1n2c; k++) {         /* ~ there are only as many S-elements as k2 elements ~ */
//
//    SE0[k]       = ( one - (ge0 * k2[k]));       /* ~ S's are initialized. See documentation           ~ */
//    SE1[k]       = ( one - (ge1 * k2[k]));  
//
//    if (model.compare("hall") == 0) {
//
//      SE2[k]     = ( one - (ge2 * k2[k]));  
//      SE3[k]     = ( one - (ge3 * k2[k]));  
//
//    }
//
//    SI0[k]       = one / ( one + (gi0 * k2[k]));
//    SI1[k]       = one / ( one + (gi1 * k2[k]));
//
//    if (model.compare("hall") == 0) {
//
//      SI2[k]     = one / ( one + (gi2 * k2[k]));
//      SI3[k]     = one / ( one + (gi3 * k2[k]));
//
//    }
//  }
//
///* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
//
////  if (str_step.compare("correct") == 0 && rank == 1) {
//
////  unsigned kstart = n1n2c;
////  unsigned kstop  = (2 * kstart) - 1;
//
////   std::cout << "setS: one   = " << one    << "\n\n";
//
////   std::cout << "setS: dt    = " << dt     << "\n\n";
////   std::cout << "setS: pfrac = " << pfrac  << "\n\n";
//
////   std::cout << "setS: s0    = " << s0     << "\n\n";
////   std::cout << "setS: s1    = " << s1     << "\n\n";
//
////   std::cout << "setS: qs0   = " << qs0    << "\n\n";
////   std::cout << "setS: qs1   = " << qs1    << "\n\n";
//
////   std::cout << "setS: gi0   = " << gi0    << "\n\n";
////   std::cout << "setS: gi1   = " << gi1    << "\n\n";
//
////   unsigned kstart = 0;
////   unsigned kstop  = kstart + n1n2c;
////    
////   for (unsigned k = kstart; k < kstop; k++) {
//
////     std::cout << std::setw(24) << std::right << std::setprecision(16) << std::scientific << SI1[k] << std::endl;
//
// //   }
////  }
//
///* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
//
////  SE0.assign(nc, one);                         /* ~ retained for testing                             ~ */
////  SE1.assign(nc, one);
////  SE2.assign(nc, one);
////  SE3.assign(nc, one);
////
////  SI0.assign(nc, one);
////  SI1.assign(nc, one);
////  SI2.assign(nc, one);
////  SI3.assign(nc, one);

//}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void redhallmhd::setB( std::string str_step, stack& run, lcsolve& solve) {

//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  int n1n2;
//  run.stack_data.fetch( "n1n2" ,&n1n2  );
//  int n1n2c; 
//  run.stack_data.fetch( "n1n2c",&n1n2c );
//  int iu2;
//  run.stack_data.fetch( "iu2"  ,&iu2   );
//
//  unsigned kstop = n1n2c * iu2;
//
//  double rho;
//  physics_data.fetch(  "rho"  , &rho   );                     /* ~ gyro-radius parameter                                      ~ */
//  double beta, rtbeta;                           
//  run.palette.fetch(   "beta" , &beta  );                     /* ~ zero'th-order plasma-beta                                  ~ */
//  rtbeta = sqrt(beta);
//
//  ComplexArray BrKt;
//  BrKt.reserve(n1n2c * iu2);
//
//  RealArray d1x, d1y;
//  RealArray d2x, d2y;
//  RealArray d3x, d3y;
//
//  d1x.reserve(n1n2 * iu2);
//  d1y.reserve(n1n2 * iu2);
//
//  d2x.reserve(n1n2 * iu2);
//  d2y.reserve(n1n2 * iu2);
//
//  d3x.reserve(n1n2 * iu2);
//  d3y.reserve(n1n2 * iu2);
//
//  ComplexArray& B0 = solve.B0;
//  ComplexArray& B1 = solve.B1;
//  ComplexArray& B2 = solve.B2;
//  ComplexArray& B3 = solve.B3;
//
//  ComplexArray& O  = solve.U0;
//  ComplexArray& H  = solve.U1;
//  ComplexArray& Z  = solve.U2;
//  ComplexArray& V  = solve.U3;
//
//  ComplexArray& tO = solve.tU0;
//  ComplexArray& tH = solve.tU1;
//  ComplexArray& tZ = solve.tU2;
//  ComplexArray& tV = solve.tU3;
//
//  std::string model;
//  run.palette.fetch("model", &model);
//
//  solve.partialsInXandY(  run, P,      d1x, d1y);                /* ~ d1x, d1y hold real-space partials in x and y of P         ~ */
//
//  if (str_step.compare("predict") == 0 && rank == 0) {
//    for( unsigned k = n1n2; k < 2*n1n2; ++k ) {
//
//      std::cout << std::setw(16) << std::right << std::setprecision(8) << std::scientific << d1x[k] << std::endl;
//
//    }
//  }
//
//  if (str_step.compare("predict"     ) == 0) {
//    solve.partialsInXandY(run, O,      d2x, d2y);                /* ~ d2x, d2y hold real-space partials in x and y of O         ~ */
//    maxU[0] = solve.maxdU(             d2x, d2y);                /* ~                                                           ~ */
//  }
//  else if (str_step.compare("correct") == 0) {
//    solve.partialsInXandY(run, tO,     d2x, d2y);                /* ~ d2x, d2y hold real-space partials in x and y of tO        ~ */
//  }
//  solve.bracket(run, BrKt, d1x, d1y, d2x, d2y);                  /* ~ calculate [phi, Omega]                                    ~ */
//  for (unsigned k = 0; k < kstop; k++) { B0[k] = - BrKt[k]; }    /* ~ place result in B0                                        ~ */
//
//  if (model.compare("hall") == 0 ) { 
//
//     if (str_step.compare("predict"     ) == 0) {
//       solve.partialsInXandY(run, Z,      d2x, d2y);             /* ~ d2x, d2y hold real-space partials in x and y of Z         ~ */
//       maxU[2] = solve.maxdU(             d2x, d2y);             /* ~                                                           ~ */
//     }
//     else if (str_step.compare("correct") == 0) {
//       solve.partialsInXandY(run, tZ,     d2x, d2y);             /* ~ d2x, d2y hold real-space partials in x and y of tZ        ~ */
//     }
//     solve.bracket(run, BrKt, d1x, d1y, d2x, d2y);               /* ~ calculate [phi, Z]                                        ~ */
//     for (unsigned k = 0; k < kstop; k++) { B2[k] = - BrKt[k]; } /* ~ place result in B2                                        ~ */
//
//  }
//
//  solve.averageAcrossLayers( run, -1, d1x, d1y );                /* ~ calculate averages of phi_x & phi_y across adj't layers   ~ */
//  if (str_step.compare("predict"     ) == 0) {
//    solve.partialsInXandY(run, H,      d3x, d3y);                /* ~ d3x, d3y hold real-space partials in x and y of H         ~ */
//    maxU[1] = solve.maxdU(             d3x, d3y);                /* ~                                                           ~ */
//  }
//  else if (str_step.compare("correct") == 0) {
//    solve.partialsInXandY(run, tH,     d3x, d3y);                /* ~ d3x, d3y hold real-space partials in x and y of tH        ~ */
//  }
//  solve.bracket(run, BrKt, d1x, d1y, d3x, d3y);                  /* ~ calculate [phibar, H]                                     ~ */
//  for (unsigned k = 0; k < kstop; k++) { B1[k] = - BrKt[k]; }    /* ~ place result in B1                                        ~ */
//
//  if (model.compare("hall") == 0 ) {
//
//    if (str_step.compare("predict"     ) == 0) {
//      solve.partialsInXandY(run, V,      d3x, d3y);             /* ~ d3x, d3y hold real-space partials in x and y of V          ~ */
//      maxU[3] = solve.maxdU(             d3x, d3y);             /* ~                                                            ~ */
//    }
//    else if (str_step.compare("correct") == 0) {
//      solve.partialsInXandY(run, tV,     d3x, d3y);             /* ~ d3x, d3y hold real-space partials in x and y of tV         ~ */
//    }
//    solve.bracket(run, BrKt, d1x, d1y, d3x, d3y);               /* ~ calculate [phibar, V]                                      ~ */
//    for (unsigned k = 0; k < kstop; k++) { B3[k] = -BrKt[k]; }  /* ~ place result in B3                                         ~ */
//    solve.averageAcrossLayers( run, -1, d2x, d2y );             /* ~ calculate averages of Z_x & Z_y across adjacent layers     ~ */
//
//    solve.partialsInXandY(run, A,      d1x, d1y);               /* ~ d1x, d1y hold real-space partials in x and y of A          ~ */
//
////     maxu? = solve.maxdU(               d1x, d1y);            /* ~ might be interesting to do this calculation                ~ */
//    solve.bracket(run, BrKt, d1x, d1y, d2x, d2y);               /* ~ calculate [A, Zbar]                                        ~ */
//    for (unsigned k = 0; k < kstop; k++) { 
//      B1[k] = B1[k] - (rho  * BrKt[k]); 
//    }                                                           /* ~ B1 = - [phibar, H] - rho * [A, Zbar]                       ~ */
//    for (unsigned k = 0; k < kstop; k++) { 
//      B3[k] = B3[k] + (half * rtbeta * BrKt[k]); 
//    }                                                           /* ~ B3 = - [phibar, V] - (1/2) * sqrt{beta} * [A, Zbar]        ~ */
//
//  }
//  else if(model.compare("rmhd") == 0 ) {
//    solve.partialsInXandY(run, A,      d1x, d1y);               /* ~ d1x, d1y hold real-space partials in x and y of A          ~ */
//  }
//  solve.averageAcrossLayers( run, +1, d1x, d1y );             /* ~ calculate averages of A_x & A_y across adjacent layers       ~ */
//
//  if (model.compare("hall") == 0 ) { 
//
//    solve.averageAcrossLayers( run, +1, d3x, d3y );             /* ~ calculate averages of V_x & V_y across adjacent layers     ~ */
//    solve.bracket(run, BrKt, d1x, d1y, d3x, d3y);               /* ~ calculate [Abar, Vbar]                                     ~ */
//    for (unsigned k = 0; k < kstop; k++) { 
//      B2[k] = B2[k] + rtbeta * BrKt[k];                         /* ~ add result to B2                                           ~ */ 
//    }
//
//  }
//
//  solve.partialsInXandY(run, J,       d3x, d3y );             /* ~ d3x, d3y hold real-space partials in x and y of J          ~ */
//  solve.averageAcrossLayers( run, +1, d3x, d3y );             /* ~ calculate averages of J_x & J_y across adjacent layers     ~ */
//  solve.bracket(run, BrKt, d1x, d1y,  d3x, d3y );             /* ~ calculate [Abar, Jbar]                                     ~ */
//
//  for (unsigned k = 0; k < kstop; k++) {                      /* ~ B0 = -[phi, Omega] + [Abar, Jbar]                          ~ */
//    B0[k] = B0[k] + BrKt[k];
//  }
// 
//  if (model.compare("hall") == 0 ) { 
//
//    for (unsigned k = 0; k < kstop; k++) { 
//      B2[k] = B2[k] - (two * rho * BrKt[k]);                    /* ~ B2 = -[phi, Z]+sqrt{beta}*[Abar, Vbar]-2*rho*[Abar, Jbar]  ~ */
//    }
//
//  }
//
//  int n_flds;                                                 /* ~ set rank 0 maxU                                            ~ */
//  run.stack_data.fetch("iu3", &n_flds);
//
//  if (str_step.compare("predict") == 0) { MPI_Allreduce(MPI_IN_PLACE, &maxU.front(), n_flds, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);}
//
//  d1y.resize(0);
//  d2x.resize(0);
//  d2y.resize(0);
//  d3x.resize(0);
//  d3y.resize(0);
//
//  /* ~ bracket collection order ~ */ 
//
//  /* ~ -> [phi, Omega] ( for O equation    - B0    ) ... done with omega containers ~ (+2 - 1)     [2]  (phi    | Omega ) */
//  /* ~ -> [phi, Z    ] ( for Z equation    - B2    )                                ~ ( 1 + 1)     [2]  (phi, Z         ) */ 
//  /* ~ -> [phibar, H ] ( for H equation    - B1    ) ... done with H containers     ~ ( 2 + 1 - 1) [3]  (phi, Z |  H    ) */
//  /* ~ -> [phibar, V ] ( for V equation    - B3    ) ... done with phi containers   ~ ( 2 + 1 - 1) [3]  (V,   Z |  phi  ) */ 
//  /* ~ -> [A, Zbar   ] ( for H/V equations - B1/B3 ) ... done with Z containers     ~ ( 2 + 1 - 1) [3]  (V,   A |  Z    ) */ 
//  /* ~ -> [Abar, Vbar] ( for Z equation    - B2    ) ... done with V containers     ~ ( 2 - 1    ) [2]  (A      |  V    ) */ 
//
//  /* ~ -> [Abar, Jbar] ( for O equation    - B0    ) ... done with all containers   ~ ( 1 + 1 - 2) [2]  (       |  A, J ) */ 
// 
//
//  /* ~ Sequence: ~ */
//
//  /* ~  1.) calculate phi_x, phi_y, Omega_x, Omega_y ~ */ 
//  /* ~  2.) get maxb from phi (Omega?)               ~ */
//  /* ~  3.) calculate [phi, Omega]                   ~ */
//  /* ~  4.) place result of 3 in B0                  ~ */
//  /* ~  5.) get Z_x, Z_y                             ~ */
//  /* ~  6.) get maxb from Z                          ~ */
//  /* ~  7.) calculate [phi, Z]                       ~ */
//  /* ~  8.) place result of  7 in B2                 ~ */
//  /* ~  9.) get phibar_x, phibar_y, H_x, H_y         ~ */
//  /* ~ 10.) get maxb from H (A?)                     ~ */
//  /* ~ 11.) calculate [phibar,H]                     ~ */
//  /* ~ 12.) place result of 11 in B1                 ~ */
//  /* ~ 13.) get V_x, V_y                             ~ */
//  /* ~ 14.) get maxb from V                          ~ */
//  /* ~ 15.) calculate [phibar, V]                    ~ */
//  /* ~ 16.) place result of 15 in B3                 ~ */
//  /* ~ 17.) get Zbar_x, Zbar_y, A_x, A_y             ~ */
//  /* ~ 18.) calculate [A, Zbar]                      ~ */
//  /* ~ 19.) add result of 18 to B1 and B3            ~ */
//  /* ~ 20.) get Abar_x, Abar_y, Vbar_x, Vbar_y       ~ */
//  /* ~ 21.) calculate [Abar, Vbar]                   ~ */
//  /* ~ 22.) add result of 21 to B2                   ~ */
//  /* ~ 23.) get Jbar_x, Jbar_y                       ~ */
//  /* ~ 24.) calculate [Abar, Jbar]                   ~ */
//  /* ~ 25.) add result of 24 to B0 & B2              ~ */
//
//  /* ~ fxy sequence ~ */
//
//  /* ~ a. ) for field f - multiply f by i            ~ */
//  /* ~ b1.) for f_x multiply result of a by kx       ~ */
//  /* ~ b2.) then inverse transform result            ~ */
//  /* ~ c1.) for f_y multiply result a by ky          ~ */
//  /* ~ c2.) then inverse transform result            ~ */
//
///* ~ retain for testing                              ~ */
//
//// B0.assign(n1n2c, czero);
//// B1.assign(n1n2c, czero);
//// B2.assign(n1n2c, czero);
//// B3.assign(n1n2c, czero);

//}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void redhallmhd::setD( std::string str_step, stack& run, lcsolve& solve) {

//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  double dz;
//  run.stack_data.fetch("dz"   , &dz    );         /* ~ inter-layer width along z                        ~ */
//  double rho;
//  physics_data.fetch(  "rho"  , &rho   );         /* ~ gyro-radius parameter                            ~ */
//  double beta;                           
//  run.palette.fetch(   "beta" , &beta  );         /* ~ zero'th-order plasma-beta                        ~ */
//
//  double dzm1      = one / dz;
//  double rtbeta    = sqrt(beta);
//
//  ComplexArray& Z  = solve.U2;
//  ComplexArray& V  = solve.U3;
//
//  ComplexArray& tZ = solve.tU2;
//  ComplexArray& tV = solve.tU3;
//  
//  ComplexArray& D0 = solve.D0;                    /* ~ subtract across  i and i + 1                    ~ */
//  ComplexArray& D1 = solve.D1;                    /* ~ subtract across  i and i - 1                    ~ */
//  ComplexArray& D2 = solve.D2;                    /* ~ subtract across  i and i + 1                    ~ */
//  ComplexArray& D3 = solve.D3;                    /* ~ subtract across  i and i - 1                    ~ */
//
//  complex<double> deltaZ, deltaV, deltaJ, deltaP; /* ~ to aid in differentiating between the           ~ */
//                                                  /* ~ predictor and corrector cases                   ~ */
//
//  unsigned kdxp1,  kdxm1;                         /* ~ neighbor - layer indices                        ~ */
//  unsigned kstart, kstop;                         /* ~ limits on k looop                               ~ */
//
//  int n1n2c; 
//  run.stack_data.fetch( "n1n2c", &n1n2c );        /* ~ number of complex elements per layer            ~ */
//  int iu2;
//  run.stack_data.fetch( "iu2"  , &iu2   );        /* ~ number of layers                                ~ */
//
//  kstart    = n1n2c;                              /* ~ D's are calculated for layers 1,2,3,..n3        ~ */
//  kstop     = n1n2c * (iu2 - 1);                  /* ~ layer 1 needs layer 0 and layer n3 needs layer  ~ */
//                                                  /* ~ iu2 - 1                                         ~ */
//
//  std::string model;
//  run.palette.fetch("model", &model);
//
//  for (unsigned kdx = kstart; kdx < kstop; kdx++) {
//
//    kdxm1      = kdx - n1n2c;                       /* ~ adjacent lower layer index                      ~ */
//    kdxp1      = kdx + n1n2c;                       /* ~ adjacent upper layer index                      ~ */
//
//    deltaJ     = J[ kdxp1 ] - J[ kdx   ];           /* ~ P's and J's are updated every half-step         ~ */
//    deltaP     = P[ kdx   ] - P[ kdxm1 ];           /* ~ see updatePAJ. Note use of kdxp1 & kdxm1.       ~ */
//    
//    if ( model.compare("hall") == 0) {
//      if (     str_step.compare("predict") == 0) {
//      
//        deltaZ = Z[ kdx   ] - Z[ kdxm1 ];           /* ~ Z and V must retain un-updated values           ~ */
//        deltaV = V[ kdxp1 ] - V[ kdx   ];           /* ~ until corrector step. Note use of kdxm1 & kdxp1 ~ */
//
//      }
//      else if (str_step.compare("correct") == 0) {
//
//        deltaZ = tZ[ kdx   ] - tZ[ kdxm1 ];         /* ~ using results of predictor step here            ~ */
//        deltaV = tV[ kdxp1 ] - tV[ kdx   ];         /* ~ note use of kdxm1 & kdxp1                       ~ */
//
//      }
//    }
//    else if( model.compare("rmhd") == 0 ) {
//
//        deltaZ = zero;
//        deltaV = zero;
//
//    }
//
//    D0[kdx]    =  deltaJ                                          * dzm1;  /* ~ i.e. Delta J / Delta z                                        ~ */
//    D1[kdx]    = (          deltaP  - (       rho    *  deltaZ )) * dzm1;  /* ~ i.e. Delta F / Delta z, where F = phi - rhobar * Z            ~ */
//
//    if ( model.compare("hall") == 0) {
//
//      D2[kdx]  = ((rtbeta * deltaV) - (two  * rho    *  deltaJ )) * dzm1;  /* ~ i.e. Delta G / Delta z, where G = rtbeta * V - 2 * rhobar * J ~ */
//      D3[kdx]  =                      (half * rtbeta * deltaZ  )  * dzm1;  /* ~ i.e. 1/2 * sqrt{beta} * Delta Z delta z                       ~ */
//
//    }
//
///* ~ retained for testing ~ */
//
////  D0[kdx] =                 (J[kdxp1] - J[kdx  ] ) * dzm1;
////  D1[kdx] = (               (P[kdx  ] - P[kdxm1] ) -       rho * ( Z[kdx  ] - Z[kdxm1] ) ) * dzm1;  /* ~ F = phi - rhobar * Z  ~*/
////  D2[kdx] = (      rtbeta * (V[kdxp1] - V[kdx  ] ) - two * rho * ( J[kdxp1] - J[kdx  ] ) ) * dzm1;  /* ~ G = rtbeta * V - 2rhobar J ~ */
////  D3[kdx] = half * rtbeta * (Z[kdx  ] - Z[kdxm1] ) * dzm1;
//
//  }
//  
///* ~ retained for testing ~ */
//
////  D0.assign(nc, czero);
////  D1.assign(nc, czero);
////  D2.assign(nc, czero);
////  D3.assign(nc, czero);

//}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void redhallmhd::setAi( stack& run ) {

//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  int n1n2c;
//  run.stack_data.fetch( "n1n2c", &n1n2c );         /* ~ number of complex elements per layer                 ~ */
//  int iu2;
//  run.stack_data.fetch( "iu2",   &iu2 );           /* ~ number of layers                                     ~ */
//
//  ComplexArray& A0 = run.A0;
//  ComplexArray& A1 = run.A1;                     /* ~ I really only need this one                          ~ */
//  ComplexArray& A2 = run.A2;
//  ComplexArray& A3 = run.A3;
//
//  std::string model;
//  run.palette.fetch("model", &model);
//
//  A0.assign((n1n2c * iu2), czero);
//  A1.assign((n1n2c * iu2), czero);                 /* ~ to be set to  -eta * ssqd * k2 * A = -eta * ssqd * J ~ */
//
//  if (model.compare("hall") == 0 ) {
//
//    A2.assign((n1n2c * iu2), czero);
//    A3.assign((n1n2c * iu2), czero);
//
//    double eta; 
//    run.palette.fetch(  "eta", &eta);
//    double ssqd;
//    physics_data.fetch("ssqd", &ssqd);
//  
//    unsigned kstart  = 0;
//    unsigned kstop   = n1n2c * iu2;
//  
//    for (unsigned k  = kstart; k < kstop; k++) {
//  
//      A1[k]          = -( eta * ssqd * J[k] );
//  
//    }
//
//  }
//}

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

  double ssqd;                                   /* ~ parameter sigma^2 relating A to H           ~ */
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
    double dt;
    run.palette.fetch("dt", &dt);
    double dz;
    run.stack_data.fetch("dz", &dz);
 
    double q1, q2, qp;

    run.palette.fetch("q1", &q1);
    run.palette.fetch("q2", &q2);
    run.palette.fetch("qp", &qp);

    double dtvb;
    double dtr;

    int i_red;
    i_red = MPI_Reduce(&maxU[0], &glbMaxU[0], maxU.size(), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
    int i_brd;
    i_brd = MPI_Bcast(&glbMaxU[0], maxU.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (unsigned i_f = 0; i_f < maxU.size(); i_f++) { maxU[i_f] = glbMaxU[i_f]; }

    dtvb        = zero;

    for (unsigned i_f = 0; i_f < maxU.size(); i_f++) { dtvb = dtvb + sqrt(maxU[i_f]); }

    dtvb        = dtvb * sqrt( (double) (n1n2) );
    dtvb        = one / dtvb;
    dtr         = dtvb / dt;

    if (( dt < (q1 * dtvb) ) && ( dt <  (half * qp * dz) ) ) {

      dt = two * dt;
      run.palette.reset("dt", dt);

    }
    else if ( dt > (q2 * dtvb) ) {

      dt = half * dt;
      run.palette.reset("dt", dt);

    }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~ Destructor ~ */

redhallmhd:: ~redhallmhd() {

  fftw.fftwFinalize();

}
