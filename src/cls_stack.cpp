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

  int srun;
  palette.fetch("srun", &srun);
  writeParameters(srun - 1);

  allocUi();
  initxyz();

//  kInit();
//  rtInit();

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//  if (rank == 0) {

//  int n1;                                        /* ~ number of coordinates in x                               ~ */
//  stack_data.fetch("n1",    &n1);
//  int n2; 
//  stack_data.fetch("n2",    &n2);                /* ~ number of coordinates in y                               ~ */
//  int n2h    = (((int)(half*n2)) + 1);
//  int ndx;
 
//   std::cout << "kx: " << std::endl << std::endl;
 
//    for (int i = 0; i < n1; ++i) {
//      for (int j = 0; j < n2/2 + 1; ++j) {
// 
//        ndx = (i * n2h) + j;
//        if (rt[ndx] == zero){
//          std::cout << std::setw(3)  << std::right  << "kx(" << std::setw(4) << std::right << i << ","; 
//          std::cout << std::setw(4)  << std::right  << j     << std::setw(4) << std::right << ") = ";
//          std::cout << std::setw(10) << std::setprecision(4)                 << kx[ndx]/two_pi << " ";
//
//          std::cout << std::setw(3)  << std::right  << "ky(" << std::setw(4) << std::right << i << ","; 
//          std::cout << std::setw(4)  << std:: right << j     << std::setw(4) << std::right << ") = ";
//          std::cout << std::setw(10) << std::setprecision(4)                 << ky[ndx]/two_pi << " ";
//
//          std::cout << std::setw(3)  << std::right << "k2(" << std::setw(4)  << std::right << i << ","; 
//          std::cout << std::setw(4)  << std::right << j     << std::setw(4)  << std::right << ") = ";
//          std::cout << std::setw(10) << std::setprecision(4)                 << k2[ndx]/(two_pi*two_pi) << " ";
//
//          std::cout << std::setw(3)  << std::right << "rt(" << std::setw(4)  << std::right << i << ","; 
//          std::cout << std::setw(4)  << std::right << j     << std::setw(4)  << std::right << ") = ";
//          std::cout << std::setw(10) << std::setprecision(4)                 << rt[ndx] << std::endl;
//       }
// 
//      }
//    }
//  }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

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

  double dz        = zl/double((n3*np));         /* ~ layer separation in z                     ~ */

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

//  std::cout << "iu1 = " << iu1 << std::endl;
//  std::cout << "iu2 = " << iu2 << std::endl;
//  std::cout << "iu3 = " << iu3 << std::endl;


  U = new double**[iu1];                         /* ~ allocate U dynamically using dimensions    ~ */

  for (int i = 0; i< iu1; ++i) {                 /* ~ determined in init_stack_data              ~ */
    U[i] = new double*[iu2];
    for (int j = 0; j < iu2; ++j) {
      U[i][j] = new double[iu3];

    }
  }

  int n1n2c;
  stack_data.fetch("n1n2c", &n1n2c);
  std::string model;
  stack_data.fetch("model", &model);
  U0.reserve(n1n2c * iu2);
  U1.reserve(n1n2c * iu2);
  if (model.compare("hall") == 0) { 
    U2.reserve(n1n2c * iu2);
    U3.reserve(n1n2c * iu2);
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

//void stack::rtInit( ) {                          /* ~ initialize de-aliasing array              ~ */
//
//  int n1;                                        /* ~ number of x-coordinates                   ~ */
//  stack_data.fetch("n1",    &n1);
//  int n2;
//  stack_data.fetch("n2",    &n2);                /* ~ number of y-coordinates                   ~ */
//  int n1n2c;
//  stack_data.fetch("n1n2c", &n1n2c);             /* ~ number of Fourier space points per layer  ~ */
//
//  RealArray::size_type nc;                       /* ~ a vector size-type version of n1n2c       ~ */
//  nc  = n1n2c;
//
//  int n1h = (((int)(half*n1)) + 1);
//  int n2h = (((int)(half*n2)) + 1);
//
//// double threshold = two_thirds * sqrt(two * ((( half * n1 - one) * two_pi) * ((half * n1 - one) * two_pi)));
//
//  double kx_max = two_pi * half * n1;
//  double ky_max = two_pi * half * n2;
//  double k_max;
//
//  if (kx_max >= ky_max) { k_max = ky_max;}
//  else                  { k_max = kx_max;}
//
//  double threshold = two_thirds * k_max;
//
////  std::cout << "rtInit: kx_max  = "   << kx_max /two_pi    << std::endl;
////  std::cout << "rtInit: ky_max  = "   << ky_max /two_pi    << std::endl;
////  std::cout << "rtInit: k_max   = "   << k_max /two_pi     << std::endl;
////  std::cout << "rtInit: threshold = " << threshold /two_pi << std::endl;
//
////; std::cout << "rtInit: threshold = " << threshold / (two_pi*two_pi) << std::endl;
//
//  rt.reserve(nc);                                /* ~ create space in rt                        ~ */
//
//  RealArray::size_type ndx;
//
//  for (unsigned i = 0; i < n1; i++) {            /* ~ initialization                            ~ */
//    for (unsigned j = 0; j < n1h; j++) {
//
//      ndx = i * n2h + j;
//
//      if (  ky[ndx] == zero) { 
//        if ( abs(kx[ndx]) < k_max) {
//          rt[ndx] = one;
//        }
//        else {
//          rt[ndx] = zero;
//        }
//      }
//      else if ( abs(ky[ndx]) >  threshold  || abs(kx[ndx]) >  threshold  ) {
//        rt[ndx] = zero;
//      }
//      else {
//        rt[ndx] = one;
//      }
//
////      if ( sqrt(std::abs(k2[ndx])) < threshold ) rt[ndx] = one;
////
////      else                                       rt[ndx] = zero;
//
////    }
//  }
// }
//}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void stack::rtFree() {   /* ~ this is a bit excessive no? maybe just put in destructor or something ? ~ */
//
// rt.resize(0);
//
//}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//void stack::kInit() {                            /* ~ initialize the wave number arrays kx, ky, k2, and inv_k2 ~ */
//
//  int n1;                                        /* ~ number of coordinates in x                               ~ */
//  stack_data.fetch("n1",    &n1);
//  int n2; 
//  stack_data.fetch("n2",    &n2);                /* ~ number of coordinates in y                               ~ */
//  int nc;
//  stack_data.fetch("n1n2c", &nc);                /* ~ number of Fourier space points per layer                 ~ */
//
//  int n1h    = (((int)(half*n1)) + 1);
//  int n2h    = (((int)(half*n2)) + 1);
//
//  RealArray kp;                                  /* ~ temporary definition of convenience                      ~ */
//
//  kp.reserve(n1);
//
//      kx.reserve(nc);                            /* ~ these are all lcstack members and carry the wave-number  ~ */
//      ky.reserve(nc);                            /* ~ information                                              ~ */
//      k2.reserve(nc);
//  inv_k2.reserve(nc);
//
//                                                 /* ~ the initialization as done in gpu port of old code       ~ */
//  int ndx;
//
//  for (int i = 0;   i < n1h; i++ ) kp[i] = ((double) i) * two_pi;
//  for (int i = n1h; i < n1;  i++ ) kp[i] = -kp[n1 - i];
//
//  for (int i = 0; i < n1; i++) {
//    for (int j = 0; j < n2h; j++) {
//   
//      ndx     = (i * n2h) + j;
//      kx[ndx] = kp[j];
//      ky[ndx] = kp[i];
//    }
//  }
//
//  kp.resize(0);
//
//  for (int i = 0; i < n1; i++) {
//    for (int j = 0; j < n1h; j++) {
//
//      ndx     = (i * n2h) + j;
//      k2[ndx] = (kx[ndx] * kx[ndx]) +  (ky[ndx] * ky[ndx]);
//
//      if (std::abs(k2[ndx]) > teensy) inv_k2[ndx] = one/k2[ndx];
//      else inv_k2[ndx] = zero;
//
//    }
//  }
//
///* ~ save for debugging ~ */
//
////int rank;
////MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 
////if (rank == 0) {
// 
////  std::cout << "kx: " << std::endl << std::endl;
// 
////  for (int i = 0; i < n1; ++i) {
////    for (int j = 0; j < n1/2 + 1; ++j) {
// 
////      ndx = (i * n2h) + j;
////      std::cout << "kx[" << i << "," << j << "] = "
////      std::cout << std::setw(10) << std::right << std::setprecision(4) << std::scientific << kx[ndx]/two_pi << std::endl;
// 
////    }
////  }
////}
//  
//
///* ~ save for debugging ~ */
//
//}

///* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//
//void stack::kFree() {                            /* ~ also possibly a bit excessive ~ */
//
//    kx.resize(0);
//    ky.resize(0);
//    k2.resize(0);
//    inv_k2.resize(0);
//}


  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* ~ Initialization            ~ */
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::writeUData() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank  );

  int    srun;
  palette.fetch("srun",    &srun  );

  std::cout << "writeUData: srun = " << srun << std::endl;

  std::string data_file;
  data_file   = getLastDataFilename(srun-1);

  std::cout << "writeUData: writing to file - " << data_file << std::endl;
  const char  *c_data_file;
  c_data_file = data_file.c_str();

  std::ofstream ofs;
  ofs.open( c_data_file, std::ios::out );

  if ( ofs.good() ) {

//    InputOutputArray& U = un.U;
    
    int iu3; 
    stack_data.fetch("iu3",  &iu3);
    int n1; 
    stack_data.fetch("n1",   &n1);
    int n2; 
    stack_data.fetch("n2",   &n2);
    int n3; 
    stack_data.fetch("n3",   &n3);
    int n1n2;
    stack_data.fetch("n1n2", &n1n2);

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
    
    while ( slab_index < n3 + 1 ) {

//    U[to_row_maj_idx][slab_index][0] = next_p;
//    U[to_row_maj_idx][slab_index][1] = next_a;

      ++point_count;
      next_p = U[to_row_maj_idx][slab_index][0];
      ++point_count;
      next_a = U[to_row_maj_idx][slab_index][1];

      if(iu3 > 2) {

//        U[to_row_maj_idx][slab_index][2] = next_bz;
//        U[to_row_maj_idx][slab_index][3] = next_vz;

        ++point_count;
        next_bz = U[to_row_maj_idx][slab_index][2];
        ++point_count;
        next_vz = U[to_row_maj_idx][slab_index][3];

      }

//      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_p << " ";
//      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_a << " ";
//
      ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_p << " ";
      ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_a << " ";

      if (iu3 > 2)  {

//      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_bz << " ";
//      ofs << std::setw(19) << std::right << std::setprecision(11) << std::scientific << next_vz << " ";

        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_bz << " ";
        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_vz << " ";

      }

      ofs << std::endl;
     
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

    ofs.close();

  }
}
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

std::string stack::getLastDataFilename(int srun) {

  int rank;
  run_data.fetch("rank",      &rank   );

  std::string prefix;
  palette.fetch("prefix",     &prefix );

  std::string run_label;
  palette.fetch("run_label",  &run_label);

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

//  data_file             = prefix + "_" + res_str + "." + rnk_str + ".ots" + srn_str;
  data_file             = prefix + "_" + res_str + "." + rnk_str + ".o" + run_label + srn_str;

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

void stack::writeParameters() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  if (rank == 0) {palette.report("coronos.in"); }

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::writeParameters(int srun) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  std::string prefix;
  palette.fetch(   "prefix",  &prefix );

  std::string res_str;
  stack_data.fetch("res_str", &res_str);

  if (rank == 0) {palette.report(prefix + '_' + res_str, srun ); }

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

//  if (rank == 1 ) {
//    for (int i = 0; i < iu2 - 1; ++i) {
//
//      std::cout << "Z[" <<  i <<  "] = " <<  z[i]  << endl;
//    }
//  }

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

  U0.resize(0);
  U1.resize(0);
  U2.resize(0);
  U3.resize(0);

  tU0.resize(0);
  tU1.resize(0);
  tU2.resize(0);
  tU3.resize(0);

//  rtFree();
//  kFree();
  deallocUi();

}
