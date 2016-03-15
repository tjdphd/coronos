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
 *        FILE: Implementation of class "stack"
 *
 * DESCRIPTION: Longcope-type staggered stack of slabs
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

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

#ifndef HAVE_CUDA_H

  allocUi( );
  allocAUX();
  initxyz( );

#endif

}

/* ~~~~~~~~~~~~~~~~ */
/* ~ initializers ~ */
/* ~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::init_stack_data() {                     /* ~ gather/infer information to be           ~ */
                                                    /* ~ included in stack_data container         ~ */

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* ~ incoming parameters from palette            ~ */
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  std::string model;                                /* ~ reduced mhd or hall-mrhd                 ~ */
  palette.fetch("model",&model);
  int p1;                                           /* ~ power of 2 specifying resolution in x    ~ */
  palette.fetch("p1"   ,&p1   );
  int p2;                                           /* ~ power of 2 specifying resolution in y    ~ */
  palette.fetch("p2"   ,&p2   );
  int p3;                                           /* ~ total number of layers in z              ~ */
  palette.fetch("p3"   ,&p3   );
  int np;                                           /* ~ np number of processes                   ~ */
  palette.fetch("np"   ,&np   );
  RealVar zl;                                       /* ~ length in z of computational domain      ~ */
  palette.fetch("zl"   ,&zl   );

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* ~ to be made parameters of stack_data         ~ */
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  std::string resolution;                           /* ~ full 3-d resolution string                ~ */

  int n1           = (int) pow(2.0, p1);            /* ~ number of x-coordinates in a layer        ~ */
  int n2           = (int) pow(2.0, p2);            /* ~ number of y-coordinates in a layer        ~ */
  int n3           =                p3 ;            /* ~ number of interior layers per process     ~ */

  int n1n2         = n1*n2;                         /* ~ total number points on a (real) layer     ~ */
  int n1n2c        = n1 * (((int)(half * n2)) + 1); /* ~ total number points on a (Fourier) layer  ~ */

  int iu1          = n1n2;                          /* ~ dimension of first index of U             ~ */
  int iu2          = n3+2;                          /* ~ dimension of second index of U            ~ */

                                                    /* ~ Ui's should have dimensions:              ~ */
                                                    /* ~ Ui[n1n2c * iu2]                           ~ */
                                                    /* ~                                           ~ */
                                                    /* ~ NOTE: relative to old codes:              ~ */
                                                    /* ~ pbot(:) = U0[0:(n1n2c - 1)]               ~ */
                                                    /* ~ atop(:) = U1[n1n2c*(iu2-1):(n1ncc*iu2)-1] ~ */

  int iu3;                                          /* ~ dimension of third index of U             ~ */
                                                    /* ~ counts number of fields in plasma model   ~ */

  if (model.compare("rmhd") == 0) iu3 = 2;          /* ~ fix number of field variables             ~ */
  if (model.compare("inhm") == 0) iu3 = 2;
  if (model.compare("hall") == 0) iu3 = 4;

  RealVar dz       = zl/((RealVar)(n3*np));         /* ~ layer separation in z                     ~ */

  int izres        = (int) (n3 * np)/zl;            /* ~ integer effective resolution in z         ~ */

  std::string xres = static_cast<std::ostringstream*>( &(std::ostringstream() << n1   ) ) -> str();
  std::string yres = static_cast<std::ostringstream*>( &(std::ostringstream() << n2   ) ) -> str();
  std::string zres = static_cast<std::ostringstream*>( &(std::ostringstream() << izres) ) -> str();

  if (xres.compare(yres) == 0 ) resolution.assign(xres + "_" + zres);
  else             resolution.assign(xres + "_" + yres + "_" + zres);

  std::string pname;                                /* ~ for containing parameter names            ~ */
  std::string padjust;                              /* ~ for specifying adjustability              ~ */
  
  padjust.assign("rfx"  );                          /* ~ assigning parameters that are run-fixed   ~ */

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

#ifndef HAVE_CUDA_H

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~ Allocators / De-allocators ~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::allocUi() {              /* ~ U is the input/output array for the fields ~ */
                                     /* ~ NOTE: U holds the real-space fields        ~ */
  int iu1, iu2, iu3;                 /* ~       and is read and written to dataframe ~ */
                                     /* ~       files for each process               ~ */

  stack_data.fetch("iu1", &iu1);
  stack_data.fetch("iu2", &iu2);
  stack_data.fetch("iu3", &iu3);

  U            = new RealVar**[iu1]; /* ~ allocate U dynamically using dimensions    ~ */

  for (int i   = 0; i< iu1; ++i) {   /* ~ determined in init_stack_data              ~ */
    U[i]       = new RealVar*[iu2];
    for (int j = 0; j < iu2; ++j) {
      U[i][j]  = new RealVar[iu3];
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

void stack::allocAUX() {              /* ~ AUX is the input/output array for auxiliary fields ~ */

  int iu1, iu2, iu3;

  stack_data.fetch("iu1", &iu1);      /* ~ coordinates per layer                              ~ */
  stack_data.fetch("iu2", &iu2);      /* ~ layers including bots and tops                     ~ */
  stack_data.fetch("iu3", &iu3);      /* ~ number of fields                                   ~ */

  int iaux1, iaux2, iaux3;

  iaux1         = iu1;
  iaux2         = iu2;
  iaux3         = iu3;

  AUX           = new RealVar**[iaux1];

  for ( int i = 0; i < iu1; ++i) {
    AUX[i]      = new RealVar*[iaux1]; 
    for ( int j = 0; j < iaux2; ++j) {
      AUX[i][j] = new RealVar[iaux3];
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

void stack::deallocAUX() {

  int iu1, iu2;

  stack_data.fetch("iu1", &iu1);
  stack_data.fetch("iu2", &iu2);

  int iaux1, iaux2;

  iaux1 = iu1;
  iaux2 = iu2;

  for (int i = 0; i< iaux1; ++i) {

    for (int j = 0; j < iaux2; ++j) delete [] AUX[i][j];
    delete [] AUX[i];

  }

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::initAUX() {


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::writeUData() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank  );

  int srun;
  palette.fetch("srun",    &srun  );

  std::string data_file;
  data_file                = getLastDataFilename(srun-1);

  const char *c_data_file;
  c_data_file              = data_file.c_str();

  std::ofstream ofs;
  ofs.open( c_data_file, std::ios::out );

  if ( ofs.good() ) {

    int iu3; 
    stack_data.fetch("iu3" , &iu3 );
    int n1; 
    stack_data.fetch("n1"  , &n1  );
    int n2; 
    stack_data.fetch("n2"  , &n2  );
    int n3; 
    stack_data.fetch("n3"  , &n3  );
    int n1n2;
    stack_data.fetch("n1n2", &n1n2);

    int n_slab_points      = n1n2 * iu3;
    int point_count        = 0;
    int slab_index         = 1;
    int from_col_maj_idx   = 0;
    int to_row_maj_idx     = 0;
    int i                  = 0;
    int j                  = 0;

    RealVar next_p;
    RealVar next_a;
    RealVar next_bz;
    RealVar next_vz;

    RealVar next_o;
    RealVar next_j;
    
    while ( slab_index < n3 + 1 ) {

      ++point_count;
      next_p               =   U[to_row_maj_idx][slab_index][0];
      next_o               = AUX[to_row_maj_idx][slab_index][0];
      ++point_count;
      next_a               =   U[to_row_maj_idx][slab_index][1];
      next_j               = AUX[to_row_maj_idx][slab_index][1];

      if(iu3 > 2) {

        ++point_count;
        next_bz            =   U[to_row_maj_idx][slab_index][2];
        ++point_count;
        next_vz            =   U[to_row_maj_idx][slab_index][3];

      }

      ofs   << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_p << " ";
      ofs   << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_a << " ";
      if (iu3 < 3)  {
        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_o << " ";
        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_j << " ";
      }

      if (iu3 > 2)  {

        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_bz << " ";
        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_vz << " ";
        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_o  << " ";
        ofs << std::setw(24) << std::right << std::setprecision(16) << std::scientific << next_j  << " ";

      }

      ofs   << std::endl;
     
      if (from_col_maj_idx < n1n2) {

        ++from_col_maj_idx;
        if (from_col_maj_idx % n2 != 0) ++j;
        else {
                 j         = 0;
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

void stack::reportEnergyQs( double t_cur ) {

  int rank;
  run_data.fetch("rank",      &rank   );
  std::string prefix;
  palette.fetch("prefix",     &prefix );
  std::string run_label;
  palette.fetch("run_label",  &run_label);
  std::string res_str;
  stack_data.fetch("res_str", &res_str);

  std::string energy_data_file;
  energy_data_file = prefix + "_" + res_str + ".o" + run_label;

  const char *c_energy_data_file;
  c_energy_data_file = energy_data_file.c_str();
  std::ofstream ofs;

  int esize;
  esize = EnergyQs.size();

  if (rank == 0 ) {

    ofs.open( c_energy_data_file, std::ios::out | std::ios::app );
    if (ofs.good() ) {
      ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << t_cur << " ";
      for (unsigned k = 0; k < esize; k++){
        ofs << std::setw(24) << std::right << std::setprecision(12) << std::scientific << EnergyQs[k] << " ";
      }
      ofs << std::endl;
    }
    else {
      std::cout << "reportEnergyQs: ofs is not good! " << std::endl;
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

  data_file             = prefix + "_" + res_str + "." + rnk_str + ".o" + run_label + srn_str;

  return data_file;

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  std::string stack::getNextDataFilename() {

  std::string data_file = "not_this_one";

  std::string prefix;
  palette.fetch(   "prefix" , &prefix );

  std::string run_label;
  palette.fetch("run_label",  &run_label);

  int rank; 
  run_data.fetch(  "rank"   , &rank   );
  std::string rnk_str;
  rnk_str           = static_cast<std::ostringstream*>( &(std::ostringstream() << rank ) ) -> str();

  int srun;
  palette.fetch(   "srun"   , &srun   );
  std::string srn_str;
  srn_str           = static_cast<std::ostringstream*>( &(std::ostringstream() << srun ) ) -> str();

  std::string res_str;
  stack_data.fetch("res_str", &res_str);

  int rnk_len       = rnk_str.length();

  switch(rnk_len) {

  case(1) : rnk_str = "00" + rnk_str;
            break;
  case(2) : rnk_str =  "0" + rnk_str;
            break;
  default : std::cout << "this can't be right " << std::endl;

  }

  data_file         = prefix + "_" + res_str + "." + rnk_str + ".o" + run_label + srn_str;

  return data_file;

  }

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void stack::initxyz() {                     /* ~ Calculate x- and y-coordinates of layers ~ */

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1; 
  stack_data.fetch("n1",    &n1 );
  x.reserve(n1);

  RealVar dx = one / ((RealVar) n1);


  RealVar next_x;
  for (int i = 0; i < n1; ++i) {

    next_x   = ( (RealVar) i) * dx; 
    x.push_back(next_x);

  }


  int n2;
  stack_data.fetch("n2",    &n2 );
  y.reserve(n2);

  RealVar dy = one / ((RealVar) n2);

  RealVar next_y;

  for (int j = 0; j < n2; ++j) {

    next_y   = ( (RealVar) j) * dy;
    y.push_back(next_y);

  }


  int np;
  palette.fetch(   "np"  , &np );

  int p3;
  palette.fetch(   "p3"  , &p3 );
  int    iu2;
  stack_data.fetch("iu2" , &iu2 );
  RealVar dz;
  stack_data.fetch("dz"  , &dz );
  RealVar zl;
  palette.fetch(   "zl"  , &zl );

//  z.reserve(iu2);
  z.reserve(p3);

  RealVar next_z;

  for (int i = 0; i < p3-1; ++i) {

    next_z   =  ((RealVar) (rank)) * (zl / ((RealVar) (np))) + (((RealVar) i) * dz);
    z.push_back(next_z);
 
  }
}

#endif

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

  std::string run_label;
  palette.fetch("run_label",  &run_label);

  if (rank == 0) {palette.report(prefix + '_' + res_str, run_label, srun ); }

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~ */
/* ~ Destructor  ~ */
/* ~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

stack::~stack() {
#ifndef HAVE_CUDA_H
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

  deallocUi();
#endif
}

