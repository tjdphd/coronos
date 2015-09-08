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

  redhallmhd physics ( run );

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

    passAdjacentLayers("predict", run );

    physics.updatePAJ( "predict", run );            /* ~ P, A, and J contain un-updated/corrector-updated values ~ */
    physics.applyBC(              run );

    setS(              "predict", run, physics );   /* ~ set predictor S's                                       ~ */
    setB(              "predict", run, physics );   /* ~ set predictor Brackets                                  ~ */
    setD(              "predict", run, physics );   /* ~ set predictor finite differences                        ~ */
    setAi(                        run, physics );   /* ~ set predictor A's                                       ~ */
    Step(              "predict", run );            /* ~ execute predictor update                                ~ */

    passAdjacentLayers("correct", run );
    physics.updatePAJ( "correct", run );            /* ~ P, A, and J now contain predictor-updated values        ~ */

    setS(              "correct", run, physics );   /* ~ set corrector S's                                       ~ */
    setB(              "correct", run, physics );   /* ~ set corrector Brackets                                  ~ */
    setD(              "correct", run, physics );   /* ~ set corrector finite differences                        ~ */
    setAi(                        run, physics );   /* ~ set corrector A's                                       ~ */
    Step(              "correct", run );            /* ~ execute corrector update                                ~ */

  run.palette.fetch("dt", &dt);
  t_cur       = t_cur + dt;

  physics.updateTimeInc(          run );

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

  physics.updatePAJ("predict", run        );   /* ~ P, A, and J contain final corrector-updated values ~ */

  run.palette.reset(   "tstart", t_cur    );

  std::string prefix;
  run.palette.fetch(   "prefix", &prefix  );

  std::string res_str;
  run.stack_data.fetch("res_str", &res_str);

  int srun;
  run.palette.fetch(   "srun"   , &srun   );
  

  if (rank == 0) { run.palette.report(prefix + '_' +  res_str, srun); }

  physics.writeUData( run );

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

  ComplexArray& tO    = run.tU0;                   /* ~ for corrector case                          ~ */
  ComplexArray& tH    = run.tU1;                   /* ~ results from predictor step are transferred ~ */
  ComplexArray& tZ    = run.tU2;
  ComplexArray& tV    = run.tU3;

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

void lcsolve::setS( std::string str_step, stack& run, redhallmhd& physics ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string model;
  run.palette.fetch("model", &model);

  double  s0,  s1,  s2,   s3;                    /* ~ implicit fraction s - parameters                 ~ */
  run.palette.fetch(  "s0",   &s0 );
  run.palette.fetch(  "s1",   &s1 );

  if (model.compare("hall") == 0) {

    run.palette.fetch("s2",   &s2 );
    run.palette.fetch("s3",   &s3 );

  }

  double qs0, qs1, qs2,  qs3;                    /* ~ dissipation 'q' - parameters see documentation   ~ */

  physics.physics_data.fetch("qs0", &qs0);
  physics.physics_data.fetch("qs1", &qs1);

  if (model.compare("hall") == 0) {

    physics.physics_data.fetch("qs2", &qs2);
    physics.physics_data.fetch("qs3", &qs3);

  }

  double ge0, ge1, ge2,  ge3;
  double gi0, gi1, gi2,  gi3;

  double  dt;                                    /* ~ the current time increment                       ~ */
  run.palette.fetch("dt",   &dt );

  double  pfrac;                                 /* ~ fraction of dt to use in predictor step          ~ */
  run.palette.fetch("pfrac", &pfrac);

  if ( str_step.compare("predict") == 0 ) {      /* ~ use partial time-step for predictor case         ~ */
    dt           = dt * pfrac;                   /* ~ dt is local to setS so no harm done here         ~ */
  }

  ge0            = (one - s0) * qs0 * dt;        /* ~ th g^(ex)'s  see documentation                   ~ */
  ge1            = (one - s1) * qs1 * dt;

  if (model.compare("hall") == 0) {

    ge2          = (one - s2) * qs2 * dt;
    ge3          = (one - s3) * qs3 * dt;

  }

  gi0            =        s0  * qs0 * dt;        /* ~ th g^(im)'s  see documentation                   ~ */
  gi1            =        s1  * qs1 * dt;

  if (model.compare("hall") == 0) {

    gi2          =        s2  * qs2 * dt;
    gi3          =        s3  * qs3 * dt;

  }

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

  RealArray& k2  = run.k2;                       /* ~ square magnitude of k-space vectors              ~ */

  int n1n2c;                                     /* ~ number of complex elements in layer              ~ */
  run.stack_data.fetch( "n1n2c", &n1n2c );

  RealArray::size_type nc = SE0.capacity();      /* ~ S's should already by sized. This is a check     ~ */
  assert (nc     == n1n2c);

  for (unsigned k = 0; k < n1n2c; k++) {         /* ~ there are only as many S-elements as k2 elements ~ */

    SE0[k]       = ( one - (ge0 * k2[k]));       /* ~ S's are initialized. See documentation           ~ */
    SE1[k]       = ( one - (ge1 * k2[k]));  

    if (model.compare("hall") == 0) {

      SE2[k]     = ( one - (ge2 * k2[k]));  
      SE3[k]     = ( one - (ge3 * k2[k]));  

    }

    SI0[k]       = one / ( one + (gi0 * k2[k]));
    SI1[k]       = one / ( one + (gi1 * k2[k]));

    if (model.compare("hall") == 0) {

      SI2[k]     = one / ( one + (gi2 * k2[k]));
      SI3[k]     = one / ( one + (gi3 * k2[k]));

    }
  }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

//  if (str_step.compare("correct") == 0 && rank == 1) {

//  unsigned kstart = n1n2c;
//  unsigned kstop  = (2 * kstart) - 1;

//   std::cout << "setS: one   = " << one    << "\n\n";

//   std::cout << "setS: dt    = " << dt     << "\n\n";
//   std::cout << "setS: pfrac = " << pfrac  << "\n\n";

//   std::cout << "setS: s0    = " << s0     << "\n\n";
//   std::cout << "setS: s1    = " << s1     << "\n\n";

//   std::cout << "setS: qs0   = " << qs0    << "\n\n";
//   std::cout << "setS: qs1   = " << qs1    << "\n\n";

//   std::cout << "setS: gi0   = " << gi0    << "\n\n";
//   std::cout << "setS: gi1   = " << gi1    << "\n\n";

//   unsigned kstart = 0;
//   unsigned kstop  = kstart + n1n2c;
//    
//   for (unsigned k = kstart; k < kstop; k++) {

//     std::cout << std::setw(24) << std::right << std::setprecision(16) << std::scientific << SI1[k] << std::endl;

 //   }
//  }

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

//  SE0.assign(nc, one);                         /* ~ retained for testing                             ~ */
//  SE1.assign(nc, one);
//  SE2.assign(nc, one);
//  SE3.assign(nc, one);
//
//  SI0.assign(nc, one);
//  SI1.assign(nc, one);
//  SI2.assign(nc, one);
//  SI3.assign(nc, one);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void lcsolve::setB( std::string str_step, stack& run, redhallmhd& physics ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2;
  run.stack_data.fetch( "n1n2" ,&n1n2  );
  int n1n2c; 
  run.stack_data.fetch( "n1n2c",&n1n2c );
  int iu2;
  run.stack_data.fetch( "iu2"  ,&iu2   );

  unsigned kstop = n1n2c * iu2;

  double rho;
  physics.physics_data.fetch(  "rho", &rho );               /* ~ gyro-radius parameter                                      ~ */
  double beta, rtbeta;                           
  run.palette.fetch( "beta" , &beta  );                   /* ~ zero'th-order plasma-beta                                  ~ */
  rtbeta = sqrt(beta);

  ComplexArray BrKt;
  BrKt.reserve(n1n2c * iu2);

  RealArray d1x, d1y;
  RealArray d2x, d2y;
  RealArray d3x, d3y;

  d1x.reserve(n1n2 * iu2);
  d1y.reserve(n1n2 * iu2);

  d2x.reserve(n1n2 * iu2);
  d2y.reserve(n1n2 * iu2);

  d3x.reserve(n1n2 * iu2);
  d3y.reserve(n1n2 * iu2);

//  ComplexArray& B0 = solve.B0;
//  ComplexArray& B1 = solve.B1;
//  ComplexArray& B2 = solve.B2;
//  ComplexArray& B3 = solve.B3;

  ComplexArray& O  = run.U0;
  ComplexArray& H  = run.U1;
  ComplexArray& Z  = run.U2;
  ComplexArray& V  = run.U3;

  ComplexArray& tO = run.tU0;
  ComplexArray& tH = run.tU1;
  ComplexArray& tZ = run.tU2;
  ComplexArray& tV = run.tU3;

  ComplexArray& P  = physics.P;
  ComplexArray& A  = physics.A;
  ComplexArray& J  = physics.J;

  RealArray& maxU  = physics.maxU;

  std::string model;
  run.palette.fetch("model", &model);

  partialsInXandY( run, physics, P, d1x, d1y);                /* ~ d1x, d1y hold real-space partials in x and y of P         ~ */

//  if (str_step.compare("predict") == 0 && rank == 0) {
//    for( unsigned k = n1n2; k < 2*n1n2; ++k ) {

 //     std::cout << std::setw(16) << std::right << std::setprecision(8) << std::scientific << d1x[k] << std::endl;

//    }
//  }

if (str_step.compare("predict"     ) == 0) {
  partialsInXandY( run, physics, O, d2x, d2y);                /* ~ d2x, d2y hold real-space partials in x and y of O         ~ */
  maxU[0] = maxdU(         d2x, d2y);                /* ~                                                           ~ */
}
else if (str_step.compare("correct") == 0) {
  partialsInXandY( run, physics, tO,     d2x, d2y);  /* ~ d2x, d2y hold real-space partials in x and y of tO        ~ */
}
bracket( run, physics, BrKt, d1x, d1y, d2x, d2y);                  /* ~ calculate [phi, Omega]                                    ~ */
for (unsigned k = 0; k < kstop; k++) { B0[k] = - BrKt[k]; } /* ~ place result in B0                                 ~ */

if (model.compare("hall") == 0 ) {

   if (str_step.compare("predict"     ) == 0) {
     partialsInXandY( run, physics, Z, d2x, d2y);   /* ~ d2x, d2y hold real-space partials in x and y of Z         ~ */
     maxU[2] = maxdU(        d2x, d2y);             /* ~                                                           ~ */
   }
   else if (str_step.compare("correct") == 0) {
     partialsInXandY( run, physics, tZ, d2x, d2y);  /* ~ d2x, d2y hold real-space partials in x and y of tZ        ~ */
   }
   bracket( run, physics, BrKt, d1x, d1y, d2x, d2y);                     /* ~ calculate [phi, Z]                                        ~ */
   for (unsigned k = 0; k < kstop; k++) { B2[k] = - BrKt[k]; } /* ~ place result in B2                                  ~ */
}

averageAcrossLayers( run, -1, d1x, d1y );                       /* ~ calculate averages of phi_x & phi_y across adj't layers   ~ */
if (str_step.compare("predict"     ) == 0) {
  partialsInXandY( run, physics, H, d3x, d3y);            /* ~ d3x, d3y hold real-space partials in x and y of H         ~ */
  maxU[1] = maxdU(         d3x, d3y);                     /* ~                                                           ~ */
}
else if (str_step.compare("correct") == 0) {
  partialsInXandY( run, physics, tH, d3x, d3y);           /* ~ d3x, d3y hold real-space partials in x and y of tH        ~ */
}
bracket( run, physics, BrKt, d1x, d1y, d3x, d3y);                       /* ~ calculate [phibar, H]                                     ~ */
for (unsigned k = 0; k < kstop; k++) { B1[k] = - BrKt[k]; }    /* ~ place result in B1                                        ~ */

if (model.compare("hall") == 0 ) {

  if (str_step.compare("predict"     ) == 0) {
    partialsInXandY( run, physics, V, d3x, d3y);             /* ~ d3x, d3y hold real-space partials in x and y of V          ~ */
    maxU[3] = maxdU(             d3x, d3y);             /* ~                                                            ~ */
  }
  else if (str_step.compare("correct") == 0) {
    partialsInXandY( run, physics, tV, d3x, d3y);             /* ~ d3x, d3y hold real-space partials in x and y of tV         ~ */
  }
  bracket( run, physics, BrKt, d1x, d1y, d3x, d3y);               /* ~ calculate [phibar, V]                                      ~ */
  for (unsigned k = 0; k < kstop; k++) { B3[k] = -BrKt[k]; }  /* ~ place result in B3                                         ~ */
  averageAcrossLayers( run,  -1, d2x, d2y );             /* ~ calculate averages of Z_x & Z_y across adjacent layers     ~ */

  partialsInXandY( run, physics, A, d1x, d1y);               /* ~ d1x, d1y hold real-space partials in x and y of A          ~ */

//     maxu? = maxdU(               d1x, d1y);            /* ~ might be interesting to do this calculation                ~ */

  bracket( run, physics, BrKt, d1x, d1y, d2x, d2y);               /* ~ calculate [A, Zbar]                                        ~ */
  for (unsigned k = 0; k < kstop; k++) { 
    B1[k] = B1[k] - (rho  * BrKt[k]); 
  }                                                           /* ~ B1 = - [phibar, H] - rho * [A, Zbar]                       ~ */
  for (unsigned k = 0; k < kstop; k++) { 
    B3[k] = B3[k] + (half * rtbeta * BrKt[k]); 
  }                                                           /* ~ B3 = - [phibar, V] - (1/2) * sqrt{beta} * [A, Zbar]        ~ */

}
else if(model.compare("rmhd") == 0 ) {
  partialsInXandY( run, physics, A, d1x, d1y);               /* ~ d1x, d1y hold real-space partials in x and y of A          ~ */
}
averageAcrossLayers( run,  +1, d1x, d1y );             /* ~ calculate averages of A_x & A_y across adjacent layers       ~ */

if (model.compare("hall") == 0 ) { 

  averageAcrossLayers( run,  +1, d3x, d3y );             /* ~ calculate averages of V_x & V_y across adjacent layers     ~ */
  bracket( run, physics, BrKt, d1x, d1y, d3x, d3y);               /* ~ calculate [Abar, Vbar]                                     ~ */
  for (unsigned k = 0; k < kstop; k++) { 
    B2[k] = B2[k] + rtbeta * BrKt[k];                         /* ~ add result to B2                                           ~ */ 
  }

}

partialsInXandY( run, physics, J, d3x, d3y );             /* ~ d3x, d3y hold real-space partials in x and y of J          ~ */
averageAcrossLayers( run,  +1, d3x, d3y );             /* ~ calculate averages of J_x & J_y across adjacent layers     ~ */
bracket( run, physics, BrKt, d1x, d1y,  d3x, d3y );             /* ~ calculate [Abar, Jbar]                                     ~ */

for (unsigned k = 0; k < kstop; k++) {                      /* ~ B0 = -[phi, Omega] + [Abar, Jbar]                          ~ */
  B0[k] = B0[k] + BrKt[k];
}

if (model.compare("hall") == 0 ) { 

  for (unsigned k = 0; k < kstop; k++) { 
    B2[k] = B2[k] - (two * rho * BrKt[k]);                    /* ~ B2 = -[phi, Z]+sqrt{beta}*[Abar, Vbar]-2*rho*[Abar, Jbar]  ~ */
  }

}

int n_flds;                                                 /* ~ set rank 0 maxU                                            ~ */
run.stack_data.fetch("iu3", &n_flds);

if (str_step.compare("predict") == 0) { MPI_Allreduce(MPI_IN_PLACE, &maxU.front(), n_flds, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);}

d1y.resize(0);
d2x.resize(0);
d2y.resize(0);
d3x.resize(0);
d3y.resize(0);

  /* ~ bracket collection order ~ */ 

  /* ~ -> [phi, Omega] ( for O equation    - B0    ) ... done with omega containers ~ (+2 - 1)     [2]  (phi    | Omega ) */
  /* ~ -> [phi, Z    ] ( for Z equation    - B2    )                                ~ ( 1 + 1)     [2]  (phi, Z         ) */ 
  /* ~ -> [phibar, H ] ( for H equation    - B1    ) ... done with H containers     ~ ( 2 + 1 - 1) [3]  (phi, Z |  H    ) */
  /* ~ -> [phibar, V ] ( for V equation    - B3    ) ... done with phi containers   ~ ( 2 + 1 - 1) [3]  (V,   Z |  phi  ) */ 
  /* ~ -> [A, Zbar   ] ( for H/V equations - B1/B3 ) ... done with Z containers     ~ ( 2 + 1 - 1) [3]  (V,   A |  Z    ) */ 
  /* ~ -> [Abar, Vbar] ( for Z equation    - B2    ) ... done with V containers     ~ ( 2 - 1    ) [2]  (A      |  V    ) */ 

  /* ~ -> [Abar, Jbar] ( for O equation    - B0    ) ... done with all containers   ~ ( 1 + 1 - 2) [2]  (       |  A, J ) */ 
 

  /* ~ Sequence: ~ */

  /* ~  1.) calculate phi_x, phi_y, Omega_x, Omega_y ~ */ 
  /* ~  2.) get maxb from phi (Omega?)               ~ */
  /* ~  3.) calculate [phi, Omega]                   ~ */
  /* ~  4.) place result of 3 in B0                  ~ */
  /* ~  5.) get Z_x, Z_y                             ~ */
  /* ~  6.) get maxb from Z                          ~ */
  /* ~  7.) calculate [phi, Z]                       ~ */
  /* ~  8.) place result of  7 in B2                 ~ */
  /* ~  9.) get phibar_x, phibar_y, H_x, H_y         ~ */
  /* ~ 10.) get maxb from H (A?)                     ~ */
  /* ~ 11.) calculate [phibar,H]                     ~ */
  /* ~ 12.) place result of 11 in B1                 ~ */
  /* ~ 13.) get V_x, V_y                             ~ */
  /* ~ 14.) get maxb from V                          ~ */
  /* ~ 15.) calculate [phibar, V]                    ~ */
  /* ~ 16.) place result of 15 in B3                 ~ */
  /* ~ 17.) get Zbar_x, Zbar_y, A_x, A_y             ~ */
  /* ~ 18.) calculate [A, Zbar]                      ~ */
  /* ~ 19.) add result of 18 to B1 and B3            ~ */
  /* ~ 20.) get Abar_x, Abar_y, Vbar_x, Vbar_y       ~ */
  /* ~ 21.) calculate [Abar, Vbar]                   ~ */
  /* ~ 22.) add result of 21 to B2                   ~ */
  /* ~ 23.) get Jbar_x, Jbar_y                       ~ */
  /* ~ 24.) calculate [Abar, Jbar]                   ~ */
  /* ~ 25.) add result of 24 to B0 & B2              ~ */

  /* ~ fxy sequence ~ */

  /* ~ a. ) for field f - multiply f by i            ~ */
  /* ~ b1.) for f_x multiply result of a by kx       ~ */
  /* ~ b2.) then inverse transform result            ~ */
  /* ~ c1.) for f_y multiply result a by ky          ~ */
  /* ~ c2.) then inverse transform result            ~ */

/* ~ retain for testing                              ~ */

// B0.assign(n1n2c, czero);
// B1.assign(n1n2c, czero);
// B2.assign(n1n2c, czero);
// B3.assign(n1n2c, czero);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::setD( std::string str_step, stack& run, redhallmhd& physics ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double dz;
  run.stack_data.fetch("dz"   , &dz    );         /* ~ inter-layer width along z                        ~ */
  double rho;
  physics.physics_data.fetch(  "rho"  , &rho    );         /* ~ gyro-radius parameter                            ~ */
  double beta;                           
  run.palette.fetch(   "beta" , &beta  );         /* ~ zero'th-order plasma-beta                        ~ */

  double dzm1      = one / dz;
  double rtbeta    = sqrt(beta);

  ComplexArray& Z  = run.U2;
  ComplexArray& V  = run.U3;

  ComplexArray& tZ = run.tU2;
  ComplexArray& tV = run.tU3;
  
  ComplexArray& P  = physics.P;
//  ComplexArray& A  = physics.A;
  ComplexArray& J  = physics.J;

  complex<double> deltaZ, deltaV, deltaJ, deltaP; /* ~ to aid in differentiating between the           ~ */
                                                  /* ~ predictor and corrector cases                   ~ */

  unsigned kdxp1,  kdxm1;                         /* ~ neighbor - layer indices                        ~ */
  unsigned kstart, kstop;                         /* ~ limits on k looop                               ~ */

  int n1n2c; 
  run.stack_data.fetch( "n1n2c", &n1n2c );        /* ~ number of complex elements per layer            ~ */
  int iu2;
  run.stack_data.fetch( "iu2"  , &iu2   );        /* ~ number of layers                                ~ */

  kstart    = n1n2c;                              /* ~ D's are calculated for layers 1,2,3,..n3        ~ */
  kstop     = n1n2c * (iu2 - 1);                  /* ~ layer 1 needs layer 0 and layer n3 needs layer  ~ */
                                                  /* ~ iu2 - 1                                         ~ */

  std::string model;
  run.palette.fetch("model", &model);

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    kdxm1      = kdx - n1n2c;                       /* ~ adjacent lower layer index                      ~ */
    kdxp1      = kdx + n1n2c;                       /* ~ adjacent upper layer index                      ~ */

    deltaJ     = J[ kdxp1 ] - J[ kdx   ];           /* ~ P's and J's are updated every half-step         ~ */
    deltaP     = P[ kdx   ] - P[ kdxm1 ];           /* ~ see updatePAJ. Note use of kdxp1 & kdxm1.       ~ */
    
    if ( model.compare("hall") == 0) {
      if (     str_step.compare("predict") == 0) {
      
        deltaZ = Z[ kdx   ] - Z[ kdxm1 ];           /* ~ Z and V must retain un-updated values           ~ */
        deltaV = V[ kdxp1 ] - V[ kdx   ];           /* ~ until corrector step. Note use of kdxm1 & kdxp1 ~ */

      }
      else if (str_step.compare("correct") == 0) {

        deltaZ = tZ[ kdx   ] - tZ[ kdxm1 ];         /* ~ using results of predictor step here            ~ */
        deltaV = tV[ kdxp1 ] - tV[ kdx   ];         /* ~ note use of kdxm1 & kdxp1                       ~ */

      }
    }
    else if( model.compare("rmhd") == 0 ) {

        deltaZ = czero;
        deltaV = czero;

    }

    D0[kdx]    =  deltaJ                                          * dzm1;  /* ~ i.e. Delta J / Delta z                                        ~ */
    D1[kdx]    = (          deltaP  - (       rho    *  deltaZ )) * dzm1;  /* ~ i.e. Delta F / Delta z, where F = phi - rhobar * Z            ~ */

    if ( model.compare("hall") == 0) {

      D2[kdx]  = ((rtbeta * deltaV) - (two  * rho    *  deltaJ )) * dzm1;  /* ~ i.e. Delta G / Delta z, where G = rtbeta * V - 2 * rhobar * J ~ */
      D3[kdx]  =                      (half * rtbeta * deltaZ  )  * dzm1;  /* ~ i.e. 1/2 * sqrt{beta} * Delta Z delta z                       ~ */

    }

/* ~ retained for testing ~ */

//  D0[kdx] =                 (J[kdxp1] - J[kdx  ] ) * dzm1;
//  D1[kdx] = (               (P[kdx  ] - P[kdxm1] ) -       rho * ( Z[kdx  ] - Z[kdxm1] ) ) * dzm1;  /* ~ F = phi - rhobar * Z  ~*/
//  D2[kdx] = (      rtbeta * (V[kdxp1] - V[kdx  ] ) - two * rho * ( J[kdxp1] - J[kdx  ] ) ) * dzm1;  /* ~ G = rtbeta * V - 2rhobar J ~ */
//  D3[kdx] = half * rtbeta * (Z[kdx  ] - Z[kdxm1] ) * dzm1;

  }
  
/* ~ retained for testing ~ */

//  D0.assign(nc, czero);
//  D1.assign(nc, czero);
//  D2.assign(nc, czero);
//  D3.assign(nc, czero);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::setAi( stack& run, redhallmhd& physics ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2c;
  run.stack_data.fetch( "n1n2c", &n1n2c );         /* ~ number of complex elements per layer                 ~ */
  int iu2;
  run.stack_data.fetch( "iu2",   &iu2 );           /* ~ number of layers                                     ~ */

  ComplexArray& J = physics.J;

  std::string model;
  run.palette.fetch("model", &model);

  A0.assign((n1n2c * iu2), czero);
  A1.assign((n1n2c * iu2), czero);                 /* ~ to be set to  -eta * ssqd * k2 * A = -eta * ssqd * J ~ */

  if (model.compare("hall") == 0 ) {

    A2.assign((n1n2c * iu2), czero);
    A3.assign((n1n2c * iu2), czero);

    double eta; 
    run.palette.fetch(  "eta", &eta);
    double ssqd;
    physics.physics_data.fetch("ssqd", &ssqd);
  
    unsigned kstart  = 0;
    unsigned kstop   = n1n2c * iu2;
  
    for (unsigned k  = kstart; k < kstop; k++) {
  
      A1[k]          = -( eta * ssqd * J[k] );
  
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
  
    run.tU0.reserve(n1n2c * iu2);                  /* ~ for predictor-step results           ~ */
    run.tU1.reserve(n1n2c * iu2);                  /* ~ Note: U0, U1, U2, & U3 are defined   ~ */
                                                   /* ~       on the stack.                  ~ */
  
    if (model.compare("hall") == 0 ) {
      run.tU2.reserve(n1n2c * iu2);
      run.tU3.reserve(n1n2c * iu2);
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

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::partialsInXandY( stack& run, redhallmhd& physics, ComplexArray& U, RealArray& Ux, RealArray& Uy) {

  unsigned usize  = U.capacity();

  RealArray& kx   = run.kx;
  RealArray& ky   = run.ky;
  
  unsigned kx_size = kx.capacity();

  ComplexArray U_tmp(usize, czero);

  int      n1n2c;
  unsigned idk;

  run.stack_data.fetch("n1n2c", &n1n2c   );              /* ~ number of complex elements in a layer       ~ */

  for (unsigned k = 0; k < usize; k++) { 

    if ( k % n1n2c == 0 ) { idk = 0; }
    U_tmp[k] =  iunit * kx[idk] * U[k];
    ++idk;

  } /* ~ dU/dx -> ik_x U         ~ */

  physics.fftw.fftwReverseRaw( run, U_tmp,  Ux);                                          /* ~ transform to real space ~ */
  for (unsigned k = 0; k < usize; k++) { U_tmp[k] =  iunit * ky[k] * U[k]; }             /* ~ dU/dy -> ik_y U ~ */
  physics.fftw.fftwReverseRaw( run, U_tmp,  Uy);                                          /* ~ transform to real space ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::bracket( stack& run, redhallmhd& physics, ComplexArray& BrKt, RealArray& d1x, RealArray& d1y, RealArray& d2x, RealArray& d2y) {

  unsigned dsize = d1x.capacity();

  RealArray B_tmp(dsize, zero);

  for ( unsigned k = 0; k < dsize; k++ ) { B_tmp[k] = (d1y[k] * d2x[k]) - (d1x[k] * d2y[k]); }

  physics.fftw.fftwForwardRaw( run, B_tmp, BrKt);

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

void lcsolve::averageAcrossLayers( stack& run, int shift_sign, RealArray& dx, RealArray& dy) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2;
  run.stack_data.fetch("n1n2",  &n1n2 );       /* ~ number of real elements in a layer       ~ */
  int iu2;
  run.stack_data.fetch("iu2"  , &iu2  );       /* ~ number of layers in stack                ~ */

  int dsize = dx.capacity();

  assert(dsize == n1n2 * iu2);

  RealArray d_tmp_x(dsize, zero);
  RealArray d_tmp_y(dsize, zero);

  unsigned kstart    = n1n2;
  unsigned kstop     = n1n2 * (iu2 - 1);
  unsigned kshift    = shift_sign * n1n2;

  unsigned idx       = 0;

    for (unsigned k = kstart; k < kstop; k++) {

          idx        =  k + kshift;

        d_tmp_x[k]   = half * (dx[k] + dx[idx]) ;
        d_tmp_y[k]   = half * (dy[k] + dy[idx]) ;

    }

    for (unsigned k = kstart; k < kstop; k++) {

      dx[k] = d_tmp_x[k];
      dy[k] = d_tmp_y[k];

    }

}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::Step( std::string str_step, stack& run ) {

  std::string model;
  run.palette.fetch(   "model", &model );          /* ~ either "hall" or "rmhd"              ~ */

  int n1n2c;
  run.stack_data.fetch("n1n2c" ,&n1n2c );          /* ~ number of complex elements per layer ~ */
  int iu2;
  run.stack_data.fetch("iu2"   ,&iu2   );          /* ~ number layers                        ~ */

  double dt; 
  run.stack_data.fetch("dt"    ,&dt    );          /* ~ current time increment               ~ */
  double pfrac; 
  run.stack_data.fetch("pfrac" ,&pfrac );          /* ~ fraction of dt to use for predictor- ~ */
                                                   /* ~ step                                 ~ */
  if (str_step.compare("predict") == 0 ) {         /* ~ use partial step in predictor case   ~ */
    dt             = pfrac * dt;                   /* ~ dt is local so no problem here       ~ */
  }

  int kstart       = n1n2c;                        /* ~ stepping is only done for layers 1   ~ */
  int kstop        = n1n2c * (iu2 - 1);            /* ~ through iu2 - 2                      ~ */
                                                   /* ~ layers 0 and iu2 - 1 are for the     ~ */
                                                   /* ~ boundaries and overlaps              ~ */
  int idx;                                         /* ~ an index for the S's which are       ~ */
                                                   /* ~ field-independent and thus the same  ~ */
                                                   /* ~ across layers                        ~ */

  ComplexArray& U0    = run.U0;                    /* ~ for predictor case                          ~ */
  ComplexArray& U1    = run.U1;                    /* ~ un-updated values are transferred           ~ */
  ComplexArray& U2    = run.U2;
  ComplexArray& U3    = run.U3;

  ComplexArray& tU0   = run.tU0;                   /* ~ for corrector case                          ~ */
  ComplexArray& tU1   = run.tU1;                   /* ~ results from predictor step are transferred ~ */
  ComplexArray& tU2   = run.tU2;
  ComplexArray& tU3   = run.tU3;

  for (unsigned k  = kstart; k < kstop; k++) {

    if (k % kstart == 0 ) { idx = 0; }             /* ~ reset idx when starting new layer    ~ */

      if (     str_step.compare("predict") == 0) { /* ~ the predictor case                   ~ */

        tU0[k]     = (SE0[idx] * U0[k] + (dt * (B0[k] + D0[k] + A0[k]))) * SI0[idx];
        tU1[k]     = (SE1[idx] * U1[k] + (dt * (B1[k] + D1[k] + A1[k]))) * SI1[idx];

        if ( model.compare("hall") == 0 ) {

          tU2[k]   = (SE2[idx] * U2[k] + (dt * (B2[k] + D2[k] + A2[k]))) * SI2[idx];
          tU3[k]   = (SE3[idx] * U3[k] + (dt * (B3[k] + D3[k] + A3[k]))) * SI3[idx];

        }
      }
      else if (str_step.compare("correct") == 0) { /* ~ the corrector case                   ~ */

         U0[k]     = (SE0[idx] * U0[k] + (dt * (B0[k] + D0[k] + A0[k]))) * SI0[idx];
         U1[k]     = (SE1[idx] * U1[k] + (dt * (B1[k] + D1[k] + A1[k]))) * SI1[idx];

        if ( model.compare("hall") == 0 ) {

           U2[k]   = (SE2[idx] * U2[k] + (dt * (B2[k] + D2[k] + A2[k]))) * SI2[idx];
           U3[k]   = (SE3[idx] * U3[k] + (dt * (B3[k] + D3[k] + A3[k]))) * SI3[idx];

        }
      }

    ++idx;

  }

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

   /* ~ Destructor  ~ */

  lcsolve::~lcsolve( ) {

     destroyFields();

  }
