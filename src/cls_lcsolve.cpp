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
 *        FILE: Implementation of class "lcsolve"
 *
 * DESCRIPTION: solver class for longcope solver designed to  work on
 *              canvas class type "stack"
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "cls_lcsolve.hpp"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~ */
/* ~ Constructors ~ */
/* ~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

lcsolve::lcsolve( stack& run ) {

#ifndef HAVE_CUDA_H
    createFields( run );
#endif

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef HAVE_CUDA_H

void lcsolve::Loop( stack& run ) {

  redhallmhd physics ( run );

  int     rank  ; MPI_Comm_rank(MPI_COMM_WORLD   , &rank   );
  int     n1n2  ; run.stack_data.fetch( "n1n2"   , &n1n2   );
  int     n1n2c ; run.stack_data.fetch( "n1n2c"  ,  &n1n2c );
  int     ndt   ; run.palette.fetch(    "ndt"    , &ndt    );
  int     iptest; run.palette.fetch(    "iptest" , &iptest );
  int     nw    ; run.palette.fetch(    "nw"     , &nw     );
  int     srun  ; run.palette.fetch(    "srun"   , &srun   );

  RealVar t_cur ; run.palette.fetch(    "tstart" , &t_cur  );
  RealVar dt    ;

  for (unsigned l = 0; l < ndt;l++) {

  /* ~ iptest conditional goes here              ~ */

    passAdjacentLayers( "predict", run            );
    physics.applyBC(    "predict", run            );
    physics.updatePAJ(  "predict", run            );   /* ~ P, A, and J contain un-updated/corrector-updated values ~ */

    physics.trackEnergies(     t_cur, run         );
    physics.trackQtyVsZ(       t_cur, run         );
    physics.trackPowerSpectra( t_cur, run         );
    if ( l % nw == 0 ) { run.reportEnergyQs( t_cur ); }
  
    setS(               "predict", run, physics   );   /* ~ set predictor S's                                       ~ */
    setB(               "predict", run, physics   );   /* ~ set predictor Brackets                                  ~ */
    setD(               "predict", run, physics   );   /* ~ set predictor finite differences                        ~ */
    setAi(                         run, physics   );   /* ~ set predictor A's                                       ~ */

    Step(               "predict", run            );   /* ~ execute predictor update                                ~ */

    passAdjacentLayers( "correct", run            );
    physics.updatePAJ(  "correct", run            );   /* ~ P, A, and J now contain predictor-updated values        ~ */

    setS(               "correct", run, physics   );   /* ~ set corrector S's                                       ~ */
    setB(               "correct", run, physics   );   /* ~ set corrector Brackets                                  ~ */
    setD(               "correct", run, physics   );   /* ~ set corrector finite differences                        ~ */
    setAi(                         run, physics   );   /* ~ set corrector A's                                       ~ */
    Step(               "correct", run            );   /* ~ execute corrector update                                ~ */
  
    run.palette.fetch(         "dt",   &dt        );
    t_cur = t_cur + dt;
    physics.physics_data.reset("t_cur", t_cur     );
  
    physics.updateTimeInc(         run            );

  }

  physics.applyBC(   "finalize",   run            );
  physics.updatePAJ(  "predict",   run            );   /* ~ P, A, and J contain final corrector-updated values      ~ */

  physics.fftw.fftwReverseAll(     run, physics.J );

  physics.PfromO (                 run            );   /* ~ O still in U0. Replacing with P for Primary data output ~ */
  physics.AfromH (                 run            );   /* ~ H still in U1. Replacing with A for Primary data output ~ */

  physics.physicsFinalize(         run            );

  physics.fftw.fftwReverseAll(     run            );

  physics.reportQtyVsZ( t_cur,     run            );
  physics.reportPowerSpectra(      run            );

  run.palette.reset(   "tstart", t_cur            );

  ++srun;
  run.palette.reset(   "srun"  , srun             );

  run.writeUData(                                 );
  run.writeParameters(                            );

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::passAdjacentLayers( std::string str_step, stack& run ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status * status = 0;

  ComplexArray& O     = run.U0;                           /* ~ for predictor case                          ~ */
  ComplexArray& H     = run.U1;                           /* ~ un-updated values are transferred           ~ */
  ComplexArray& Z     = run.U2;
  ComplexArray& V     = run.U3;

  ComplexArray& tO    = run.tU0;                          /* ~ for corrector case                          ~ */
  ComplexArray& tH    = run.tU1;                          /* ~ results from predictor step are transferred ~ */
  ComplexArray& tZ    = run.tU2;
  ComplexArray& tV    = run.tU3;

  int np;       run.palette.fetch(   "np",    &np      ); /* ~ number of processes                         ~ */
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   ); /* ~ number of complex elements in a layer       ~ */
  int n_layers; run.stack_data.fetch("n3"   , &n_layers); /* ~ number of layers in stack                   ~ */

  unsigned n3_idx     =   n_layers       * n1n2c;         /* ~ starting index for n3'th layer              ~ */
  unsigned atop_idx   = ( n_layers + 1 ) * n1n2c;         /* ~ starting index for top boundary layer       ~ */
  unsigned abot_idx   =                    n1n2c;         /* ~ starting index for first layer              ~ */

  std::string model; run.palette.fetch("model", &model);

  if (str_step.compare("predict") == 0) {                 /* ~ predictor case                              ~ */

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

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string model;
  run.palette.fetch("model", &model);

  RealVar  s0,  s1,  s2,   s3;                    /* ~ implicit fraction s - parameters                 ~ */
  run.palette.fetch(  "s0",   &s0 );
  run.palette.fetch(  "s1",   &s1 );

  if (model.compare("hall") == 0) {

    run.palette.fetch("s2",   &s2 );
    run.palette.fetch("s3",   &s3 );

  }

  RealVar qs0, qs1, qs2,  qs3;                    /* ~ dissipation 'q' - parameters see documentation   ~ */

  physics.physics_data.fetch("qs0", &qs0);
  physics.physics_data.fetch("qs1", &qs1);

  if (model.compare("hall") == 0) {

    physics.physics_data.fetch("qs2", &qs2);
    physics.physics_data.fetch("qs3", &qs3);

  }

  RealVar ge0, ge1, ge2,  ge3;
  RealVar gi0, gi1, gi2,  gi3;

  RealVar  dt;                                    /* ~ the current time increment                       ~ */
  run.palette.fetch("dt",   &dt );

  RealVar  pfrac;                                 /* ~ fraction of dt to use in predictor step          ~ */
  run.palette.fetch("pfrac", &pfrac);

  if ( str_step.compare("predict") == 0 ) {       /* ~ use partial time-step for predictor case         ~ */
    dt           = dt * pfrac;                    /* ~ dt is local to setS so no harm done here         ~ */
  }

  ge0            = (one - s0) * qs0 * dt;         /* ~ th g^(ex)'s  see documentation                   ~ */
  ge1            = (one - s1) * qs1 * dt;

  if (model.compare("hall") == 0) {

    ge2          = (one - s2) * qs2 * dt;
    ge3          = (one - s3) * qs3 * dt;

  }

  gi0            =        s0  * qs0 * dt;         /* ~ th g^(im)'s  see documentation                   ~ */
  gi1            =        s1  * qs1 * dt;

  if (model.compare("hall") == 0) {

    gi2          =        s2  * qs2 * dt;
    gi3          =        s3  * qs3 * dt;

  }

  RealArray& k2  = run.k2;                        /* ~ square magnitude of k-space vectors              ~ */

  int n1n2c;                                      /* ~ number of complex elements in layer              ~ */
  run.stack_data.fetch( "n1n2c", &n1n2c );

  RealArray::size_type nc = SE0.capacity();       /* ~ S's should already by sized. This is a check     ~ */
  assert (nc     == n1n2c);

  for (unsigned k = 0; k < n1n2c; k++) {          /* ~ there are only as many S-elements as k2 elements ~ */

    SE0[k]       = ( one - (ge0 * k2[k]));        /* ~ S's are initialized. See documentation           ~ */
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
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::setB( std::string str_step, stack& run, redhallmhd& physics ) {

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2;  run.stack_data.fetch( "n1n2" , &n1n2  );
  int n1n2c; run.stack_data.fetch( "n1n2c", &n1n2c );
  int iu2;   run.stack_data.fetch( "iu2"  , &iu2   );

  unsigned kstop   = n1n2c * iu2;

  std::string model;
  run.palette.fetch("model", &model);

  RealVar rho;  
  if (model.compare("hall") == 0) { physics.physics_data.fetch(  "rho", &rho ); } /* ~ gyro-radius parameter                                      ~ */
  else { rho       = zero;                                                      }
  
  RealVar beta; run.palette.fetch(          "beta", &beta );  /* ~ zero'th-order plasma-beta                                  ~ */
  RealVar rtbeta   = sqrt(beta);

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


  partialsInXandY( run, physics, P, d1x, d1y);                  /* ~ d1x, d1y hold real-space partials in x and y of P          ~ */

  maxU[0]          = maxdU(         d1x, d1y, -1, iu2);

if (str_step.compare("predict"     ) == 0) {
  partialsInXandY( run, physics, O, d2x, d2y);                  /* ~ d2x, d2y hold real-space partials in x and y of O          ~ */
}
else if (str_step.compare("correct") == 0) {
  partialsInXandY( run, physics, tO,     d2x, d2y);             /* ~ d2x, d2y hold real-space partials in x and y of tO         ~ */
}

bracket( run, physics, BrKt, d1x, d1y, d2x, d2y);               /* ~ calculate [phi, Omega]                                     ~ */

for (unsigned k = 0; k < kstop; k++) { B0[k] = - BrKt[k]; }     /* ~ place result in B0                                         ~ */

if (model.compare("hall") == 0 ) {

   if (str_step.compare("predict"     ) == 0) {
     partialsInXandY( run, physics, Z, d2x, d2y);               /* ~ d2x, d2y hold real-space partials in x and y of Z          ~ */
     maxU[2]       = maxdU(        d2x, d2y, -1, iu2);          /* ~                                                            ~ */
   }
   else if (str_step.compare("correct") == 0) {
     partialsInXandY( run, physics, tZ, d2x, d2y);              /* ~ d2x, d2y hold real-space partials in x and y of tZ         ~ */
   }
   bracket( run, physics, BrKt, d1x, d1y, d2x, d2y);            /* ~ calculate [phi, Z]                                         ~ */
   for (unsigned k = 0; k < kstop; k++) { B2[k] = - BrKt[k]; }  /* ~ place result in B2                                         ~ */
}

averageAcrossLayers( run, -1, d1x, d1y );                       /* ~ calculate ave's of phi_x & phi_y across adj't lyr's        ~ */

if (str_step.compare("predict"     ) == 0) {

  partialsInXandY( run, physics, H, d3x, d3y);                  /* ~ d3x, d3y hold real-space partials in x and y of H          ~ */

  maxU[1]          = maxdU(         d3x, d3y, 1, iu2);
}
else if (str_step.compare("correct") == 0) {
  partialsInXandY( run, physics, tH, d3x, d3y);                 /* ~ d3x, d3y hold real-space partials in x and y of tH         ~ */
}

bracket( run, physics, BrKt, d1x, d1y, d3x, d3y);               /* ~ calculate [phibar, H]                                      ~ */
for (unsigned k = 0; k < kstop; k++) { B1[k] = - BrKt[k]; }     /* ~ place result in B1                                         ~ */

if (model.compare("hall") == 0 ) {

  if (str_step.compare("predict"     ) == 0) {
    partialsInXandY( run, physics, V, d3x, d3y);                /* ~ d3x, d3y hold real-space partials in x and y of V          ~ */
    maxU[3]        = maxdU(             d3x, d3y, 1, iu2);      /* ~                                                            ~ */
  }
  else if (str_step.compare("correct") == 0) {
    partialsInXandY( run, physics, tV, d3x, d3y);               /* ~ d3x, d3y hold real-space partials in x and y of tV         ~ */
  }
  bracket( run, physics, BrKt, d1x, d1y, d3x, d3y);             /* ~ calculate [phibar, V]                                      ~ */
  for (unsigned k = 0; k < kstop; k++) { B3[k] = -BrKt[k]; }    /* ~ place result in B3                                         ~ */
  averageAcrossLayers( run,  -1, d2x, d2y );                    /* ~ calculate ave's of Z_x & Z_y across adjacent lyr's         ~ */

  partialsInXandY( run, physics, A, d1x, d1y);                  /* ~ d1x, d1y hold real-space partials in x and y of A          ~ */

//     maxu? = maxdU(               d1x, d1y);                  /* ~ might be interesting to do this calculation                ~ */

  bracket( run, physics, BrKt, d1x, d1y, d2x, d2y);             /* ~ calculate [A, Zbar]                                        ~ */
  for (unsigned k = 0; k < kstop; k++) { 
    B1[k]          = B1[k] - (rho  * BrKt[k]);
  }                                                             /* ~ B1 = - [phibar, H] - rho * [A, Zbar]                       ~ */
  for (unsigned k = 0; k < kstop; k++) { 
    B3[k] = B3[k] + (half * rtbeta * BrKt[k]); 
  }                                                             /* ~ B3 = - [phibar, V] - (1/2) * sqrt{beta} * [A, Zbar]        ~ */

}
else if(model.compare("rmhd") == 0 || model.compare("inhm") == 0) {
  partialsInXandY( run, physics, A, d1x, d1y);                  /* ~ d1x, d1y hold real-space partials in x and y of A          ~ */
}

averageAcrossLayers( run,  +1, d1x, d1y );                      /* ~ calculate ave's of A_x & A_y across adjacent lyr's         ~ */

if (model.compare("hall") == 0 ) { 

  averageAcrossLayers( run,  +1, d3x, d3y );                    /* ~ calculate ave's of V_x & V_y across adjacent layers        ~ */
  bracket( run, physics, BrKt, d1x, d1y, d3x, d3y);             /* ~ calculate [Abar, Vbar]                                     ~ */
  for (unsigned k = 0; k < kstop; k++) {
    B2[k]          = B2[k] + rtbeta * BrKt[k];                  /* ~ add result to B2                                           ~ */ 
  }
}

partialsInXandY( run, physics, J, d3x, d3y );                   /* ~ d3x, d3y hold real-space partials in x and y of J          ~ */
averageAcrossLayers( run,  +1, d3x, d3y );                      /* ~ calculate ave's of J_x & J_y across adjacent layers        ~ */

bracket( run, physics, BrKt, d1x, d1y,  d3x, d3y );             /* ~ calculate [Abar, Jbar]                                     ~ */

for (unsigned k = 0; k < kstop; k++) {                          /* ~ B0 = -[phi, Omega] + [Abar, Jbar]                          ~ */
  B0[k] = B0[k] + BrKt[k];
}

if (model.compare("hall") == 0 ) { 

  for (unsigned k = 0; k < kstop; k++) { 
    B2[k]          = B2[k] - (two * rho * BrKt[k]);             /* ~ B2 = -[phi, Z]+sqrt{beta}*[Abar, Vbar]-2*rho*[Abar, Jbar] ~ */
  }

}

int n_flds;                                                     /* ~ set rank 0 maxU                                            ~ */
run.stack_data.fetch("iu3", &n_flds);

if (str_step.compare("predict") == 0) { MPI_Allreduce(MPI_IN_PLACE, &maxU.front(), n_flds, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);}

d1y.resize(0);
d2x.resize(0);
d2y.resize(0);
d3x.resize(0);
d3y.resize(0);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::setD( std::string str_step, stack& run, redhallmhd& physics ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  RealVar dz;
  run.stack_data.fetch("dz"       , &dz    );                                  /* ~ inter-layer width along z                   ~ */
  RealVar rho;
//physics.physics_data.fetch("rho", &rho   );                                  /* ~ gyro-radius parameter                       ~ */

  std::string model;
  run.palette.fetch("model"       , &model );

  if (model.compare("hall") == 0) { physics.physics_data.fetch(  "rho", &rho ); } /* ~ gyro-radius parameter                                      ~ */
  else { rho       = zero;                                                      }

  RealVar beta;                           
  run.palette.fetch("beta"        , &beta  );                                  /* ~ zero'th-order plasma-beta                   ~ */

  RealVar dzm1     = one / dz;
  RealVar rtbeta   = sqrt(beta);

  ComplexArray& Z  = run.U2;
  ComplexArray& V  = run.U3;

  ComplexArray& tZ = run.tU2;
  ComplexArray& tV = run.tU3;
  
  ComplexArray& P  = physics.P;
  ComplexArray& J  = physics.J;

//  ComplexArray& O  = run.U0;
//  ComplexArray& A  = run.U1;

// ComplexArray& tO = run.tU0;
// ComplexArray& tA = run.tU1;

  ComplexVar deltaJ, deltaP;                                                   /* ~ to aid in differentiating between the       ~ */
//  ComplexVar deltaO, deltaA;                                                 /* ~ predictor and corrector cases               ~ */
  ComplexVar deltaZ, deltaV;

  unsigned kdxp1,  kdxm1;                                                      /* ~ neighbor - layer indices                    ~ */
  unsigned kstart, kstop;                                                      /* ~ limits on k looop                           ~ */

  int n1n2c; 
  run.stack_data.fetch( "n1n2c"   , &n1n2c );                                  /* ~ number of complex elements per layer        ~ */
  int iu2;
  run.stack_data.fetch( "iu2"     , &iu2   );                                  /* ~ number of layers                            ~ */

  int l_idx        = -1;
  kstart           = n1n2c;                                                    /* ~ D's are calculated for layers 1,2,3,..n3    ~ */
  kstop            = n1n2c * (iu2 - 1);                                        /* ~ layer 1 needs layer 0 and layer n3          ~ */
                                                                               /* ~ needs layer iu2 - 1                         ~ */

  RealArray& valfven = physics.valfven;

  for (unsigned kdx = kstart; kdx < kstop; kdx++) {

    if ( kdx % n1n2c  == 0 ) { ++l_idx; }                                      /* ~ next layer                                  ~ */

    kdxm1          = kdx - n1n2c;                                              /* ~ adjacent lower layer index                  ~ */
    kdxp1          = kdx + n1n2c;                                              /* ~ adjacent upper layer index                  ~ */

    deltaJ         = valfven[l_idx] * (J[ kdxp1 ] - J[ kdx   ]);               /* ~ P's and J's are updated every half-step     ~ */
    deltaP         = valfven[l_idx] * (P[ kdx   ] - P[ kdxm1 ]);               /* ~ see updatePAJ. Note use of kdxp1 & kdxm1.   ~ */

//   if ( model.compare("inhm") == 0) { 
//     if (     str_step.compare("predict") == 0) {

//       deltaO     = O[kdxp1]  - O[kdxm1]
//       deltaA     = A[kdxp1]  - A[kdxm1] /*~ WARNING - must deal with lowest kdxm1 ~*/

//     }
//     else if (str_step.compare("correct") == 0) {

//       deltaO     = tO[kdxp1] - tO[kdxm1]
//       deltaA     = tA[kdxp1] - tA[kdxm1]/*~ WARNING - must deal with lowest kdxm1 ~*/

//     }
//   }
    
    if ( model.compare("hall") == 0) {
      if (     str_step.compare("predict") == 0) {
      
        deltaZ     = Z[ kdx   ] - Z[ kdxm1 ];                                  /* ~ Z and V must retain un-updated values until ~ */
        deltaV     = V[ kdxp1 ] - V[ kdx   ];                                  /* ~ corrector step. Note use of kdxm1 & kdxp1   ~ */

      }
      else if (str_step.compare("correct") == 0) {

        deltaZ     = tZ[ kdx   ] - tZ[ kdxm1 ];                                /* ~ using results of predictor step here        ~ */
        deltaV     = tV[ kdxp1 ] - tV[ kdx   ];                                /* ~ note use of kdxm1 & kdxp1                   ~ */

      }
    }
    else if(model.compare("rmhd") == 0 || model.compare("inhm") == 0) {

        deltaZ     = czero;
        deltaV     = czero;

    }

    D0[kdx]        =  deltaJ                                          * dzm1;  /* ~ i.e. Delta J / Delta z                      ~ */
    D1[kdx]        = (          deltaP  - (       rho    *  deltaZ )) * dzm1;  /* ~ i.e. Delta F / Delta z,                     ~ */
                                                                               /* ~ where F = phi - rhobar * Z                  ~ */
    if ( model.compare("hall") == 0) {

      D2[kdx]      = ((rtbeta * deltaV) - (two  * rho    *  deltaJ )) * dzm1;  /* ~ i.e. Delta G / Delta z,                     ~ */
                                                                               /* ~ where G = rtbeta * V - 2 * rhobar * J       ~ */
      D3[kdx]      =                      (half * rtbeta * deltaZ  )  * dzm1;  /* ~ i.e. 1/2 * sqrt{beta} * Delta Z delta z     ~ */

    }
  }

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::setAi( stack& run, redhallmhd& physics ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2c;
  run.stack_data.fetch( "n1n2c", &n1n2c );         /* ~ number of complex elements per layer                 ~ */
  int iu2;
  run.stack_data.fetch( "iu2",   &iu2 );           /* ~ number of layers                                     ~ */

  unsigned usize    = n1n2c * iu2;

  ComplexArray& J   = physics.J;                   /* ~ required for all models                              ~ */

  ComplexArray& P   = physics.P;                   /* ~ required for inhm                                    ~ */

  ComplexArray Jbar;
  ComplexArray Pbar;

  RealArray& h11    = physics.h11;
  RealArray& h12    = physics.h12;
  RealArray& h21    = physics.h21;
  RealArray& h22    = physics.h22;

  std::string model;
  run.palette.fetch("model", &model);

  A0.assign((usize), czero);
  A1.assign((usize), czero);                 /* ~ to be set to  -eta * ssqd * k2 * A = -eta * ssqd * J ~ */

  if (model.compare("inhm") == 0 ) {
    
    Jbar.assign(usize, czero);
    Pbar.assign(usize, czero);

    for (unsigned k = 0; k < usize; k++) {Jbar[k]  = J[k];}
    for (unsigned k = 0; k < usize; k++) {Pbar[k]  = P[k];}

    averageAcrossLayers( run, -1, Pbar );
    averageAcrossLayers( run, +1, Jbar );

    int l_idx       = -1;
    unsigned kstart = n1n2c;
    unsigned kstop  = n1n2c * (iu2 - 1);

    for (unsigned kdx = kstart; kdx < kstop; kdx++) {
   
      if (kdx % n1n2c == 0 ) { ++l_idx; }

      A0[ kdx ]     = h12[l_idx] * Jbar[kdx];           /* ~ A0 - L1 = h11(z) * Omega + h12(z)*jbar                       ~ */
      A1[ kdx ]     = h21[l_idx] * Pbar[kdx];           /* ~ A1 - L2 = h21(z)*phibar  + h22(z)*A                          ~ */

    }

     Jbar.erase(Jbar.begin(),Jbar.end());
     Pbar.erase(Pbar.begin(),Pbar.end());

  }

  if (model.compare("hall") == 0 ) {

    A2.assign((n1n2c * iu2), czero);
    A3.assign((n1n2c * iu2), czero);

    RealVar eta; 
    run.palette.fetch(  "eta", &eta);
    RealVar ssqd;
    run.palette.fetch("ssqd", &ssqd);
  
    unsigned kstart = 0;
    unsigned kstop  = n1n2c * iu2;
  
    for (unsigned k = kstart; k < kstop; k++) { A1[k]         = -( eta * ssqd * J[k] ); }

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

  unsigned usize   = U.capacity();

  RealArray& kx    = run.kx;
  RealArray& ky    = run.ky;
  
  unsigned kx_size = kx.capacity();

  ComplexArray U_tmp(usize, czero);

  int      n1n2c;
  run.stack_data.fetch("n1n2c", &n1n2c   );         /* ~ number of complex elements in a layer ~ */

  unsigned idk;

  for (unsigned k  = 0; k < usize; k++) { 

    if ( k % n1n2c == 0 ) { idk = 0; }
    U_tmp[k]       =  iunit * kx[idk] * U[k];
    ++idk;

  } /* ~ dU/dx -> ik_x U         ~ */

  physics.fftw.fftwReverseRaw( run, U_tmp,  Ux);    /* ~ transform to real space               ~ */

  for (unsigned k  = 0; k < usize; k++) { 
    if ( k % n1n2c == 0 ) { idk = 0; }
    U_tmp[k]       =  iunit * ky[idk] * U[k]; 
    ++idk;
    
  } /* ~ dU/dy -> ik_y U ~ */

  physics.fftw.fftwReverseRaw( run, U_tmp,  Uy);    /* ~ transform to real space               ~ */

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::bracket( stack& run, redhallmhd& physics, ComplexArray& BrKt, RealArray& d1x, RealArray& d1y, RealArray& d2x, RealArray& d2y) {

  unsigned dsize   = d1x.capacity();

  RealArray B_tmp(dsize, zero);

  for ( unsigned k = 0; k < dsize; k++ ) { B_tmp[k] = (d1y[k] * d2x[k]) - (d1x[k] * d2y[k]); }

  physics.fftw.fftwForwardRaw( run, B_tmp, BrKt);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

RealVar lcsolve::maxdU( RealArray& dx, RealArray& dy, int i_grid, int i_layers) {

  unsigned dsize   = dx.capacity();
  unsigned k_size, k_start, k_stop;
  unsigned n3;

  k_size           = dsize / i_layers;
  n3               =  i_layers - 2;
  
  if ( i_grid == - 1)    {    k_start = 0;   }                       /* ~ lower grid ~ */
  else if ( i_grid == 1) { k_start = k_size; }                       /* ~ upper grid ~ */
  else {
         std::cout << "maxdU: ERROR - invalid grid" << std::endl; 
         k_start   =  dsize;
       }                                                             /* ~ error      ~ */

  k_stop  = k_start + (n3 + 1) * k_size;

  RealVar test     = 0;
  RealVar max      = 0;

   for (unsigned k = k_start; k < k_stop; k++ ) {

     test          = ((dx[k] * dx[k]) + (dy[k] * dy[k]));

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

  int dsize          = dx.capacity();

  assert(dsize       == n1n2 * iu2);

  RealArray d_tmp_x(dsize, zero);
  RealArray d_tmp_y(dsize, zero);

  unsigned kstart    = n1n2;
  unsigned kstop     = n1n2 * (iu2 - 1);
  unsigned kshift    = shift_sign * n1n2;

  unsigned idx       = 0;

    for (unsigned k  = kstart; k < kstop; k++) {

          idx        =  k + kshift;

        d_tmp_x[k]   = half * (dx[k] + dx[idx]) ;
        d_tmp_y[k]   = half * (dy[k] + dy[idx]) ;

    }

    for (unsigned k  = kstart; k < kstop; k++) {

      dx[k]          = d_tmp_x[k];
      dy[k]          = d_tmp_y[k];

    }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::averageAcrossLayers( stack& run, int shift_sign, ComplexArray& dx ) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n1n2c;
  run.stack_data.fetch("n1n2c",  &n1n2c );       /* ~ number of complex elements in a layer    ~ */
  int iu2;
  run.stack_data.fetch("iu2"  , &iu2    );       /* ~ number of layers in stack                ~ */

  int dsize          = dx.capacity();

  assert(dsize       == n1n2c * iu2);

  ComplexArray d_tmp_x(dsize, czero);

  unsigned kstart    = n1n2c;
  unsigned kstop     = n1n2c * (iu2 - 1);
  unsigned kshift    = shift_sign * n1n2c;

  unsigned idx       = 0;

    for (unsigned k  = kstart; k < kstop; k++) {
          idx        =  k + kshift;
        d_tmp_x[k]   = half * ( dx[k] + dx[idx] ) ;
    }

    for (unsigned k  = kstart; k < kstop; k++) {
      dx[k]          = d_tmp_x[k];
    }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void lcsolve::Step( std::string str_step, stack& run ) {

  std::string model; run.palette.fetch(   "model" , &model ); /* ~ one of "hall" or "rmhd", or "inhm"          ~ */
  int n1n2c;         run.stack_data.fetch("n1n2c" , &n1n2c ); /* ~ number of complex elements per layer        ~ */
  int iu2;           run.stack_data.fetch("iu2"   , &iu2   ); /* ~ number layers                               ~ */
  RealVar dt;        run.palette.fetch(   "dt"    , &dt    ); /* ~ current time increment                      ~ */
  RealVar pfrac;     run.palette.fetch(   "pfrac" , &pfrac ); /* ~ fraction of dt to use for predictor-        ~ */
                                                              /* ~ step                                        ~ */
  if (str_step.compare("predict") == 0 ) {                    /* ~ use partial step in predictor case          ~ */
    dt              = pfrac * dt;                             /* ~ dt is local so no problem here              ~ */
  }

  int kstart        = n1n2c;                                  /* ~ stepping is only done for layers 1          ~ */
  int kstop         = n1n2c * (iu2 - 1);                      /* ~ through iu2 - 2                             ~ */
                                                              /* ~ layers 0 and iu2 - 1 are for the            ~ */
                                                              /* ~ boundaries and overlaps                     ~ */
  int idx;                                                    /* ~ an index for the S's which are              ~ */
                                                              /* ~ field-independent and thus the same         ~ */
                                                              /* ~ across layers                               ~ */

  ComplexArray& U0  = run.U0;                                 /* ~ for predictor case                          ~ */
  ComplexArray& U1  = run.U1;                                 /* ~ un-updated values are transferred           ~ */
  ComplexArray& U2  = run.U2;
  ComplexArray& U3  = run.U3;

  ComplexArray& tU0 = run.tU0;                                /* ~ for corrector case                          ~ */
  ComplexArray& tU1 = run.tU1;                                /* ~ results from predictor step are transferred ~ */
  ComplexArray& tU2 = run.tU2;
  ComplexArray& tU3 = run.tU3;

  for (unsigned k   = kstart; k < kstop; k++) {

    if (k % kstart  == 0 ) { idx = 0; }                      /* ~ reset idx when starting new layer            ~ */

      if (     str_step.compare("predict") == 0) {           /* ~ the predictor case                           ~ */

          tU0[k]    = (SE0[idx] * U0[k] + (dt * (B0[k] + D0[k] + A0[k]))) * SI0[idx];
          tU1[k]    = (SE1[idx] * U1[k] + (dt * (B1[k] + D1[k] + A1[k]))) * SI1[idx];

        if ( model.compare("hall") == 0 ) {

          tU2[k]    = (SE2[idx] * U2[k] + (dt * (B2[k] + D2[k] + A2[k]))) * SI2[idx];
          tU3[k]    = (SE3[idx] * U3[k] + (dt * (B3[k] + D3[k] + A3[k]))) * SI3[idx];

        }
      }
      else if (str_step.compare("correct") == 0) { /* ~ the corrector case                          ~ */

         U0[k]      = (SE0[idx] * U0[k] + (dt * (B0[k] + D0[k] + A0[k]))) * SI0[idx];
         U1[k]      = (SE1[idx] * U1[k] + (dt * (B1[k] + D1[k] + A1[k]))) * SI1[idx];

        if ( model.compare("hall") == 0 ) {

           U2[k]    = (SE2[idx] * U2[k] + (dt * (B2[k] + D2[k] + A2[k]))) * SI2[idx];
           U3[k]    = (SE3[idx] * U3[k] + (dt * (B3[k] + D3[k] + A3[k]))) * SI3[idx];

        }
      }

    ++idx;

  }
}

#endif
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

   /* ~ Destructor  ~ */

  lcsolve::~lcsolve( ) {

#ifndef HAVE_CUDA_H

     destroyFields();

#endif

  }
