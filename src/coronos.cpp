/* program coronos
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 */

#include "mpi.h"
#include "cls_stack.hpp"
#include "cls_redhallmhd.hpp"

// #include <iostream>
// #include "nsp_constants.hpp"
// #include "cls_lcsolve.hpp"
// #include "cls_parameter_map.hpp"

int main(void) {

  MPI::Init();

  stack       run(     "coronos.in");

  redhallmhd  physics( run         );
  lcsolve     solve(   run         );

  solve.Loop (         run         );

//  physics.finalize(    run, solve );

//    pm = run.palette;
//    pm.report("palette");

//    parameter_map &pm = run.run_data;

//    pm.report("run_data");

    MPI::Finalize();

}
