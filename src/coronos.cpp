#include <iostream>
#include "nsp_constants.hpp"
#include "cls_stack.hpp"
#include "cls_lcsolve.hpp"
#include "cls_redhallmhd.hpp"
#include "cls_parameter_map.hpp"

int main(void) {

  MPI::Init();

  stack       run(     "coronos.in");
  redhallmhd  physics( run );
//  lcsolve     solve(   run );
 
//  physics.Loop (       run, solve );
//  physics.finalize(    run, solve );

//    pm = run.palette;
//    pm.report("palette");

//    parameter_map &pm = run.run_data;

//    pm.report("run_data");

    MPI::Finalize();

}
