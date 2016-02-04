/* program coronos
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 */

#include "mpi.h"
#include "cls_stack.hpp"
#include "cls_lcsolve.hpp"

int main(void) {

  MPI::Init();

  stack       run(     "coronos.in");
  lcsolve     solve(   run         );
#ifndef HAVE_CUDA_H
  solve.Loop (         run         );
#endif

  MPI::Finalize();

}
