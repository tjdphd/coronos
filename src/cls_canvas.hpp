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
 *        FILE: definition of class "canvas"
 *
 * DESCRIPTION: child class of run_instance  contains conmplete information 
 *              about problem and solution in its member classes
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef CLS_CANVAS
#define CLS_CANVAS

#include "cls_run_instance.hpp"

#ifdef HAVE_CUDA_H
  #include "cls_canvas_cuda_ext.hpp"
#endif

class canvas : public run_instance
{
   private:


   public:

   parameter_map palette;

   /* ~ Constructor ~ */

   canvas();
   canvas(std::string coronos_in);

   /* ~ Destructor  ~ */
	
   ~canvas();
};

#endif
