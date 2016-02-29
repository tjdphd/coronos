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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 *
 *        FILE: implementation of class "run_instance"
 *
 * DESCRIPTION: To be added prior to public-release
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "cls_canvas.hpp"

/* ~ Constructors ~ */

canvas::canvas() {

  /* ~ default ~ */

  std::cout << "creating empty canvas" << std::endl;

}

canvas::canvas(std::string coronos_in) {

  parameter_map cv_map(coronos_in);

  palette = cv_map;

}

/* ~ Destructor ~ */

canvas::~canvas() {

}
