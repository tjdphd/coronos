/* class run_instance (implementation)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * For ..
 *
 */

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
