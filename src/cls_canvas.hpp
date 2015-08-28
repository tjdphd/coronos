/* class canvas (definition)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * child class of run instance
 * contains conmplete information
 * about problem and solution
 * in its member classes
 *
 */

#ifndef CLS_CANVAS
#define CLS_CANVAS

#include "cls_run_instance.hpp"

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
