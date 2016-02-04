/* class run_instance (definition)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * Parent class for a run instance of coronos
 *
 */

#ifndef CLS_RUN_INSTANCE
#define CLS_RUN_INSTANCE


 #include "cls_parameter_map.hpp"
 #include <ctime>
 #include <unistd.h>

 #ifdef HAVE_CUDA_H
   #include "cls_run_instance_cuda_ext.hpp"
 #endif

class run_instance
{
   private:

   std::string getTime();
   std::string getNode();

   public:

   parameter_map run_data;

   /* ~ Constructor ~ */

   run_instance();

   /* ~ Destructor  ~ */

   ~run_instance();

};

#endif
