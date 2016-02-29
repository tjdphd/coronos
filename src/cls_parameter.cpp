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
 *        FILE: Implementation of class "parameter"
 *
 * DESCRIPTION: For the initialization and management of the values of run 
 *              parameters
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "cls_parameter.hpp"

#include <iostream>
#include <string>

/* ~ Constructors ~ */

parameter::parameter() {

  std::string empty = "empty";

           name.assign(empty);
           type.assign(empty);
  adjustability.assign(empty);
          value.assign(empty);

}

parameter::parameter(std::string par_name, std::string par_val, std::string par_adj) {

  std::string str_val;

  str_val.assign(par_val);

           name.assign(par_name);
          value.assign(str_val );
           type.assign("str"   );
  adjustability.assign(par_adj );

}

parameter::parameter(std::string par_name, int         par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
          value.assign(str_val );
           type.assign("int"   );
  adjustability.assign(par_adj );
}

parameter::parameter(std::string par_name, float       par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
          value.assign(str_val );
           type.assign("flt"   );
  adjustability.assign(par_adj );

}

parameter::parameter(std::string par_name, double      par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
           value.assign(str_val);
           type.assign("dbl"   );
  adjustability.assign(par_adj );

}

parameter::parameter(std::string par_name, long double  par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
           value.assign(str_val);
           type.assign("dbl"   );
  adjustability.assign(par_adj );

}

parameter::parameter(std::string par_name, bool        par_val, std::string par_adj) {

  std::string str_val;

  if (par_val == true) {

    str_val.assign("true");

  }
  else {

    str_val.assign("false");

  }
           name.assign(par_name  );
             value.assign(str_val);
           type.assign("log"     );
  adjustability.assign(par_adj   );

}


// ~ Initializers ~ //

// ~ type ~ //

void parameter::setName(std::string par_name) {

  name.assign(par_name);

}

void parameter::setType(std::string par_type) {

       if (par_type.compare("str") == 0 || \
           par_type.compare("int") == 0 || \
           par_type.compare("flt") == 0 || \
           par_type.compare("dbl") == 0 || \
           par_type.compare("log") == 0)
       {
           type.assign(par_type);
       }
       else
       {
         type.assign("xxx");
       }
}

// ~ adjustability ~ //

void parameter::setAdjustability(std::string par_adj) {

  if ( par_adj.compare("adj") == 0 || \
       par_adj.compare("rfx") == 0 || \
       par_adj.compare("sfx") == 0)
  {
     adjustability.assign(par_adj);
  }
  else 
  {
     // should put some real error-handling here;
     adjustability.assign("xxx");
  }
}

// ~ value ~ //

void parameter::setValue(std::string par_val) {

    value.assign(par_val);
     
}

bool parameter::resetValue(std::string str_val) {

     bool l_reset = false;
    
     if (!adjustability.compare("rfx")       == 0)   {

       if (adjustability.compare("sfx")      == 0) {

         adjustability.assign("sfa");
         l_reset  = true;

       }
       else if (adjustability.compare("adj") == 0) {

         l_reset  = true;

       }

       if (l_reset) value.assign(str_val);

     }

     return l_reset;
}

bool parameter::resetValue(int         int_val) {

     bool l_reset        = false;
     std::string str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << int_val) ) -> str();

     if (!adjustability.compare("rfx")       == 0)   {

       if (adjustability.compare("sfx")      == 0) {

         adjustability.assign("sfa");
         l_reset  = true;
       }
       else if (adjustability.compare("adj") == 0) {

         l_reset  = true;

       }

       if(l_reset) value.assign(str_val);

     }

     return l_reset;
}

bool parameter::resetValue(float       flt_val) {

     bool l_reset        = false;
     std::string str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << flt_val) ) -> str();

     if (!adjustability.compare("rfx") == 0) {

       if (adjustability.compare("sfx") == 0) {

         adjustability.assign("sfa");
         l_reset  = true;

       }
       else if (adjustability.compare("adj") == 0) {

         l_reset  = true;

       }

       if (l_reset) value.assign(str_val);

     }

     return l_reset;
}

bool parameter::resetValue(double      dbl_val) {

     bool l_reset = false;
     std::string str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << dbl_val) ) -> str();

     if (!adjustability.compare("rfx")       == 0)   {

       if (adjustability.compare("sfx")      == 0) {

         adjustability.assign("sfa");
         l_reset  = true;
       }
       else if (adjustability.compare("adj") == 0) {

         l_reset  = true;

       }

       if (l_reset) value.assign(str_val);

     }

     return l_reset;
}

bool parameter::resetValue(long double  dbl_val) {

     bool l_reset = false;
     std::string str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << dbl_val) ) -> str();

     if (!adjustability.compare("rfx")       == 0)   {

       if (adjustability.compare("sfx")      == 0) {

         adjustability.assign("sfa");
         l_reset  = true;
       }
       else if (adjustability.compare("adj") == 0) {

         l_reset  = true;

       }

       if (l_reset) value.assign(str_val);

     }

     return l_reset;
}

bool parameter::resetValue(bool        log_val) {

     bool l_reset = false;
     std::string str_val;

     if (log_val) {
       str_val.assign("true");
     }
     else {
       str_val.assign("false");
     }

     if (!adjustability.compare("rfx")       == 0)   {

       if (adjustability.compare("sfx")      == 0) {

         adjustability.assign("sfa");
         l_reset  = true;
       }
       else if (adjustability.compare("adj") == 0) {

         l_reset  = true;

       }

       if (l_reset) value.assign(str_val);

     }

     return l_reset;
}

void parameter::reAssign(std::string par_name, std::string par_val, std::string par_adj) {

  std::string str_val;

  str_val.assign(par_val);

           name.assign(par_name);
          value.assign(str_val );
           type.assign("str"   );
  adjustability.assign(par_adj );

}

void parameter::reAssign(std::string par_name, int         par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
          value.assign(str_val );
           type.assign("int"   );
  adjustability.assign(par_adj );

}

void parameter::reAssign(std::string par_name, float       par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
          value.assign(str_val );
           type.assign("flt"   );
  adjustability.assign(par_adj );

}

void parameter::reAssign(std::string par_name, double      par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
           value.assign(str_val);
           type.assign("dbl"   );
  adjustability.assign(par_adj );

}

void parameter::reAssign(std::string par_name, long double par_val, std::string par_adj) {

  std::string str_val;

  str_val = static_cast<std::ostringstream*>( &(std::ostringstream() << par_val) ) -> str();

           name.assign(par_name);
           value.assign(str_val);
           type.assign("dbl"   );
  adjustability.assign(par_adj );

}

void parameter::reAssign(std::string par_name, bool        par_val, std::string par_adj) {

  std::string str_val;

  if (par_val == true) {

    str_val.assign("true");

  }
  else {

    str_val.assign("false");

  }
           name.assign(par_name  );
             value.assign(str_val);
           type.assign("log"     );
  adjustability.assign(par_adj   );

}


/* ~ Destructor ~ */

parameter::~parameter() {

  /* ~ not implemented ~ */
}

