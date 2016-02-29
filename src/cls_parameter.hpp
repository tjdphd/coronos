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
 *        FILE: Definition of class "parameter"
 *
 * DESCRIPTION: For the initialization and management of the values of run 
 *              parameters
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef CLS_PARAMETER
#define CLS_PARAMETER

#include <string>
#include <sstream>

class parameter
{

  friend class parameter_map;

  private:

  std::string  name;
  std::string  type;
  std::string  adjustability;
  std::string  value;

  /* ~ Element Initializers ~ */

  void setName(std::string par_name);
  void setType(std::string par_type);
  void setAdjustability(std::string par_adj);
  void setValue(std::string par_val);

  /* ~ Resetting  Functions  ~ */

  bool resetValue(std::string str_val );
  bool resetValue(int         int_val );
  bool resetValue(float       flt_val );
  bool resetValue(double      dbl_val );
  bool resetValue(long double ldbl_val);
  bool resetValue(bool        log_val );

  public:

  /* ~ Constructors ~ */

  parameter();

       parameter(std::string par_name, std::string par_val, std::string par_adj);
       parameter(std::string par_name, int         par_val, std::string par_adj);
       parameter(std::string par_name, float       par_val, std::string par_adj);
       parameter(std::string par_name, double      par_val, std::string par_adj);
       parameter(std::string par_name, long double par_val, std::string par_adj);
       parameter(std::string par_name, bool        par_val, std::string par_adj);

  void reAssign( std::string par_name, std::string par_val, std::string par_adj);
  void reAssign( std::string par_name, int         par_val, std::string par_adj);
  void reAssign( std::string par_name, float       par_val, std::string par_adj);
  void reAssign( std::string par_name, double      par_val, std::string par_adj);
  void reAssign( std::string par_name, long double par_val, std::string par_adj);
  void reAssign( std::string par_name, bool        par_val, std::string par_adj);

  /* ~ Destructor (not implemented) ~ */

  ~parameter();

};

#endif
