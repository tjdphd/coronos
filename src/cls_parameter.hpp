/* class parameter (definition)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * For the initialization and management of the values of run parameters
 *
 */

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

  bool resetValue(std::string str_val);
  bool resetValue(int         int_val);
  bool resetValue(float       flt_val);
  bool resetValue(double      dbl_val);
  bool resetValue(bool        log_val);

  public:

  /* ~ Constructors ~ */

  parameter();

  parameter(std::string par_name, std::string par_val, std::string par_adj);
  parameter(std::string par_name, int         par_val, std::string par_adj);
  parameter(std::string par_name, float       par_val, std::string par_adj);
  parameter(std::string par_name, double      par_val, std::string par_adj);
  parameter(std::string par_name, bool        par_val, std::string par_adj);

  void reAssign(std::string par_name, std::string par_val, std::string par_adj);
  void reAssign(std::string par_name, int         par_val, std::string par_adj);
  void reAssign(std::string par_name, float       par_val, std::string par_adj);
  void reAssign(std::string par_name, double      par_val, std::string par_adj);
  void reAssign(std::string par_name, bool        par_val, std::string par_adj);


  /* ~ Destructor (not implemented) ~ */

  ~parameter();

};

#endif
