
/* class parameter_map (definition)
 *
 * Timothy J. Dennis
 * tdennis@gi.alaska.edu
 * copyright 2014
 *
 * For the initialization and maintainance of a set or "list" of elements
 * of type parameter. Parameter type is defined in the file 'cls_parameter.h'
 *
 */

#ifndef CLS_PARAMETER_MAP
#define CLS_PARAMETER_MAP

/* #include <algorithm> 
 * #include <vector>
 */

#include "cls_parameter.hpp"
#include<map>
#include "mpi.h"

#include <fstream>
#include <iostream>
#include<iomanip>
#include<stdlib.h>
#include<sstream>

class parameter_map
{
  private:

  std::map<std::string,parameter> parameters;

  public:

  /* ~ Constructors: ~ */

  parameter_map();
  parameter_map(std::string parameter_file);

  /* ~ emplacing     ~ */

/*  void emplace(parameter::parameter pmtr); */

  void emplace(parameter pmtr);

  void emplace(std::string par_name, std::string par_val, std::string par_adj);
  void emplace(std::string par_name, int         par_val, std::string par_adj);
  void emplace(std::string par_name, float       par_val, std::string par_adj);
  void emplace(std::string par_name, double      par_val, std::string par_adj);
  void emplace(std::string par_name, bool        par_val, std::string par_adj);

  /* ~ resetting     ~ */

  bool reset(std::string str_key, std::string str_val);
  bool reset(std::string str_key, int         int_val);
  bool reset(std::string str_key, float       flt_val);
  bool reset(std::string str_key, double      dbl_val);
  bool reset(std::string str_key, bool        log_val);

  /* ~ fetching      ~ */

  bool fetch(std::string par_name, std::string *val);
  bool fetch(std::string par_name, int         *val);
  bool fetch(std::string par_name, float       *val);
  bool fetch(std::string par_name, double      *val);
  bool fetch(std::string par_name, bool        *val);

  /* ~ reporting     ~ */

  void report(std::string out_file                 );
  void report(std::string out_file_prefix, int srun);

  /* ~ Destructor (not implemented) ~ */

  ~parameter_map();

};

#endif
