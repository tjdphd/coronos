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
 *        FILE: Implementation of class "parameter_map"
 *
 * DESCRIPTION: For the initialization and maintainance of a set of elements 
 *              of type parameter. Parameter type is defined in the file 
 *              'cls_parameter.h'
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "cls_parameter_map.hpp"

// ~ constructors ~ //

parameter_map::parameter_map()
{
   parameters.clear();
}

parameter_map::parameter_map(std::string parameter_file)
{

  std::ifstream infile;                 //  file stream pointer for parameter_file
  std::string   file_line;              //  holds a line from parameter_file
  std::string   field;                  //  holds a space-separated field from line

  const char *to_be_read;               //  infile requires a const char
  to_be_read        = parameter_file.c_str();

  infile.open(to_be_read, std::ios::in);

  size_t non_blank;                     //  for locating field
  size_t blank;                         //  for locating blanks between fields
  size_t field_length;
  size_t str_length;                    //  location of line end

  int field_count;                      //  tracks number of fields encountered in file line
  int line_count       = 0;              //  counts number of file-lines read

  parameter   pmtr;
  std::string pmtr_key;

  while(getline(infile, file_line))     // Get next line from infile
  {
    ++line_count;
    field_count        = 0;
    str_length         = file_line.length();
    non_blank          = file_line.find_first_not_of(" ");            // locate first non-blank character in line
    blank              = file_line.find_first_of(" ", non_blank + 1); // locate first blank after beginning of first

    while(blank <= str_length)
    {
      field_length     = blank - non_blank;                           // Number of characters in first field
      field            = file_line.substr(non_blank, field_length);   // The first field in the line
      ++field_count;

      switch(field_count)
      {
        case 1 : pmtr.setName(field);
                 pmtr_key.assign(field);
                 break;
        case 2 : pmtr.setValue(field);
                 break;
        case 3 : pmtr.setType(field);
                 break;
        case 4 : pmtr.setAdjustability(field);
                 break;
      }

      non_blank        = file_line.find_first_not_of(" ", blank + 1);  // location of the next field's first character

      if (non_blank   == file_line.npos) non_blank = str_length;       // if non_blank is location of 
                                                                     // last character of last field
         
          if (blank   != str_length && non_blank != str_length)        // we're not at the end of the line
          {
            blank      = file_line.find_first_of(" ", non_blank + 1);// where the next field ends
            if (blank == file_line.npos) blank = str_length;         // in case there's a blank after the last
          }                                                          // field
          else blank=str_length + 1;                                 // We are at the end of the line
         
     }

     parameters.insert(std::pair<std::string,parameter>(pmtr_key, pmtr));
     
  }

  infile.close();

}

  // ~ emplacing     ~ //

void parameter_map::emplace(parameter pmtr) {

  std::string pmtr_key;
  pmtr_key = pmtr.name;
  parameters.insert(std::pair<std::string,parameter>(pmtr_key, pmtr));
}

void parameter_map::emplace(std::string par_name, std::string par_val, std::string par_adj) {

  parameter pmtr(par_name, par_val, par_adj);
  parameters.insert(std::pair<std::string,parameter>(par_name, pmtr));

}

void parameter_map::emplace(std::string par_name, int         par_val, std::string par_adj) {

  parameter pmtr(par_name, par_val, par_adj);
  parameters.insert(std::pair<std::string,parameter>(par_name, pmtr));
}

void parameter_map::emplace(std::string par_name, float       par_val, std::string par_adj) {

  parameter pmtr(par_name, par_val, par_adj);
  parameters.insert(std::pair<std::string,parameter>(par_name, pmtr));
}

void parameter_map::emplace(std::string par_name, double      par_val, std::string par_adj) {

  parameter pmtr(par_name, par_val, par_adj);
  parameters.insert(std::pair<std::string,parameter>(par_name, pmtr));
}

void parameter_map::emplace(std::string par_name, long double  par_val, std::string par_adj) {

  parameter pmtr(par_name, par_val, par_adj);
  parameters.insert(std::pair<std::string,parameter>(par_name, pmtr));
}

void parameter_map::emplace(std::string par_name, bool         par_val, std::string par_adj) {

  parameter pmtr(par_name, par_val, par_adj);
  parameters.insert(std::pair<std::string,parameter>(par_name, pmtr));

}


/* ~ resetting ~ */

bool parameter_map::reset(std::string str_key, std::string str_val) {

     bool          l_success         = false;
     parameter     pmtr              = parameters[str_key];

     l_success = parameters[str_key].resetValue(str_val);
   
     return l_success;
}

bool parameter_map::reset(std::string str_key, int         int_val) {

     bool          l_success         = false;
     parameter     pmtr              = parameters[str_key];

     l_success = parameters[str_key].resetValue(int_val);

     return l_success;
}

bool parameter_map::reset(std::string str_key, float       flt_val) {

     bool          l_success         = false;
     parameter     pmtr              = parameters[str_key];

     l_success = parameters[str_key].resetValue(flt_val);
   
     return l_success;
}
bool parameter_map::reset(std::string str_key, double      dbl_val) {

     bool          l_success         = false;
     parameter     pmtr              = parameters[str_key];

     l_success = parameters[str_key].resetValue(dbl_val);
   
     return l_success;
}

bool parameter_map::reset(std::string str_key, long double      dbl_val) {

     bool          l_success         = false;
     parameter     pmtr              = parameters[str_key];

     l_success = parameters[str_key].resetValue(dbl_val);
   
     return l_success;
}

bool parameter_map::reset(std::string str_key, bool        log_val) {

     bool          l_success         = false;
     parameter     pmtr              = parameters[str_key];

     l_success = parameters[str_key].resetValue(log_val);
   
     return l_success;
}

// ~ fetching ~ //

bool parameter_map::fetch(std::string par_name, std::string *val) {

     parameter   pmtr;
     std::string val_str;

     pmtr    = parameters[par_name];
     if (pmtr.name.compare("empty") == 0) 
     { std::cout << "fetch: WARNING - pmtr is empty for parameter " << par_name << std::endl;}
     val_str = pmtr.value;
        *val = val_str;

     return true;
}

bool parameter_map::fetch(std::string par_name, int         *val) {

     parameter       pmtr;
     std::string  val_str;

     pmtr    = parameters[par_name];
     if (pmtr.name.compare("empty") == 0) 
     { std::cout << "fetch: WARNING - pmtr is empty for parameter " << par_name << std::endl;}
     val_str = pmtr.value;
     *val    = atoi(val_str.c_str());

     return true;
}

bool parameter_map::fetch(std::string par_name, float       *val) {

     parameter       pmtr;
     std::string  val_str;

     pmtr    = parameters[par_name];
     if (pmtr.name.compare("empty") == 0) 
     { std::cout << "fetch: WARNING - pmtr is empty for parameter " << par_name << std::endl;}
     val_str = pmtr.value;
     *val    = atof(val_str.c_str());

     return true;

}
bool parameter_map::fetch(std::string par_name, double      *val) {

     parameter       pmtr;
     std::string  val_str;

     pmtr    = parameters[par_name];
     if (pmtr.name.compare("empty") == 0) 
     { std::cout << "fetch: WARNING - pmtr is empty for parameter " << par_name << std::endl;}
     val_str = pmtr.value;
     *val    = atof(val_str.c_str());

     return true;
}

bool parameter_map::fetch(std::string par_name, long double      *val) {

     parameter       pmtr;
     std::string  val_str;

     pmtr    = parameters[par_name];
     if (pmtr.name.compare("empty") == 0) 
     { std::cout << "fetch: WARNING - pmtr is empty for parameter " << par_name << std::endl;}
     val_str = pmtr.value;
     *val    = atof(val_str.c_str());

     return true;
}

bool parameter_map::fetch(std::string par_name, bool        *val) {

     parameter       pmtr;
     std::string  val_str;

     pmtr    = parameters[par_name];
     if (pmtr.name.compare("empty") == 0) 
     { std::cout << "fetch: WARNING - pmtr is empty for parameter " << par_name << std::endl;}
     val_str = pmtr.value;

     if (val_str.compare("true")      == 0) {
       *val  = true;
     }
     else if(val_str.compare("false") == 0) {
       *val  = false;
     }
     else {
       *val  = false;
     }

     return true;
}

/* ~   Reporting     ~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void parameter_map::report(std::string out_file_prefix, std::string out_file_run_label, int srun) {

     std::map<std::string, parameter>::iterator it;

     parameter   pmtr;


     std::string str_val;
     int         int_val;
     float       flt_val;
     long double dbl_val;
     bool        log_val;

     int         str_key_len;
     int         max_str_key_len;

     int         str_val_len;
     int         max_str_val_len;

     max_str_key_len     = 0;
     max_str_val_len     = 0;

     for (it = parameters.begin(); it != parameters.end(); ++it) {

       pmtr              = parameters[it->first];

       str_key_len       = pmtr.name.length();
       if (str_key_len > max_str_key_len) max_str_key_len = str_key_len;
       
       if (pmtr.type.compare("str") == 0) {
 
         str_val_len     = pmtr.value.length();
         if (str_val_len > max_str_val_len) max_str_val_len = str_val_len;

       }
     }
     
     if (max_str_val_len < 24) max_str_val_len = 24;

     int rank;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

     std::string    str_srun               = static_cast<std::ostringstream*>( &(std::ostringstream() << srun) ) -> str();
     std::string    str_rank               = static_cast<std::ostringstream*>( &(std::ostringstream() << rank) ) -> str();
     int rnk_len                           = str_rank.length();
     switch(rnk_len) {

     case(1) : str_rank = "0" + str_rank;
     case(2) : ;
     default : ;

     }

     std::string    out_file               = out_file_prefix + "." + str_rank + ".o" + out_file_run_label + str_srun;
     const char    *c_str_out_file         = out_file.c_str();
     std::fstream   ofs;

     ofs.open(c_str_out_file, std::fstream::out);

     for (it = parameters.begin(); it     != parameters.end(); ++it) {

        pmtr                               = parameters[it->first];

        ofs  << std::setw(max_str_key_len + 2) << std::left << pmtr.name;

        if      (pmtr.type.compare("str") ==       0)
        {
          fetch(pmtr.name, &str_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << str_val;
        }
        else if (pmtr.type.compare("int") ==       0)
        {
          fetch(pmtr.name, &int_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << int_val;
        }
        else if (pmtr.type.compare("flt") ==       0) 
        {
          fetch(pmtr.name, &flt_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << std::setprecision(15) << std::scientific;
          ofs << flt_val;
        }
        else if (pmtr.type.compare("dbl") ==       0) 
        {
          fetch(pmtr.name, &dbl_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << std::setprecision(15) << std::scientific;
          ofs <<  dbl_val;
        }
        else if (pmtr.type.compare("log") ==       0) 
        {
          fetch(pmtr.name, &log_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << log_val;
        }
        
        ofs   << std::setw(4)                   << pmtr.type;
        ofs   << std::setw(4)                   << pmtr.adjustability;
        ofs   << std::endl;

     }

     ofs.close();
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void parameter_map::report(std::string out_file) {

     std::map<std::string, parameter>::iterator it;

     parameter   pmtr;

     std::string str_val;
     int         int_val;
     float       flt_val;
     long double dbl_val;
     bool        log_val;

     int         str_key_len;
     int         max_str_key_len;

     int         str_val_len;
     int         max_str_val_len;

     max_str_key_len     = 0;
     max_str_val_len     = 0;

     for (it = parameters.begin(); it != parameters.end(); ++it) {

       pmtr              = parameters[it->first];

       str_key_len       = pmtr.name.length();
       if (str_key_len > max_str_key_len) max_str_key_len = str_key_len;
       
       if (pmtr.type.compare("str") == 0) {
 
         str_val_len     = pmtr.value.length();
         if (str_val_len > max_str_val_len) max_str_val_len = str_val_len;

       }
     }
     
     if (max_str_val_len < 24) max_str_val_len = 24;

     const char    *c_str_out_file         = out_file.c_str();
     std::fstream   ofs;

     ofs.open(c_str_out_file, std::fstream::out);

     for (it = parameters.begin(); it     != parameters.end(); ++it) {

        pmtr                               = parameters[it->first];


        ofs  << std::setw(max_str_key_len + 2) << std::left << pmtr.name;

        if      (pmtr.type.compare("str") ==       0)
        {
          fetch(pmtr.name, &str_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << str_val;
        }
        else if (pmtr.type.compare("int") ==       0)
        {
          fetch(pmtr.name, &int_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << int_val;
        }
        else if (pmtr.type.compare("flt") ==       0) 
        {
          fetch(pmtr.name, &flt_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << std::setprecision(15) << std::scientific;
          ofs << flt_val;
        }
        else if (pmtr.type.compare("dbl") ==       0) 
        {
          fetch(pmtr.name, &dbl_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << std::setprecision(15) << std::scientific;
          ofs <<  dbl_val;
        }
        else if (pmtr.type.compare("log") ==       0) 
        {
          fetch(pmtr.name, &log_val);
          ofs << std::setw(max_str_val_len + 2) << std::left;
          ofs << log_val;
        }
        
        ofs   << std::setw(4)                   << pmtr.type;

        if (pmtr.adjustability.compare("sfa") == 0 ) { pmtr.adjustability.assign("sfx"); }
        ofs   << std::setw(4)                   << pmtr.adjustability;

        ofs   << std::endl;

     }

     ofs.close();
}

  /* ~ Destructor    ~ */

parameter_map::~parameter_map() {

  /* ~ not implemented ~ */

}
