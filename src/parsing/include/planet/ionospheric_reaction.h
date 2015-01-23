//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_IONOSPHERIC_PARSING_H
#define PLANET_IONOSPHERIC_PARSING_H

//Antioch
//Planet
#include "KineticsBranchingStructure.h"

//C++
#include <fstream>

namespace Planet
{
   void read_ionospheric_database();

   inline
   void read_ionospheric_database(const std::string &input_file)
   {
      std::ifstream database(input_file.c_str());

      database.close();
   }
}
#endif
