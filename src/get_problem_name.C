/**@file   get_problem_name.C
 * @brief  Implementation of function to get the problem name from a filename
 * @author Marc Pfetsch
 */

#include "get_problem_name.h"
#include <cassert>


std::string get_problem_name(const char* filename)
{
   // find last part of filename (used as problem name)
   int i = 0;
   while ( filename[i] != 0) // find end of string
      ++i;
   int l = i;                // end of string

   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/')) // find last part
      --i;
   assert(i > 0);

   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/'))   // find last part
	 --i;
   }

   if (filename[i] == '/')
      ++i;

   std::string problem_name;
   problem_name.reserve(l+1);
   while ( (i < l) && (filename[i] != 0) )
      problem_name.push_back(filename[i++]);

   return problem_name;
}
