// Generate a random set covering instance.
//
// The entries in the matrix are 1 with probability @a p supplied by the user.
//
// 8/2006

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>


//--------------------------------------------------------------------
// Parameters = m n p
//
// - m    : number of inequalities
// - n    : number of variables
// - p    : probability for having a 1
// - seed : seed value
int main(int argc, const char** argv)
{
   if (argc != 5)
   {
      std::cerr << "Usage: " << argv[0] << "  <n> <m> <p> <seed>" << std::endl;
      exit(1);
   }
   int n = atoi(argv[1]);
   int m = atoi(argv[2]);
   double p = atof(argv[3]);
   int seed = atoi(argv[4]);

   if ( n <= 1 || m <= 1 )
   {
      std::cerr << "Error: Input sizes are out of range." << std::endl;
      exit(1);
   }

   if ( p < 0.0 || p > 1.0 )
   {
      std::cerr << "Error: Probability value should be in the interval [0,1]." << std::endl;
      exit(1);
   }

   srand48(seed);

   //--- set up data structures
   int** A = new int * [m];
   for (int i = 0; i < m; i++)
      A[i] = new int [n];

   // generate problem:
   for (int i = 0; i < m; i++)
   {
      for (int j = 0; j < n; j++)
      {
	 double rnd = drand48();
	 if ( rnd < p )
	    A[i][j] = 1;
	 else
	    A[i][j] = 0;
      }
   }

   // check for zero rows
   bool* zeroRow = new bool [m];
   int nRows = 0;
   for (int i = 0; i < m; ++i)
   {
      int j = 0;
      for (j = 0; j < n; ++j)
      {
	 assert( A[i][j] == 0 || A[i][j] == 1 );
	 if ( A[i][j] == 1 )
	    break;
      }
      if (j >= n)
	 zeroRow[i] = true;
      else
      {
	 zeroRow[i] = false;
	 ++nRows;
      }
   }

   // output problem
   std::stringstream filename;
   filename << "sc." << nRows << "." << n << "." << p << "." << seed << ".lp" << std::ends;
   std::ofstream file(filename.str().c_str());

   file << "Minimize" << std::endl;
   file << " obj: x1";
   for (int j = 1; j < n; ++j)
      file << " + x" << j+1;
   file << std::endl;
   file << "Subject To" << std::endl;
   for (int i = 0; i < m; ++i)
   {
      if ( ! zeroRow[i] )
      {
	 bool first = true;
	 for (int j = 0; j < n; ++j)
	 {
	    if (A[i][j] != 0)
	    {
	       if ( ! first )
		  file << " + ";
	       file << "x" << j+1;
	       first = false;
	    }
	 }
	 file << " >= 1" << std::endl;
      }
   }
   /*
   file << "Bounds" << std::endl;
   for (int j = 0; j < n; ++j)
      file << "0 <= x" << j+1 << " <= 1" << std::endl;
   */
   file << "Binary" << std::endl;
   for (int j = 0; j < n; ++j)
      file << "x" << j+1 << std::endl;
   file << "End" << std::endl;

   std::cout << "Problem written to: " << filename.str() << std::endl;

   delete [] zeroRow;
   for (int i = 0; i < m; i++)
      delete [] A[i];
   delete [] A;
}
