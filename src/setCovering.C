/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program   setCovering                         */
/*                                                                           */
/*    a simple branch-and-cut algorithm to solve the set covering problem    */
/*                                                                           */
/*    Copyright (C) 2006 - 2026  Discrete Optimization Group, TU Darmstadt   */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*    Copyright (C) 2002-2026 Tobias Achterberg                              */
/*                  2002-2026 Konrad-Zuse-Zentrum                            */
/*                               fuer Informationstechnik Berlin             */
/*       mailto: scip@zib.de                                                 */
/*       SCIP is distributed under the terms of the SCIP Academic Licence,   */
/*       see file COPYING in the SCIP distribution.                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   setCovering.C
 * @brief  Implementation of top class for handling optimization data
 * @author Marc Pfetsch
 */

#include "setCovering.h"

// SCIP includes:
extern "C" {
#include "sepa_bang.h"
}
#include <lpi/lpi.h>
#include <scip/scip.h>
#include <scip/cons_logicor.h>

// own includes:
#include "macro.h"
#include "branch_most_covered.h"
#include "branch_intersection.h"
#include "branch_half_constraint.h"
#include "branch_orbital.h"
#include "symmetry.h"
#include "event_sol.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/times.h>

// default values for parameters
#define DEFAULT_PRESOLVING    FALSE     //!< Perform presolving?
#define DEFAULT_PREADDITION   TRUE      //!< Perform addition presolving?
#define DEFAULT_PREMERGE      TRUE      //!< Perform merge presolving?
#define DEFAULT_PREFIXATION   TRUE      //!< Perform fixation presolving?
#define DEFAULT_MAXFACTOR     1000.0    //!< factor for early exit in addition
#define DEFAULT_MAXPREADD     1         //!< number of rows to add for each leader/follower combination (-1 = unbounded)


/**@brief Constructor
 * @param problem_name  name of the problem
 * @param file_name     name of the file containing the file_name
 * @param scip          instance of scip
 */
SETCOVERING::SETCOVERING(const char* problemName, const char* fileName, SCIP* scip)
{
   scip_ = scip;

   CHECK_OK( SCIPaddBoolParam(scip_, "sc/presolving", "Perform presolving?",
         &presolving_, TRUE, DEFAULT_PRESOLVING, nullptr, nullptr) );

   CHECK_OK( SCIPaddBoolParam(scip_, "sc/preaddition", "Perform addition presolving?",
         &preaddition_, TRUE, DEFAULT_PREADDITION, nullptr, nullptr) );

   CHECK_OK( SCIPaddBoolParam(scip_, "sc/premerge", "Perform merge presolving?",
         &premerge_, TRUE, DEFAULT_PREMERGE, nullptr, nullptr) );

   CHECK_OK( SCIPaddBoolParam(scip_, "sc/prefixation", "Perform fixation presolving?",
         &prefixation_, TRUE, DEFAULT_PREFIXATION, nullptr, nullptr) );

   CHECK_OK( SCIPaddRealParam(scip_, "sc/maxfactor", "factor for early exit in addition",
         &maxfactor_, TRUE, DEFAULT_MAXFACTOR, 0.0, SCIP_REAL_MAX, nullptr, nullptr) );

   CHECK_OK( SCIPaddIntParam(scip_, "sc/maxpreadd", "number of rows to add for each leader/follower combination (-1 = unbounded)",
         &maxpreadd_, TRUE, DEFAULT_MAXPREADD, -1, INT_MAX, nullptr, nullptr) );

   CHECK_OK( SCIPaddBoolParam(scip_, "sc/output_optimal_cover", "whether the optimal cover should be output at the end",
         &output_optimal_cover_, TRUE, FALSE, nullptr, nullptr) );

   CHECK_OK( SCIPaddBoolParam(scip_, "sc/writeBestSolution", "whether optimal solution should be written to a file",
         &writeBestSolution_, TRUE, FALSE, nullptr, nullptr) );

   (void) problemName_.assign(problemName);

   Aorig_ = nullptr;
   morig_ = 0;
   norig_ = 0;
   A_ = nullptr;
   m_ = 0;
   n_ = 0;
   fixed_ = nullptr;
   nfixed_ = -1;
   bestCover_ = nullptr;
   bestCoverSize_ = -1;
   bestCoverValue_ = -1.0;
   origVars_ = nullptr;
   objorig_ = nullptr;
   obj_ = nullptr;
   origVars_ = nullptr;
   vars_ = nullptr;
   incidentVars_ = nullptr;

   // read input
   CHECK_OK( readProblemLP(fileName) );

   std::cout << "\nNumber of rows:\t\t " << morig_ << std::endl;
   std::cout << "Number of variables:\t " << norig_ << std::endl;
   std::cout << std::endl;
}



//! Destructor
SETCOVERING::~SETCOVERING()
{
   if ( bestCover_ != nullptr )
      delete bestCover_;

   for (int i = 0; i < morig_; ++i)
      SCIPfreeBlockMemoryArray(scip_, &Aorig_[i], norig_);
   SCIPfreeBlockMemoryArray(scip_, &Aorig_, morig_);
   SCIPfreeBlockMemoryArray(scip_, &objorig_, norig_);

   for (int i = 0; i < m_; ++i)
      SCIPfreeBlockMemoryArray(scip_, &A_[i], n_);
   SCIPfreeBlockMemoryArray(scip_, &A_, m_);
   SCIPfreeBlockMemoryArray(scip_, &obj_, n_);
   SCIPfreeBlockMemoryArrayNull(scip_, &fixed_, n_);

   SCIPfreeBlockMemoryArrayNull(scip_, &incidentVars_, n_);
}


//! include plugins
SCIP_RETCODE SETCOVERING::includePlugins()
{
   SCIP_CALL( SCIPincludeObjEventhdlr(scip_, new EVENT_SOL(this), TRUE) );

   SCIP_CALL( SCIPincludeObjBranchrule(scip_, new BRANCH_MOST_COVERED(this), TRUE) );
   SCIP_CALL( SCIPincludeObjBranchrule(scip_, new BRANCH_INTERSECTION(this), TRUE) );
   SCIP_CALL( SCIPincludeObjBranchrule(scip_, new BRANCH_HALF_CONSTRAINT(this), TRUE) );
   SCIP_CALL( SCIPincludeBranchruleOrbital(scip_) );

   // constraint handlers:
   SCIP_CALL( SCIPincludeSepaBaNg(scip_) );

   return SCIP_OKAY;
}


/**@brief Initialize everything needed for the optimization, e.g., constraints, variables, etc.
 */
SCIP_RETCODE SETCOVERING::initializeOptimization(const char* paramfile)
{
   // compute incident constraints
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &incidentVars_, n_) );
   for (int j = 0; j < n_; ++j)
   {
      int cnt = 0;
      for (int i = 0; i < m_; ++i)
      {
	 if ( A_[i][j] == 1)
	    ++cnt;
      }
      incidentVars_[j] = cnt;
   }

   // set best primal solution (= all variables)
   bestCover_ = new Bitset(n_);
   bestCoverValue_ = 0.0;
   bestCoverSize_ = 0;
   for (int i = 0; i < n_; ++i)
   {
      if ( fixed_[i] != SCfixedzero )
      {
         (void) bestCover_->insert(i);
         bestCoverValue_ += obj_[i];
         ++bestCoverSize_;
      }
   }

   //----------------------------------------------------------------------------
   // initialize SCIP
   SCIP_CALL( SCIPsetObjsense(scip_, SCIP_OBJSENSE_MINIMIZE) );


   if ( paramfile != nullptr )
   {
      SCIP_CALL( SCIPreadParams(scip_, paramfile) );
      SCIP_CALL( SCIPwriteParams(scip_, nullptr, TRUE, FALSE) );
   }

   //----------------------------------------------------------------------------
   // initialize variables
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &origVars_, n_) );

   char name[SCIP_MAXSTRLEN];
   for (int i = 0; i < n_; ++i)
   {
      SCIP_Real lb = 0.0;
      SCIP_Real ub = 1.0;

      if ( fixed_[i] == SCfixedzero )
         ub = 0.0;
      else if ( fixed_[i] == SCfixedone )
         lb = 1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL( SCIPcreateVar(scip_, &origVars_[i], const_cast<const char*>(name), lb, ub, obj_[i], SCIP_VARTYPE_BINARY,
            TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip_, origVars_[i]) );
   }

   // generate covering constraints
   SCIP_VAR** vars;
   SCIP_CALL( SCIPallocBufferArray(scip_, &vars, n_) );
   for (int i = 0; i < m_; ++i)
   {
      int cnt = 0;
      for (int j = 0; j < n_; ++j)
      {
         assert( A_[i][j] == 0 || A_[i][j] == 1 );
         if ( A_[i][j] == 1 )
            vars[cnt++] = origVars_[j];
      }
      assert( cnt > 0 );

      SCIP_CONS* cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "C%d", i);
      SCIP_CALL( SCIPcreateConsLogicor(scip_, &cons, const_cast<const char*>(name), cnt, vars, FALSE, TRUE, TRUE, TRUE, TRUE,
            FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip_, cons) );
      SCIP_CALL( SCIPreleaseCons(scip_, &cons) );
   }
   SCIPfreeBufferArray(scip_, &vars);

   //SCIP_CALL( SCIPprintOrigProblem(scip_, 0) );

   // transform solution - do this explicitly since we want to submit our initial solution
   SCIP_CALL( SCIPtransformProb(scip_) );

   // tell SCIP about initial solution
   // either the trivial solution initialized in constructor or a user defined solution
   // set up betterCover().
   SCIP_CALL( submitSolution(0) );

   return SCIP_OKAY;
}





/**@brief Output statistics and Cleanup after the optimization process
 */
SCIP_RETCODE SETCOVERING::terminateOptimization()
{
   std::cout << "\n";

   // print Scip statistics
   SCIP_CALL( SCIPprintStatistics(scip_, NULL) );

   std::cout << "\n\nThe objective value is: " << bestCoverValue_ << "\n";

   if ( output_optimal_cover_ )
   {
      std::cout << "\n";
      assert( (SCIPgetStage(scip_) == SCIP_STAGE_SOLVED) || (SCIPgetStage(scip_) == SCIP_STAGE_SOLVING) );
      if ( SCIPgetStage(scip_) == SCIP_STAGE_SOLVED )
	 std::cout << "An optimal cover is: " << std::endl;
      else
	 std::cout << "The best cover is: " << std::endl;
      for (int j = 0; j < n_; ++j)
      {
	 if ( bestCover_->contains(j) )
	    std::cout << "1 ";
	 else
	    std::cout << "0 ";
      }
      std::cout << "\n";
   }

   if ( writeBestSolution_ )
   {
      std::string filename = problemName_ + ".sol";
      writeSolution(filename.c_str());
   }

   // free Problem
   SCIP_CALL( SCIPfreeTransform(scip_) );

   // release variables
   for (int i = 0; i < n_; ++i)
   {
      if ( origVars_[i] != nullptr )
      {
         SCIP_CALL( SCIPreleaseVar(scip_, &origVars_[i]) );
      }
   }
   SCIPfreeBlockMemoryArrayNull(scip_, &origVars_, n_);
   origVars_ = nullptr;

   return SCIP_OKAY;
}

#ifdef SCIP_MORE_DEBUG
//! check whether set is already present
static
SCIP_Bool isSetPresent(
   int                   m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   const int*            S,                  //!< set to be tested
   const int**           A,                  //!< set covering 0/1-matrix
   const SC_VARFIXING*   fixed               //!< array to mark fixed variables
   )
{
   // brute fore method to check whether set S is present already
   for (int i = 0; i < m; ++i)
   {
      bool equal = true;
      for (int j = 0; j < n; ++j)
      {
         if ( fixed[j] == SCunfixed )
         {
            if ( A[i][j] != S[j] )
            {
               equal = false;
               break;
            }
         }
      }

      if ( equal )
         return TRUE;
   }

#if 0
   printf("S: ");
   for (int j = 0; j < n; ++j)
      printf("%d ", S[j]);
   printf("\n\n");
   for (int i = 0; i < m; ++i)
   {
      for (int j = 0; j < n; ++j)
         printf("%d ", A[i][j]);
      printf("\n");
   }
#endif

   return FALSE;
}
#endif

//! add set to A
static
SCIP_RETCODE addSet(
   SCIP*                 scip,               //!< SCIP instance
   int&                  m,                  //!< current number of rows
   int&                  maxm,               //!< current maximal number of rows
   int                   n,                  //!< number of columns of A (variables)
   int**&                A,                  //!< matrix
   const int*            S                   //!< set to be tested
   )
{
   assert( A != nullptr );

   if ( m >= maxm )
   {
      int newsize = SCIPcalcMemGrowSize(scip, m + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &A, maxm, newsize) );
      maxm = newsize;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &A[m], n) );
   for (int j = 0; j < n; ++j)
   {
      A[m][j] = S[j];
   }
   ++m;

   return SCIP_OKAY;
}

//! perform addition operation
static
SCIP_RETCODE performAddition(
   SCIP*                 scip,               //!< SCIP instance
   SCIP_Real             maxfactor,          //!< factor for maximal number of rows
   int&                  m,                  //!< number of rows of A
   int&                  maxm,               //!< maximal number of rows
   int                   n,                  //!< number of columns of A (variables)
   int**&                A,                  //!< set covering 0/1-matrix
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const VecInt&         leaders,            //!< current leaders
   const VecVecInt&      followers,          //!< current followers
   SCIP_HASHTABLE*       hashtable,          //!< hash table of rows
   int                   maxpreadd,          //!< number of rows to add for each leader/follower combination (-1 = unbounded)
   bool&                 changed             //!< whether we have added a set
   )
{
   int* S;
   int oldm = m;

   changed = false;

   SCIP_CALL( SCIPallocBufferArray(scip, &S, n) );

   for (size_t l = 0; l < leaders.size(); ++l)
   {
      int leader = leaders[l];
      assert( 0 <= leader && leader < n );
      assert( fixed[leader] == SCunfixed );

      for (size_t f = 0; f < followers[l].size(); ++f)
      {
         int follower = followers[l][f];
         assert( 0 <= follower && follower < n );
         assert( fixed[follower] == SCunfixed );

         int nadded = 0;
         for (int i = 0; i < m; ++i)
         {
            if ( A[i][follower] != 0 )
            {
               BMScopyMemoryArray(S, A[i], n);
               S[follower] = 0;
               S[leader] = 1;

               if ( ! SCIPhashtableExists(hashtable, (void*) S) )
               {
#ifdef SCIP_MORE_DEBUG
                  assert( ! isSetPresent(m, n, S, const_cast<const int**>(A), fixed) );
#endif

                  SCIP_CALL( addSet(scip, m, maxm, n, A, S) );
                  SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) A[m-1]) );
#if 0
                  printf("Found new set: ");
                  for (int j = 0; j < n; ++j)
                     printf("%d ", S[j]);
                  printf("\n");
#endif
                  changed = true;

                  // early exit
                  if ( maxfactor > 0.0 && m > maxfactor * oldm )
                  {
                     SCIPfreeBufferArray(scip, &S);
                     return SCIP_OKAY;
                  }

                  // possibly exit early
                  ++nadded;
                  if ( maxpreadd >= 0 && nadded >= maxpreadd )
                     break;
               }
               else
               {
#ifdef SCIP_MORE_DEBUG
                  assert( isSetPresent(m, n, S, const_cast<const int**>(A), fixed) );
#endif
               }
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &S);

   return SCIP_OKAY;
}

//! perform fixing operation
static
SCIP_RETCODE performFixation(
   SCIP*                 scip,               //!< SCIP instance
   int&                  m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   int**                 A,                  //!< set covering 0/1-matrix
   SC_VARFIXING*         fixed,              //!< array to mark fixed variables
   const VecInt&         leaders,            //!< current leaders
   const VecVecInt&      followers,          //!< current followers
   int&                  nfixed              //!< number of fixed variables
   )
{
   nfixed = 0;

   SCIP_Bool* marked;
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, n) );

   for (size_t l = 0; l < leaders.size(); ++l)
   {
      int leader = leaders[l];
      assert( 0 <= leader && leader < n );
      assert( fixed[leader] == SCunfixed );

      // mark all followers
      BMSclearMemoryArray(marked, n);
      for (size_t f = 0; f < followers[l].size(); ++f)
      {
         int follower = followers[l][f];
         assert( 0 <= follower && follower < n );
         assert( fixed[follower] == SCunfixed );

         marked[follower] = TRUE;
      }

      for (int i = 0; i < m; ++i)
      {
         int j;
         for (j = 0; j < n; ++j)
         {
            if ( j != leader && A[i][j] != 0 && ! marked[j] )
               break;
         }

         // if condition is satisfied
         if ( j >= n )
         {
            printf("Performing fixation of variable %d to 1.\n", leader);

            fixed[leader] = SCfixedone;
            ++nfixed;

            // delete rows covered by leader
            for (int k = m - 1; k >= 0; --k)
            {
               if ( A[k][leader] != 0 )
               {
                  SCIPfreeBlockMemoryArray(scip, &A[k], n);
                  A[k] = A[m-1];
                  --m;
               }
            }

            break;
         }
      }
   }

   SCIPfreeBufferArray(scip, &marked);

   return SCIP_OKAY;
}

//! perform merge
static
SCIP_RETCODE performMerge(
   SCIP*                 scip,               //!< SCIP instance
   int&                  m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   int**&                A,                  //!< set covering 0/1-matrix
   SC_VARFIXING*         fixed,              //!< array to mark fixed variables
   SCIP_Real*            obj,                //!< objective
   int&                  nfixed              //!< number of fixed variables
   )
{
   nfixed = 0;

   // check all pairs
   for (int s = 0; s < n; ++s)
   {
      if ( fixed[s] != SCunfixed )
         continue;

      for (int t = s + 1; t < n; ++t)
      {
         if ( fixed[t] != SCunfixed )
            continue;

         // check for rows
         int i;
         for (i = 0; i < m; ++i)
         {
            if ( (A[i][s] != 0 && A[i][t] == 0) || (A[i][s] == 0 && A[i][t] != 0) )
            {
               // search for different row
               SCIP_Bool found = FALSE;
               for (int k = 0; k < m; ++k)
               {
                  if ( (A[i][s] != 0 && A[k][s] == 0 && A[k][t] != 0) || (A[i][t] != 0 && A[k][s] != 0 && A[k][t] == 0) )
                  {
                     // check whether rows are equal
                     int j;
                     for (j = 0; j < n; ++j)
                     {
                        if ( fixed[j] == SCunfixed )
                        {
                           if ( j != s && j != t )
                           {
                              if ( A[i][j] != A[k][j] )
                                 break;
                           }
                        }
                     }

                     if ( j >= n )
                     {
                        found = TRUE;
                        break;
                     }
                  }
               }

               if ( ! found )
                  break;
            }
         }

         // check second condition if we found corresponding rows for every row i
         if ( i >= m )
         {
            for (i = 0; i < m; ++i)
            {
               // common row
               if ( A[i][s] != 0 && A[i][t] != 0 )
               {
                  // search for row k containing s but not t and a row l containing t but not s
                  SCIP_Bool found = FALSE;
                  for (int k = 0; k < m; ++k)
                  {
                     if ( A[k][s] != 0 && A[k][t] == 0 )
                     {
                        for (int l = 0; l < m; ++l)
                        {
                           if ( A[l][s] == 0 && A[l][t] != 0 )
                           {
                              // check condition
                              int j;
                              for (j = 0; j < n; ++j)
                              {
                                 if ( fixed[j] == SCunfixed )
                                 {
                                    if ( j != s && j != t )
                                    {
                                       if ( A[i][j] != A[k][j] || A[k][j] != A[l][j] )
                                          break;
                                    }
                                 }
                              }

                              if ( j >= n )
                              {
                                 found = TRUE;
                                 break;
                              }
                           }
                        }

                        if ( found )
                           break;
                     }
                  }

                  if ( ! found )
                     break;
               }
            }

            // perform merge
            if ( i >= m )
            {
               printf("Performing merge: (%d, %d).\n", s, t);

               fixed[t] = SCfixedzero;
               obj[s] += obj[t];
               ++nfixed;

               // delete rows
               for (i = m - 1; i >= 0; --i)
               {
                  if ( A[i][s] == 0 && A[i][t] != 0 )
                  {
                     SCIPfreeBlockMemoryArray(scip, &A[i], n);
                     A[i] = A[m-1];
                     --m;
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}



//! get the key for the element
SCIP_DECL_HASHGETKEY(rowkey)
{
   /* the key is the element itself */
   return elem;
}

//! returns TRUE iff both keys are equal
SCIP_DECL_HASHKEYEQ(roweq)
{
   SCIP_Longint n = (SCIP_Longint) userptr;
   int* row1 = (int*) key1;
   int* row2 = (int*) key2;
   for (SCIP_Longint j = 0; j < n; ++j)
   {
      if ( row1[j] != row2[j] )
         return FALSE;
   }
   return TRUE;
}

//! returns hash key
SCIP_DECL_HASHKEYVAL(rowhashval)
{
   SCIP_Longint n = (SCIP_Longint) userptr;
   int* row = (int*) key;

   uint64_t val = 0;
   for (SCIP_Longint j = 0; j < n; ++j)
   {
      if ( row[j] != 0 )
         ++val;
   }

   // formula from https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
   for (SCIP_Longint j = 0; j < n; ++j)
   {
      if ( row[j] != 0 )
         val ^= j + 0x9e3779b9 + (val << 6) + (val >> 2);
  }
  return val;
}

//! Perform symmetry presolving
SCIP_RETCODE SETCOVERING::performSymmetryPresolving()
{
   // copy original matrix
   int** A;
   int maxm;
   int m;
   int n;

   m = morig_;
   maxm = 2 * m;
   n = norig_;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &A, maxm) );
   for (int i = 0; i < m; ++i)
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip_, &A[i], Aorig_[i], n) );
   }

   // init set of fixed variables
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &fixed_, n) );
   for (int i = 0; i < n; ++i)
      fixed_[i] = SCunfixed;

   // copy objective
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip_, &obj_, objorig_, n) );

   if ( presolving_ )
   {
      SCIPinfoMessage(scip_, nullptr, "Running presolving ...\n");

      // prepare hash table
      SCIP_HASHTABLE* hashtable;
      SCIP_CALL( SCIPhashtableCreate(&hashtable, SCIPblkmem(scip_), 2 * m * maxfactor_, rowkey, roweq, rowhashval, reinterpret_cast<void*>(n)) );

      // fill hashtable
      int nduplicate = 0;
      for (int i = m - 1; i >= 0; --i)
      {
         void* res = SCIPhashtableRetrieve(hashtable, (void*) A[i]);
         if ( res == nullptr )
         {
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) A[i]) );
         }
         else
         {
            // delete duplicate row
            SCIPfreeBlockMemoryArray(scip_, &A[i], n);
            A[i] = A[m-1];
            --m;
            ++nduplicate;
         }
      }
      assert( m == SCIPhashtableGetNElements(hashtable) );
      if ( nduplicate > 0 )
         SCIPinfoMessage(scip_, nullptr, "Removed %d duplicate rows.\n", nduplicate);

      bool changed = false;
      int round = 0;
      do
      {
         SCIPinfoMessage(scip_, nullptr, "-----------------------------------------------------------\n");
         SCIPinfoMessage(scip_, nullptr, "Round %d:\n", round);
         ++round;

         changed = false;
         if ( preaddition_ || prefixation_ )
         {
            // variables for local timing
            struct tms comptimer_beg;
            struct tms comptimer_end;

            // compute symmetries
            long double groupsize;
            VecVecUInt generators;
            (void) times(&comptimer_beg);
            SCIP_CALL( computeAutomorphismsSetCovering(scip_, m, n, const_cast<const int**>(A), fixed_, objorig_, FALSE, groupsize, generators) );
            (void) times(&comptimer_end);

            SCIPinfoMessage(scip_, nullptr, "Number of generators of Aut(G):\t\t%lu\n", generators.size());
            SCIPinfoMessage(scip_, nullptr, "Size of Aut(G):\t\t\t\t%.1Lg\n", groupsize);
            double t = (comptimer_end.tms_utime - comptimer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
            SCIPinfoMessage(scip_, nullptr, "Computing symmetry time:\t\t%4.2f\n\n", t);

            // do nothing if there are no generators
            if ( generators.size() == 0 )
            {
               SCIPinfoMessage(scip_, nullptr, "No symmetry left.\n");
               break;
            }

            VecInt leaders;
            VecVecInt followers;
            SCIP_CALL( computeLeaderFollowersSSTgensFilter(scip_, n, fixed_, generators, ORBIT_RULE_MAXORBIT, leaders, followers) );

            SCIPinfoMessage(scip_, nullptr, "Number of leaders:\t\t\t%zu\n", leaders.size());
            size_t nfollowers = 0;
            for (size_t i = 0; i < followers.size(); ++i)
               nfollowers += followers[i].size();
            SCIPinfoMessage(scip_, nullptr, "Number of followers:\t\t\t%lu\n", nfollowers);

            if ( prefixation_ )
            {
               int mold = m;
               int nfixed;
               SCIP_CALL( performFixation(scip_, m, n, A, fixed_, leaders, followers, nfixed) );
               SCIPinfoMessage(scip_, nullptr, "\nFixed %d variables to 1 and deleted %d rows by fixation operation.\n", nfixed, mold - m);
               if ( nfixed > 0 )
               {
                  nfixed_ += nfixed;
                  changed = true;

                  // exit early if there is only one unfixed variable left
                  if ( nfixed >= n - 1 )
                     break;
               }
            }

            if ( preaddition_ && ! changed )
            {
               int mold = m;
               SCIP_CALL( performAddition(scip_, maxfactor_, m, maxm, n, A, fixed_, leaders, followers, hashtable, maxpreadd_, changed) );
               SCIPinfoMessage(scip_, nullptr, "\nAdded %d rows using the addition operation.\n", m - mold);

               if ( maxfactor_ > 0.0 && m > maxfactor_ * mold )
               {
                  SCIPinfoMessage(scip_, nullptr, "Stopping addition with %d rows.\n", m);
                  break;
               }
            }
         }

         if ( premerge_ )
         {
            // perform merging
            int mold = m;
            int nfixed;
            SCIP_CALL( performMerge(scip_, m, n, A, fixed_, obj_, nfixed) );
            SCIPinfoMessage(scip_, nullptr, "\nFixed %d variables and deleted %d rows using the merge operation.\n", nfixed, mold - m);

            if ( nfixed > 0 )
            {
               nfixed_ += nfixed;
               changed = true;

               // exit early if there is only one unfixed variable left
               if ( nfixed >= n - 1 )
                  break;
            }
         }
      }
      while ( changed );

      SCIPinfoMessage(scip_, nullptr, "-----------------------------------------------------------\n");

      SCIPhashtablePrintStatistics(hashtable, SCIPgetMessagehdlr(scip_));
      SCIPhashtableFree(&hashtable);
   }

   // copy data (reuse rows)
   m_ = m;
   n_ = n;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &A_, m) );
   for (int i = 0; i < m; ++i)
      A_[i] = A[i];
   SCIPfreeBlockMemoryArray(scip_, &A, maxm);

   assert( n_ == norig_ );

   return SCIP_OKAY;
}






/**@brief SCIP found a better solution
 *
 * We capture the solution and store it in our local data structures.
 * This function is called by the event handler 'event_sol'.
 */
void SETCOVERING::foundBetterSolution()
{
   // get best solution
   SCIP_SOL* sol = SCIPgetBestSol(scip_);

   assert( SCIPisLE(scip_, SCIPgetSolOrigObj(scip_, sol), bestCoverValue_) );

   // copy solution
   Bitset S(n_);

   int cnt = 0;
   SCIP_Real objval = 0.0;
   for (int j = 0; j < n_; ++j)
   {
      if ( vars_[j] != nullptr && SCIPgetSolVal(scip_, sol, vars_[j]) > 0.5 )
      {
	 (void) S.insert(j);
	 objval += obj_[j];
         ++cnt;
      }
   }

   bestCoverValue_ = objval;
   bestCoverSize_ = cnt;
   *bestCover_ = S;
}



/**@brief Submit currently best stored solution to SCIP
 * @param heur   pointer to heuristic (or 0)
 */
SCIP_RETCODE SETCOVERING::submitSolution(SCIP_HEUR* heur)
{
   SCIP_SOL* sol;
   SCIP_CALL( SCIPcreateSol(scip_, &sol, heur) );
   for (Bitset::const_iterator it = bestCover_->begin(); it != bestCover_->end(); ++it)
   {
      SCIP_CALL( SCIPsetSolVal(scip_, sol, vars_[*it], 1.0) );
   }
   SCIP_Bool stored;
   SCIP_CALL( SCIPaddSolFree(scip_, &sol, &stored) );
   assert( stored == TRUE );

   return SCIP_OKAY;
}



/**@brief Store initial solution
 * @param S      solution to be stored
 * @param size   size of solution
 * @returns @c true if the solution is better than the previous solutions
 *
 * This function can be used to store initial solutions.
 * It should be called @b before initOptimization() is called!
 * (initOptimization submits the solutions to SCIP)
 *
 * We compute the value of the solution and store it if it is better than the previously
 * stored solutions. Hence, this function can be called several times.
 */
bool SETCOVERING::initialCover(const Bitset& S, int size)
{
   // compute value
   SCIP_Real value = 0.0;

   for (Bitset::const_iterator it = S.begin(); it != S.end(); ++it)
      value += obj_[*it];

   if ( value < bestCoverValue_ )
   {
      *bestCover_ = S;
      std::cout << "Found initial cover with value: " << value << std::endl << std::endl;
      bestCoverValue_= value;
      bestCoverSize_ = size;

      return(true);
   }
   return(false);
}


/**@brief Store solution if it is better than the previously found solutions
 * @param value  value of solution
 * @param S      solution to be stored
 * @param size   size of solution
 * @param heur   pointer to heuristic
 * @returns @c true if the primal solution given by @p S is the best solution so far (@p S is stored)
 *
 * @c heur can be zero, e.g., if IIS_CONSTRAINT_HANDLER::feasible() called this function.
 */
bool SETCOVERING::betterCover(SCIP_Real value, const Bitset& S, int size, SCIP_HEUR* heur)
{
   if ( value < bestCoverValue_ )
   {
      *bestCover_ = S;
      std::cout << "Found new better cover with value: " << value << std::endl;
      bestCoverValue_= value;
      bestCoverSize_ = size;

      CHECK_OK( submitSolution(heur) );
      return(true);
   }
   return(false);
}





/**@brief Write solution file
 * @param filename	name of the input file
 *
 * The format of the file is as follows:
 * - problem name
 * - number of variables
 * - {0/1} vector
 */
void SETCOVERING::writeSolution(const char* filename)
{
   std::ofstream file(filename, std::ios::out);
   if (!file)
   {
      std::cerr << "Error: could not open " << filename << std::endl;
      abort();
   }

   std::cout << "Writing solution file " << filename << " ..." << std::endl;
   file << problemName_.c_str() << std::endl;
   file << n_ << std::endl;

   // output solution
   for (int j = 0; j < n_; ++j)
   {
      if ( bestCover_->contains(j) )
	 file << "1 ";
      else
	 file << "0 ";
   }
   file << std::endl;
}



/**@brief Get next number from string s
 * @param s	string
 */
static
int getNextNumber(char*& s)
{
   // skip whitespace
   while ( isspace((unsigned char)*s) )
      ++s;
   int tmp = atoi(s);
   // number
   while ( (*s != 0) && (! isspace((unsigned char)*s)) )
      ++s;
   return tmp;
}



/**@brief Read solution file
 * @param filename	name of the input file
 *
 * The format is as specified in writeSolution()
 */
void SETCOVERING::readSolution(const char* filename)
{
   std::ifstream file(filename, std::ios::in);
   if (!file)
   {
      std::cerr << "Error: could not open " << filename << std::endl;
      abort();
   }

   std::cout << "Reading solution file " << filename << " ..." << std::endl;

   char* buffer = new char [65535];

   // get problem name
   (void) file.getline(buffer, 65534);

   assert( ! file.eof() && ! file.fail() );
   std::string name(buffer);

   if ( name != problemName_ )
   {
      std::cerr << filename << " - Error! The name (" << name << ") given in the solution file\n";
      std::cerr << "does not match the current problem name (" << problemName_ << ")." << std::endl;
      abort();
   }

   // init own solution
   Bitset S;
   int size = 0;
   SCIP_Real value = 0;

   // read number of variables
   (void) file.getline(buffer, 65543);
   char* s = buffer;
   assert( ! file.fail() );
   int n = getNextNumber(s);
   if ( n != norig_ )
   {
      std::cerr << filename << " - Error! The number of variables in the file does not match ";
      std::cerr << "the number of variables in the instance." << std::endl;
      abort();
   }

   // read solution itself
   (void) file.getline(buffer, 65543);
   s = buffer;
   assert( ! file.fail() );
   for (int i = 0; i < norig_; ++i)
   {
      int r = getNextNumber(s);
      assert( r == 0 || r == 1 );
      if ( r == 1)
      {
	 (void) S.insert(i);
	 value += obj_[i];
	 ++size;
      }
   }
   (void) betterCover(value, S, size, 0);

   delete [] buffer;
}







/**@brief Read problem file
 * @param filename	name of the input file
 *
 * The format is:
 * <number of rows> <number of variables>
 * 0/1 matrix
 */
SCIP_RETCODE SETCOVERING::readProblem(const char* filename)
{
   std::ifstream file(filename, std::ios::in);
   if (!file)
   {
      std::cerr << "Error: could not open " << filename << std::endl;
      abort();
   }

   std::cout << "Reading file " << filename << " ..." << std::endl;

   char* buffer = new char [65535];

   // read number of constraints and variables
   (void) file.getline(buffer, 65534);
   assert( ! file.eof() && ! file.fail() );
   char* s = buffer;
   morig_ = getNextNumber(s);
   norig_ = getNextNumber(s);
   assert( morig_ > 0 );
   assert( norig_ > 0 );

   // read objective coefficients
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &objorig_, norig_) );
   for (int j = 0; j < norig_; ++j)
      objorig_[j] = getNextNumber(s);

   // read matrix
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &Aorig_, morig_) );
   for (int i = 0; i < morig_; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &Aorig_[i], norig_) );
      (void) file.getline(buffer, 65543);
      s = buffer;
      assert( ! file.fail() );

      for (int j = 0; j < norig_; ++j)
      {
	 int r = getNextNumber(s);
	 assert( r == 0 || r == 1 );
	 Aorig_[i][j] = r;
      }
   }

   delete [] buffer;

   return SCIP_OKAY;
}

/**@brief Read problem file in LP format
 * @param filename	name of the input file
 */
SCIP_RETCODE SETCOVERING::readProblemLP(const char* filename)
{
   // read from file into lpi
   std::cout << "\nReading file: " << filename << " ... " << std::endl;
   SCIP_LPI* lp;
   SCIP_CALL( SCIPlpiCreate(&lp, NULL, filename, SCIP_OBJSEN_MINIMIZE) );
   SCIP_CALL( SCIPlpiReadLP(lp, filename) );

   // get number of rows and columns:
   int numnz;
   SCIP_CALL( SCIPlpiGetNRows(lp, &morig_) );
   SCIP_CALL( SCIPlpiGetNCols(lp, &norig_) );
   SCIP_CALL( SCIPlpiGetNNonz(lp, &numnz) );
   assert( morig_ > 0 );
   assert( norig_ > 0 );

   // get objective
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &objorig_, norig_) );
   SCIP_CALL( SCIPlpiGetObj(lp, 0, norig_ - 1, objorig_) );

   int* matbeg = new int [morig_ + 1];
   int* matcnt = new int [morig_];
   int* matind = new int [numnz];
   SCIP_Real* matval = new SCIP_Real [numnz];
   SCIP_Real* lhs = new SCIP_Real [morig_];
   SCIP_Real* rhs = new SCIP_Real [morig_];

   int cnt = 0;
   SCIP_CALL( SCIPlpiGetRows(lp, 0, morig_ - 1, lhs, rhs, &cnt, matbeg, matind, matval) );
   assert(cnt == numnz);
   matbeg[morig_] = cnt; /*lint !e2662*/

   // set up matrix
   SCIP_CALL( SCIPallocBlockMemoryArray(scip_, &Aorig_, morig_) );
   for (int i = 0; i < morig_; ++i)
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip_, &Aorig_[i], norig_) );

      // check sides
      if ( ! SCIPisInfinity(scip_, rhs[i]) )
      {
	 std::cerr << "Error: input file is not a set covering problem. ";
	 std::cerr << "The right hand side of a constraint is not infinity!" << std::endl;
	 abort();
      }
      if ( ! SCIPisEQ(scip_, lhs[i], 1.0) )
      {
	 std::cerr << "Error: input file is not a set covering problem. ";
	 std::cerr << "The left hand side of a constraint is not 1!" << std::endl;
	 abort();
      }

      // get matrix
      for (int j = matbeg[i]; j < matbeg[i+1]; ++j) /*lint !e2662*/
      {
	 if ( ! SCIPisEQ(scip_, matval[j], 1.0) )
	 {
	    std::cerr << "Error: input file is not a set covering problem. ";
	    std::cerr << "A nonzero entry is different from 1!" << std::endl;
	    abort();
	 }
	 Aorig_[i][matind[j]] = 1;
      }
   }

   SCIP_CALL( SCIPlpiFree(&lp) );

   delete [] rhs;
   delete [] lhs;
   delete [] matval;
   delete [] matind;
   delete [] matcnt;
   delete [] matbeg;

   return SCIP_OKAY;
}


/**@brief Create data for the transformed problem
 *
 * We simply call the corresponding SCIP function to obtain the variables
 * of the transformed problem.
 */
SCIP_RETCODE SETCOVERING::scip_trans(SCIP* scip, ObjProbData** objprobdata, SCIP_Bool* deleteobject)
{
   vars_ = new SCIP_VAR* [n_];
   for (int v = 0; v < n_; ++v)
   {
      if ( origVars_[v] != nullptr )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip_, origVars_[v], &(vars_[v])) );
         assert( vars_[v] != 0);
      }
      else
         vars_[v] = nullptr;
   }

   *objprobdata = this;
   *deleteobject = FALSE;
   return SCIP_OKAY;
}



/**@brief Delete data of the transformed problem
 *
 * Delete the vector 'vars_'
 */
SCIP_RETCODE SETCOVERING::scip_deltrans(SCIP* scip)
{
   delete [] vars_;
   return SCIP_OKAY;
}
