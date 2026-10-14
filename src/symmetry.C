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

/**@file   symmetry.C
 * @brief  Methods to deal with symmetries
 * @author Marc Pfetsch
 *
 * This code is based on the implementation in maxkcol by Christopher Hojny.
 */

#include "symmetry.h"


/* include bliss */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4189)  // local variable is initialized but not referenced
# pragma warning(disable: 4267)  // conversion of size_t to int (at sassy/preprocessor.h:2897)
# pragma warning(disable: 4388)  // compare signed and unsigned expression
# pragma warning(disable: 4456)  // shadowed variable
# pragma warning(disable: 4430)  // missing type specifier
#endif


#ifdef SC_WITH_NAUTY
extern "C" {
#include <nauty/nauty.h>
#include <nauty/nausparse.h>
}
#endif

#ifdef SC_WITH_BLISS
/* include bliss */
#include <bliss/defs.hh>
#include <bliss/graph.hh>
#endif

#ifdef SC_WITH_SBLISS
/* include bliss */
#include <bliss/defs.hh>
#include <bliss/graph.hh>
#if SCIP_VERSION >= 1000
/* include sassy (as part of dejavu) */
#include <dejavu/dejavu.h>
#include <dejavu/tools/bliss_converter.h>
#else
/* include sassy */
#include <sassy/preprocessor.h>
#include <sassy/tools/bliss_converter.h>
#endif
#endif

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wshadow"
#endif

#ifdef _MSC_VER
# pragma warning(pop)
#endif

//! Select an orbit and its leader from a list of orbits
static
SCIP_RETCODE selectOrbitAndLeader(
   SCIP*                 scip,               //!< SCIP instance
   const VecVecInt&      orbits,             //!< vector of orbits (each orbit has size > 1)
   int&                  leader,             //!< return selected leader
   size_t&               orbitidx,           //!< return selected orbit
   ORBIT_RULE            orbitrule           //!< rule to select orbit
   )
{
   leader = -1;
   orbitidx = SIZE_MAX;

   // select an orbit and its leader
   if ( orbitrule == ORBIT_RULE_MAXORBIT )
   {
      orbitidx = 0;
      leader = orbits[0][0];
      size_t orbitsize = orbits[orbitidx].size();

      for (size_t i = 1; i < orbits.size(); ++i)
      {
         if ( orbits[i].size() > orbitsize )
         {
            orbitsize = orbits[i].size();
            orbitidx = i;
            leader = orbits[i][0];
         }
      }
   }
   else if ( orbitrule == ORBIT_RULE_MINORBIT )
   {
      orbitidx = 0;
      leader = orbits[0][0];
      size_t orbitsize = orbits[orbitidx].size();

      for (size_t i = 1; i < orbits.size(); ++i)
      {
         assert( orbits[i].size() >= 2 );
         if ( orbits[i].size() < orbitsize )
         {
            orbitsize = orbits[i].size();
            orbitidx = i;
            leader = orbits[i][0];
         }
      }
   }
   else if ( orbitrule == ORBIT_RULE_FIRST )
   {
      // select the first variable of the first orbit as leader
      leader = orbits[0][0];
      orbitidx = 0;

      for (size_t i = 1; i < orbits[0].size(); ++i)
      {
         if ( orbits[0][i] < leader )
            leader = orbits[0][i];
      }
   }
   else
      abort(); // not yet implemented

   return SCIP_OKAY;
}

//! Compute leaders and followers for SST cuts heuristically by filtering generators
SCIP_RETCODE computeLeaderFollowersSSTgensFilter(
   SCIP*                 scip,               //!< SCIP instance
   int                   nelems,             //!< number of elements the permutations are acting on
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const VecVecUInt&     gens,               //!< group generators
   ORBIT_RULE            orbitrule,          //!< rule to choose next orbit
   VecInt&               leaders,            //!< returns leaders
   VecVecInt&            followers           //!< returns followers
   )
{
   assert( fixed != nullptr );

   // initialize data structures to store which generators are still active (based on filtering)
   int ngens = (int) gens.size();
   int nactivegens = ngens;
   std::vector<SCIP_Bool> genIsActive(ngens, TRUE);

   // to keep track which elements have already been added to an orbit in an iteration for a specific stabilizer
   std::vector<SCIP_Bool> elemadded(nelems, FALSE);

   // iteratively compute orbits as long as filtered stabilizer is non-trivial
   SCIP_Bool printorbsize = TRUE;
   while ( nactivegens > 0 )
   {
      // collect the orbits of the current stabilizer
      VecVecInt orbits;
      std::vector<int> var2orbit(nelems, -1);

      // compute orbits for each element of the ground set
      for (int v = 0; v < nelems; ++v)
      {
         VecInt orbit;

         // don't compute orbit of already treated elements
         if ( elemadded[v] || fixed[v] != SCunfixed )
            continue;

         orbit.push_back(v);
         elemadded[v] = TRUE;

         // iterate over all elements in the partially built orbit and compute their images
         size_t curpos = 0;
         while ( curpos < orbit.size() )
         {
            int curelem;

            curelem = orbit[curpos];

            // compute orbits using generators
            for (int i = 0; i < ngens && nactivegens > 0; ++i)
            {
               if ( ! genIsActive[i] )
                  continue;

               int image = (int) gens[i][curelem];

               // if we have found a new element of the orbit
               if ( ! elemadded[image] && fixed[image] == SCunfixed )
               {
                  orbit.push_back(image);
                  elemadded[image] = TRUE;
               }
            }
            ++curpos;
         }

         // add the orbit if it is non-trivial
         if ( orbit.size() > 1 )
         {
            orbits.push_back(orbit);
            for (VecInt::iterator it = orbit.begin(); it != orbit.end(); ++it)
               var2orbit[*it] = (int) orbits.size() - 1;
         }
      }

      if ( orbits.size() == 0 )
         return SCIP_OKAY;

      if ( printorbsize )
      {
         // compute average orbit size and percentage of variables on which symmetry group is active
         size_t orbitsum = 0;
         for (size_t i = 0; i < orbits.size(); ++i)
            orbitsum += orbits[i].size();
         SCIP_Real avorbsize = (orbits.size() > 0) ? orbitsum / (SCIP_Real) orbits.size() : 0.0;
         SCIP_Real symvarfrac = (nelems > 0) ? orbitsum / (SCIP_Real) nelems : 0.0;
         SCIPinfoMessage(scip, nullptr, "Average orbit size:\t\t\t%.1f\n", avorbsize);
         SCIPinfoMessage(scip, nullptr, "Fraction of variables with symmetry:\t%.1f\n", symvarfrac);
         printorbsize = FALSE;
      }

      // select the orbit and its leader
      int curleader;
      size_t curorbit;
      SCIP_CALL( selectOrbitAndLeader(scip, orbits, curleader, curorbit, orbitrule) );

      VecInt curfollowers;
      for (size_t i = 0; i < orbits[curorbit].size(); ++i)
      {
         if ( orbits[curorbit][i] != curleader )
            curfollowers.push_back(orbits[curorbit][i]);
      }

      // store orbit and leader
      followers.push_back(curfollowers);
      leaders.push_back(curleader);

      // clean auxiliary data structure
      for (int i = 0; i < nelems; ++i)
         elemadded[i] = FALSE;

      // filter generators
      for (int i = 0; i < ngens && nactivegens > 0; ++i)
      {
         if ( ! genIsActive[i] )
            continue;

         if ( (int) gens[i][curleader] != curleader )
         {
            genIsActive[i] = FALSE;
            --nactivegens;
         }
      }
   }

   return SCIP_OKAY;
}

//! Compute orbits taking filtered elements into account
SCIP_RETCODE computeOrbitsGensFilter(
   SCIP*                 scip,               //!< SCIP instance
   int                   n,                  //!< number of variables
   SCIP_Bool*            filtered,           //!< array marking variables to be filtered
   const VecVecUInt&     gens,               //!< group generators
   VecVecInt&            orbits              //!< orbits to be computed
   )
{
   assert( filtered != nullptr );

   orbits.clear();

   // initialize data structures to store which generators are still active (based on filtering)
   int ngens = (int) gens.size();
   int nactivegens = ngens;
   std::vector<SCIP_Bool> genIsActive(ngens, TRUE);

   // filter generators
   for (int i = 0; i < ngens && nactivegens > 0; ++i)
   {
      for (int v = 0; v < n; ++v)
      {
         if ( filtered[v] )
         {
            int image = (int) gens[i][v];
            if ( v != image )
            {
               genIsActive[i] = FALSE;
               --nactivegens;
               break;
            }
         }
      }
   }

   if ( nactivegens > 0 )
   {
      // to keep track which elements have already been added to an orbit in an iteration for a specific stabilizer
      std::vector<SCIP_Bool> elemadded(n, FALSE);

      // compute orbits for each element of the ground set
      for (int v = 0; v < n; ++v)
      {
         // don't compute orbit of already treated elements
         if ( elemadded[v] || filtered[v] )
            continue;

         VecInt orbit;
         orbit.push_back(v);
         elemadded[v] = TRUE;

         // iterate over all elements in the partially built orbit and compute their images
         size_t curpos = 0;
         while ( curpos < orbit.size() )
         {
            int curelem = orbit[curpos];

            // compute orbits using active generators
            for (int i = 0; i < ngens; ++i)
            {
               if ( ! genIsActive[i] )
                  continue;

               int image = (int) gens[i][curelem];

               // if we have found a new element of the orbit
               if ( ! elemadded[image] )
               {
                  assert( ! filtered[image] );
                  orbit.push_back(image);
                  elemadded[image] = TRUE;
               }
            }
            ++curpos;
         }

         // add the orbit if it is non-trivial
         if ( orbit.size() > 1 )
         {
            orbits.push_back(orbit);
         }
      }
   }

   return SCIP_OKAY;
}


/* To define the callback function we use the following global variables to store the information. */
static int n_ = 0;
static std::list<std::vector<unsigned int> > generator_list_ = {};        //!< generators of the graph automorphism group


#ifdef SC_WITH_NAUTY
/** auxiliary function that nauty calls when it finds a generator
 *
 *  If a permutation is found, we put it into the global list of generators.
 */
static
void writeGenerator(
   int,
   int*                  perm,               //!< Generator that nauty found
   int*,
   int,
   int,
   int                   nNodes              //!< Number of nodes in the graph
   )
{
#ifndef NDEBUG
   for (int i = 0; i < nNodes; ++i)
   {
      assert( perm[i] >= 0 );
   }
#endif

   /* check whether permutation is identity */
   SCIP_Bool isidentity = TRUE;
   for (int j = 0; j < n_; ++j)
   {
      if ( perm[j] != j )
         isidentity = FALSE;
   }

   /* don't store identity permutations */
   if ( ! isidentity )
   {
      // use iterator constructor: copy permutation given at perm; truncate at number of variables n_
      (void)generator_list_.emplace_back(std::vector<unsigned int>(perm, perm + n_) );
   }
}

//! Compute generators of automorphism group of a set covering instance with nauty
SCIP_RETCODE computeAutomorphismsSetCoveringNauty(
   SCIP*                 scip,               //!< SCIP instance
   int                   m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   const int**           A,                  //!< set covering 0/1-matrix
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const SCIP_Real*      weights,            //!< current weights
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   )
{
   assert( A != nullptr );
   assert( fixed != nullptr );

   // make sure that number of variables is stored
   n_ = n;

   // init options
   constexpr const char* mallocString = "malloc";
   static DEFAULTOPTIONS_SPARSEGRAPH(options);
   statsblk stats;

   // init callback functions for nauty (accumulate the group generators found by nauty)
   options.writeautoms = FALSE;
   options.userautomproc = writeGenerator;
   options.defaultptn = FALSE; /* use color classes */

   // memory allocation for nauty
   DYNALLSTAT(int, lab, lab_sz);
   DYNALLSTAT(int, ptn, ptn_sz);
   DYNALLSTAT(int, orbits, orbits_sz);

   sparsegraph sg;
   SG_INIT(sg);

   // compute degrees of nodes in nauty graph
   int nnodes = n + m;
   int* degrees;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &degrees, nnodes) );
   SCIP_Real* symweights;
   SCIP_CALL( SCIPallocBufferArray(scip, &symweights, n) );

   // compute degrees of columns and set weights
   int nfixed = 0;
   int nedges = 0;
   for (int v = 0; v < n; ++v)
   {
      if ( fixed[v] == SCunfixed )
      {
         symweights[v] = weights[v];
         assert( symweights[v] >= 0.0 );

         for (int i = 0; i < m; ++i)
         {
            if ( A[i][v] != 0 )
            {
               ++degrees[v];
               ++nedges;
               ++degrees[n + i];
               ++nedges;
            }
         }
      }
      else
         symweights[v] = - (SCIP_Real) (1 + nfixed++);  // mark fixed nodes with distinct negative numbers
   }

   // init space
   SG_ALLOC(sg, nnodes, nedges, mallocString); /*lint !e774*/
   nauty_check(WORDSIZE, nedges, nnodes, NAUTYVERSIONID);

   sg.nv = nnodes;
   sg.nde = (size_t) nedges;

   DYNALLOC1(int, lab, lab_sz, nnodes, mallocString);   /*lint !e727*/
   DYNALLOC1(int, ptn, ptn_sz, nnodes, mallocString);   /*lint !e727*/

   // fill in array with colors for variables
   for (int i = 0; i < nnodes; ++i)
      lab[i] = i;

   // sort the first part according to weights
   SCIPsortRealInt(symweights, lab, n);

   // fill in degrees and color classes
   std::vector<int> pos(nnodes);
   size_t cnt = 0;
   for (int j = 0; j < n; ++j)
   {
      sg.d[j] = degrees[j];   // degree of column j
      sg.v[j] = cnt;          // position of edges for column j
      pos[j] = (int) cnt;           // also store position
      cnt += (size_t) degrees[j];

      if ( j < n - 1 && SCIPisEQ(scip, symweights[j], symweights[j+1]) )
         ptn[j] = 1;  // color class does not end
      else
         ptn[j] = 0;  // color class ends
   }
   for (int i = 0; i < m; ++i)
   {
      int k = n + i;
      sg.d[k] = degrees[k];   // degree of row
      sg.v[k] = cnt;          // position of edges for row i
      pos[k] = (int) cnt;           // also store position
      cnt += (size_t) degrees[k];
      ptn[k] = 1;  // color class does not end
   }
   ptn[n + m - 1] = 0;  // end last color class
   SCIPfreeBufferArray(scip, &symweights);
   SCIPfreeBufferArray(scip, &degrees);

   // create edges in nauty graph
   for (int v = 0; v < n; ++v)
   {
      if ( fixed[v] == SCunfixed )
      {
         // overwrite pos in the following
         for (int i = 0; i < m; ++i)
         {
            if ( A[i][v] != 0 )
            {
               sg.e[pos[n + i]++] = v;
               sg.e[pos[v]++] = n + i;
            }
         }
      }
   }

   DYNALLOC1(int, orbits, orbits_sz, nnodes, mallocString); /*lint !e727*/

   // call nauty. Calls auxiliary method writeGenerator
   sparsenauty(&sg, lab, ptn, orbits, &options, &stats, nullptr);

   // store number of generators
   assert( stats.numgenerators >= 0 );

   // put generator_list into generators and transposition generators into transpositions
   std::list<std::vector<unsigned int> >::const_iterator lit, lend = generator_list_.end();
   for (lit = generator_list_.begin(); lit != lend; ++lit)
   {
#ifndef NDEBUG
      // whether the output generator contains all numbers 0, ..., n-1
      std::vector<bool> containsZeroToN(n, false);
      for (size_t i = 0; i < lit->size(); ++i)
         containsZeroToN[ (*lit)[i] ] = true;
      for (int i = 0; i < n; ++i)
         assert( containsZeroToN[i] );
#endif

      // store permutation
      generators.push_back(*lit);
   }

   SG_FREE(sg); /*lint !e774*/

   // too many graph symmetry generators ...
   if ( generators.size() > (size_t) n * (size_t) n )
   {
      SCIPwarningMessage(scip, "The number of generators of Aut(G) is large: %lu.\n\n", generators.size());
   }

   // get the size of the group (stats.grapsize1 is a double)
   groupsize = stats.grpsize1 * exp10l(stats.grpsize2);

   return SCIP_OKAY;
}
#endif

#ifdef SC_WITH_BLISS
/** callback function for bliss */
static
void blisshook(
   void*                 user_param,         /**< parameter supplied at call to bliss */
   unsigned int          dim,                /**< size of aut vector */
   const unsigned int*   aut                 /**< automorphism */
   )
{
   assert( aut != nullptr );

   /* check whether permutation is identity */
   SCIP_Bool isidentity = TRUE;
   for (int j = 0; j < n_; ++j)
   {
      if ( aut[j] != j )
         isidentity = FALSE;
   }

   if ( ! isidentity )
   {
      (void)generator_list_.emplace_back(std::vector<unsigned int>(aut, aut + n_) );
   }
}

//! Compute generators of automorphism group of a set covering instance with bliss
static
SCIP_RETCODE computeAutomorphismsSetCoveringBliss(
   SCIP*                 scip,               //!< SCIP instance
   int                   m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   const int**           A,                  //!< set covering 0/1-matrix
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const SCIP_Real*      weights,            //!< current weights
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   )
{
   assert( A != nullptr );
   assert( fixed != nullptr );

   // make sure that number of variables is stored
   n_ = n;

   // compute degrees of nodes
   int* idx;
   SCIP_CALL( SCIPallocBufferArray(scip, &idx, n) );
   SCIP_Real* symweights;
   SCIP_CALL( SCIPallocBufferArray(scip, &symweights, n) );

   // collect weights
   int nfixed = 0;
   for (int v = 0; v < n; ++v)
   {
      idx[v] = v;

      if ( fixed[v] == SCunfixed )
      {
         symweights[v] = weights[v];
         assert( symweights[v] >= 0.0 );
      }
      else
      {
         symweights[v] = -1.0;
         ++nfixed;
      }
   }

   // sort weights
   SCIPsortRealInt(symweights, idx, n);

   // determine colors and degrees
   int nnodes = n + m;
   int* colors;
   SCIP_CALL( SCIPallocBufferArray(scip, &colors, nnodes) );
   int color = 0;
   SCIP_Real lastweight = symweights[0];
   int cnt = 0;

   // colors of fixed nodes come first with distinct values, then all other nodes with values depending on distinct weights
   for (int j = 0; j < n; ++j)
   {
      assert( idx[j] >= 0 );
      int v = idx[j];
      assert( v < n );
      if ( fixed[v] == SCunfixed )
      {
         assert( symweights[v] >= 0.0 );
         if ( ! SCIPisEQ(scip, symweights[j], lastweight) )
            ++color;
         lastweight = symweights[j];
         colors[v] = color + nfixed; // shift colors
      }
      else
      {
         assert( symweights[v] < 0.0 );
         colors[v] = cnt++;
      }
   }
   assert( cnt == nfixed );

   // fill colors for rows
   for (int i = 0; i < m; ++i)
      colors[n + i] = n;

   // create bliss graph
   bliss::Graph blissgraph(0);

   // create nodes
   for (int i = 0; i < nnodes; ++i)
   {
      (void) blissgraph.add_vertex(colors[i]);
   }

   // fill in edges
   for (int v = 0; v < n; ++v)
   {
      if ( fixed[v] == SCunfixed )
      {
         for (int i = 0; i < m; ++i)
         {
            if ( A[i][v] != 0 )
               blissgraph.add_edge(v, n + i);
         }
      }
   }

   SCIPfreeBufferArray(scip, &colors);
   SCIPfreeBufferArray(scip, &symweights);
   SCIPfreeBufferArray(scip, &idx);

   /* compute automorphisms */
   bliss::Stats stats;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   blissgraph.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   blissgraph.set_component_recursion(false);

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and pass it to the blisshook above */
   auto reportglue = [&](unsigned int dim, const unsigned int* aut) {
      blisshook(nullptr, dim, aut);
   };

   /* start search */
   blissgraph.find_automorphisms(stats, reportglue, nullptr);
#else

   /* start search */
   blissgraph.find_automorphisms(stats, blisshook, nullptr);
#endif

#ifdef SCIP_OUTPUT
   (void) stats.print(stdout);
#endif

   // determine group size
   groupsize = stats.get_group_size_approx();

   // put generator_list into generators
   std::list<std::vector<unsigned int> >::const_iterator lit, lend = generator_list_.end();
   for (lit = generator_list_.begin(); lit != lend; ++lit)
   {
#ifndef NDEBUG
      // whether the output generator contains all numbers 0, ..., n-1
      std::vector<bool> containsZeroToN(n, false);
      for (size_t i = 0; i < lit->size(); ++i)
         containsZeroToN[ (*lit)[i] ] = true;
      for (int i = 0; i < n; ++i)
         assert( containsZeroToN[i] );
#endif

      // store permutation
      generators.push_back(*lit);
   }

   return SCIP_OKAY;
}
#endif


#ifdef SC_WITH_SBLISS
/** callback function for sassy */
static
void sassyhook(
   int                   dim,                /**< size of aut vector */
   const int*            aut                 /**< automorphism */
   )
{
   assert( aut != nullptr );

   /* check whether permutation is identity */
   SCIP_Bool isidentity = TRUE;
   for (int j = 0; j < n_; ++j)
   {
      if ( aut[j] != j )
         isidentity = FALSE;
   }

   if ( ! isidentity )
   {
      generator_list_.emplace_back(std::vector<unsigned int>(aut, aut + n_) );
   }
}

//! Compute generators of automorphism group of a graph with bliss and preprocessing by sassy
static
SCIP_RETCODE computeAutomorphismsSetCoveringSBliss(
   SCIP*                 scip,               //!< SCIP instance
   int                   m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   const int**           A,                  //!< set covering 0/1-matrix
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const SCIP_Real*      weights,            //!< current weights
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   )
{
   assert( A != nullptr );
   assert( fixed != nullptr );

   // make sure that number of variables is stored
   n_ = n;

   // compute degrees of nodes
   int* idx;
   SCIP_CALL( SCIPallocBufferArray(scip, &idx, n) );
   SCIP_Real* symweights;
   SCIP_CALL( SCIPallocBufferArray(scip, &symweights, n) );

   // collect weights
   int nfixed = 0;
   for (int v = 0; v < n; ++v)
   {
      idx[v] = v;

      if ( fixed[v] == SCunfixed )
      {
         symweights[v] = weights[v];
         assert( symweights[v] >= 0.0 );
      }
      else
      {
         symweights[v] = -1.0;
         ++nfixed;
      }
   }

   // sort weights
   SCIPsortRealInt(symweights, idx, n);

   // determine colors and degrees
   int nnodes = n + m;
   int* colors;
   int* degrees;
   SCIP_CALL( SCIPallocBufferArray(scip, &colors, nnodes) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &degrees, nnodes) );
   int color = 0;
   SCIP_Real lastweight = symweights[0];
   int cnt = 0;
   int nedges = 0;

   // colors of fixed nodes come first with distinct values, then all other nodes with values depending on distinct weights
   for (int j = 0; j < n; ++j)
   {
      assert( idx[j] >= 0 );
      int v = idx[j];
      assert( v < n );
      if ( fixed[v] == SCunfixed )
      {
         assert( symweights[j] >= 0.0 );
         if ( ! SCIPisEQ(scip, symweights[j], lastweight) )
            ++color;
         lastweight = symweights[j];
         colors[v] = color + nfixed; // shift colors

         for (int i = 0; i < m; ++i)
         {
            if ( A[i][v] != 0 )
            {
               ++degrees[v];
               ++degrees[n + i];
               ++nedges;
            }
         }
      }
      else
      {
         assert( symweights[j] < 0.0 );
         colors[v] = cnt++;
      }
   }
   assert( cnt == nfixed );

   for (int i = 0; i < m; ++i)
      colors[n + i] = n;

   // create dejavu graph
#if SCIP_VERSION >= 1000
   dejavu::static_graph sassygraph;
#else
   sassy::static_graph sassygraph;
#endif

   sassygraph.initialize_graph((unsigned int) nnodes, (unsigned int) nedges);

   // create nodes
   for (int i = 0; i < nnodes; ++i)
   {
      (void) sassygraph.add_vertex(colors[i], degrees[i]);
   }

   // fill in edges of graph
   for (int v = 0; v < n; ++v)
   {
      if ( fixed[v] == SCunfixed )
      {
         for (int i = 0; i < m; ++i)
         {
            if ( A[i][v] != 0 )
               sassygraph.add_edge((unsigned int) v, (unsigned int) (n + i));
         }
      }
   }

   SCIPfreeBufferArray(scip, &degrees);
   SCIPfreeBufferArray(scip, &colors);
   SCIPfreeBufferArray(scip, &symweights);
   SCIPfreeBufferArray(scip, &idx);

#if SCIP_VERSION >= 1000
   /* set up sassy preprocessor */
   dejavu::preprocessor sassy;

   /* lambda function to have access to data and pass it to sassyhook above */
   dejavu_hook sassyglue = [&](int dim, const int* p, int nsupp, const int* suppa) {
      sassyhook(dim, p);
   };
#else
   /* set up sassy preprocessor */
   sassy::preprocessor sassy;

   /* turn off some preprocessing that generates redudant permutations */
   sassy::configstruct sconfig;
   sconfig.CONFIG_PREP_DEACT_PROBE = true;
   sassy.configure(&sconfig);

   /* lambda function to have access to data and pass it to sassyhook above */
   sassy::sassy_hook sassyglue = [&](int dim, const int* p, int nsupp, const int* suppa) /*lint --e{715}*/ {
      sassyhook(dim, p);
   };
#endif

   /* call sassy to reduce graph */
   sassy.reduce(&sassygraph, &sassyglue);

   if ( sassygraph.get_sgraph()->v_size == 0 )
   {
#if SCIP_VERSION >= 1000
      dejavu::big_number grp_sz = sassy.grp_sz;
      groupsize = grp_sz.mantissa * powl(10.0, (SCIP_Real) grp_sz.exponent);
#else
      groupsize = sassy.base * powl(10.0, sassy.exp);
#endif
   }
   else
   {
      /* create bliss graph */
      bliss::Graph blissgraph(0);

      /* convert sassy to bliss graph */
#if SCIP_VERSION >= 1000
      convert_dejavu_to_bliss(&sassygraph, &blissgraph);
#else
      convert_sassy_to_bliss(&sassygraph, &blissgraph);
#endif

      /* compute automorphisms */
      bliss::Stats stats;

      /* Prefer splitting partition cells corresponding to variables over those corresponding
       * to inequalities. This is because we are only interested in the action
       * of the automorphism group on the variables, and we need a base for this action */
      blissgraph.set_splitting_heuristic(bliss::Graph::shs_f);
      /* disable component recursion as advised by Tommi Junttila from bliss */
      blissgraph.set_component_recursion(false);

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
      auto hook = [&](unsigned int dim, const unsigned int* aut) {
         sassy.bliss_hook(dim, aut);
      };

      /* start search */
      blissgraph.find_automorphisms(stats, hook, nullptr);
#else

      /* start search */
      blissgraph.find_automorphisms(stats, sassy::preprocessor::bliss_hook, (void*) &sassy);
#endif

      // determine group size
#if SCIP_VERSION >= 1000
      /* determine log10 of symmetry group size */
      dejavu::big_number grp_sz = sassy.grp_sz;
      groupsize = stats.get_group_size_approx() * grp_sz.mantissa * powl(10.0, (SCIP_Real) grp_sz.exponent);
#else
      groupsize = stats.get_group_size_approx() * sassy.base * powl(10.0, sassy.exp);
#endif
   }

   // put generator_list into generators
   std::list<std::vector<unsigned int> >::const_iterator lit, lend = generator_list_.end();
   for (lit = generator_list_.begin(); lit != lend; ++lit)
   {
#ifndef NDEBUG
      // whether the output generator contains all numbers 0, ..., n-1
      std::vector<bool> containsZeroToN(n, false);
      for (size_t i = 0; i < lit->size(); ++i)
         containsZeroToN[ (*lit)[i] ] = true;
      for (int i = 0; i < n; ++i)
         assert( containsZeroToN[i] );
#endif

      // store permutation
      generators.push_back(*lit);
   }

   return SCIP_OKAY;
}
#endif


#ifndef NDEBUG
//! check whether generators are actually automorphisms
static
void checkAutomorphismsSetCovering(
   SCIP*                 scip,               //!< SCIP instance
   int                   m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   const int**           A,                  //!< set covering 0/1-matrix
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const VecVecUInt&     generators,         //!< returns generators
   SCIP_Bool             silent              //!< whether no output should be produced
   )
{
   assert( A != nullptr );
   assert( fixed != nullptr );

   for (size_t g = 0; g < generators.size(); ++g)
   {
      const std::vector<unsigned int>& gen = generators[g];
      assert( (int) gen.size() == n );

      // loop over generator to check that fixed nodes are invariant
      for (int j = 0; j < n; ++j)
      {
         int img = (int) gen[j];
         assert( fixed[j] == fixed[img] );

         if ( fixed[j] != SCunfixed )
         {
            // perform a weak test: count number of nonzeros
            int cnt = 0;
            for (int i = 0; i < m; ++i)
            {
               if ( A[i][j] != 0 )
                  ++cnt;
            }

            int cntimg = 0;
            for (int i = 0; i < m; ++i)
            {
               if ( A[i][img] != 0 )
                  ++cntimg;
            }
            assert( cnt == cntimg );
         }
      }
   }
   if ( ! silent )
      SCIPinfoMessage(scip, nullptr, "Symmetry check passed.\n");
}
#endif


//! Compute generators of automorphism group of a set covering instance
SCIP_RETCODE computeAutomorphismsSetCovering(
   SCIP*                 scip,               //!< SCIP instance
   int                   m,                  //!< number of rows of A
   int                   n,                  //!< number of columns of A (variables)
   const int**           A,                  //!< set covering 0/1-matrix
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const SCIP_Real*      weights,            //!< current weights
   SCIP_Bool             silent,             //!< whether no output should be produced
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   )
{
   generators.clear();
   generator_list_.clear();
   groupsize = 0.0;

#ifdef SC_WITH_NAUTY
   SCIPinfoMessage(scip, nullptr, "Using nauty to compute symmetries ...\n");
   SCIP_CALL( computeAutomorphismsSetCoveringNauty(scip, m, n, A, fixed, weights, groupsize, generators) );
#endif

#ifdef SC_WITH_BLISS
   SCIPinfoMessage(scip, nullptr, "Using bliss to compute symmetries ...\n");
   SCIP_CALL( computeAutomorphismsSetCoveringBliss(scip, m, n, A, fixed, weights, groupsize, generators) );
#endif

#ifdef SC_WITH_SBLISS
   SCIPinfoMessage(scip, nullptr, "Using sassy + bliss to compute symmetries ...\n");
   SCIP_CALL( computeAutomorphismsSetCoveringSBliss(scip, m, n, A, fixed, weights, groupsize, generators) );
#endif

#ifndef NDEBUG
   checkAutomorphismsSetCovering(scip, m, n, A, fixed, generators, silent);
#endif

   return SCIP_OKAY;
}
