// -*- C++ -*-
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

/**@file   symmetry.h
 * @brief  Methods to deal with symmetries
 * @author Marc Pfetsch
 *
 * This code is based on the implementation in maxkcol by Christopher Hojny.
 */

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "scip/scip.h"
#include "type_presolve.h"
#include "type_symmetry.h"

#include <vector>
#include <list>

// typedefs
typedef std::vector<std::vector<int> >            VecVecInt;
typedef std::vector<std::vector<unsigned int> >   VecVecUInt;
typedef std::vector<int>                          VecInt;
typedef std::list<std::pair<int,int> >            IntPairList;

//! Compute leaders and followers for SST cuts heuristically by filtering generators
SCIP_EXPORT
SCIP_RETCODE computeLeaderFollowersSSTgensFilter(
   SCIP*                 scip,               //!< SCIP instance
   int                   nelems,             //!< number of elements the permutations are acting on
   const SC_VARFIXING*   fixed,              //!< array to mark fixed variables
   const VecVecUInt&     gens,               //!< group generators
   ORBIT_RULE            orbitrule,          //!< rule to choose next orbit
   VecInt&               leaders,            //!< returns leaders
   VecVecInt&            followers           //!< returns followers
   );

//! Compute orbits taking filtered elements into account
SCIP_RETCODE computeOrbitsGensFilter(
   SCIP*                 scip,               //!< SCIP instance
   int                   n,                  //!< number of variables
   SCIP_Bool*            filtered,           //!< array marking variables to be filtered
   const VecVecUInt&     gens,               //!< group generators
   VecVecInt&            orbits              //!< orbits to be computed
   );

//! Compute generators of automorphism group of a set covering instance
SCIP_EXPORT
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
   );

#endif
