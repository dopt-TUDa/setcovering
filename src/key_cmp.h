// -*- C++ -*-
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program   min_iisc_scip                       */
/*                                                                           */
/*    an implementation of a branch-and-cut algorithm to solve the           */
/*    Min IIS Cover Problem / Maximum Feasible Subsystem Problem             */
/*                                                                           */
/*    Copyright (C) 1999 - 2026  Discrete Optimization Group, TU Darmstadt   */
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

/**@file   key_cmp.h
 * @brief  Class for comparing two indices by their keys
 * @author Marc Pfetsch
 */

#ifndef KEYCMP_H
#define KEYCMP_H


/**@class RealKeyCmp
 * @brief Class to compare two indices by their keys
 *
 * The keys are passed at construction time.
 * Note that no bound checking is performed!
 */
template <typename T>
class KeyCmp
{
public:
   KeyCmp(T* x) : x_(x) {};

   bool operator() ( const int &i, const int &j)
   {
      return (x_[i] > x_[j]);
   }

private:
   T* x_;
};


#endif
