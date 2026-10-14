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

/**@file   type_presolve.h
 * @brief  definitions for graph presolving
 * @author Marc Pfetsch
 */

#ifndef TYPE_PRESOLVE_H
#define TYPE_PRESOLVE_H

//! enum to encode node fixings
enum SC_Varfixing
{
   SCfixedzero  = -1,    /**< variable is fixed to 0 */
   SCunfixed    = 0,     /**< variable is unfixed */
   SCfixedone   = 1      /**< variable is fixed to 1 */
};
typedef enum SC_Varfixing SC_VARFIXING;

#endif
