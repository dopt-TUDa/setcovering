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

/**@file   type_symmetry.h
 * @brief  type definitions for symmetry
 * @author Marc Pfetsch
 */

#ifndef TYPE_SYMMETRY_H
#define TYPE_SYMMETRY_H

/** different types of orbit selection rule */
enum ORBIT_Rule
{
   ORBIT_RULE_MAXORBIT         = 0,                 //!< rule to select orbit of maximum size
   ORBIT_RULE_MINORBIT         = 1,                 //!< rule to select orbit of minimum size
   ORBIT_RULE_STRINGENT        = 2,                 //!< rule to select orbits to generate stringent SST cuts
   ORBIT_RULE_FIRST            = 3                  //!< select the first orbit
};
typedef enum ORBIT_Rule ORBIT_RULE;

#endif
