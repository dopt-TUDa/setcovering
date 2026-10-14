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

/**@file   branch_orbital.h
 * @brief  Implementation of a orbital branching
 * @author Marc Pfetsch
 */

#ifndef BRANCH_ORBITAL_H
#define BRANCH_ORBITAL_H

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** creates the orbital branching rule and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBranchruleOrbital(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
