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

/**@file   changesoltime.h
 * @brief  Allow to change the solving time
 * @author Marc Pfetsch
 */

#ifndef CHANGESOLTIME_H
#define CHANGESOLTIME_H

#include "scip/def.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** change solving time */
void changeSolvingTime(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_Real             newtime             /**< new solving time */
   );

#ifdef __cplusplus
}
#endif

#endif
