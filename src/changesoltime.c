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

/**@file   changesoltime.c
 * @brief  Allow to change the solving time
 * @author Marc Pfetsch
 */

#include "changesoltime.h"

#include "scip/scip_timing.h"
#include "scip/clock.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/struct_clock.h"


/** change solving time */
void changeSolvingTime(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_Real             newtime             /**< new solving time */
   )
{
   SCIPclockSetTime(scip->stat->solvingtime, newtime);
}
