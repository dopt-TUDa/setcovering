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

/**@file   sepa_bang.h
 * @brief  Balas/Ng cut separator
 * @author Marc Pfetsch
 */

#ifndef SEPA_BANG_H
#define SEPA_BANG_H


#include "scip/scip.h"


/** creates the Balas/Ng cut separator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeSepaBaNg(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif
