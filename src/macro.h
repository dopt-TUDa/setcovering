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

/**@file   macro.h
 * @brief  defines usefull macros
 * @author Marc Pfetsch
 */

#ifndef MACRO_H
#define MACRO_H


#define CHECK_OK(x) { SCIP_RETCODE _restat_; if ( (_restat_ = (x)) != SCIP_OKAY ) { printf("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_); exit(1); } }


#define CHECK_PARAM_OK(x) { SCIP_RETCODE _restat_; if ( (_restat_ = (x)) != SCIP_OKAY && (_restat_ != SCIP_PARAMETERUNKNOWN)) { printf("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_); exit(1); } }

#endif
