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

/**@file   event_sol.h
 * @brief  Event handler for new solutions
 * @author Marc Pfetsch
 */

#ifndef EVENT_SOL_H
#define EVENT_SOL_H


#include <objscip/objeventhdlr.h>


class SETCOVERING;




/**@class EVENT_SOL
 * @brief Event handler for new solutions
 */
class EVENT_SOL : public scip::ObjEventhdlr
{
public:

   //! @name Constructor/Destructor
   //@{
   EVENT_SOL(SETCOVERING* setcovering);
   virtual ~EVENT_SOL();
   //@}

   //! @name Virtual functions of SCIP to redefine
   //@{
   virtual SCIP_DECL_EVENTINIT(scip_init);
   virtual SCIP_DECL_EVENTEXEC(scip_exec);
   //@}

private:
   SETCOVERING* setcovering_;
};

#endif
