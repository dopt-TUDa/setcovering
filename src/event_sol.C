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
//#define SCIP_DEBUG
/**@file   event_sol.C
 * @brief  Implementation of the handler for new solutions
 * @author Marc Pfetsch
 */


#include "event_sol.h"

#include "setCovering.h"
#include "macro.h"


#define EVENTHDLR_NAME         "newsol"
#define EVENTHDLR_DESC         "new solutions event handler"




EVENT_SOL::EVENT_SOL(SETCOVERING* setcovering) :
   ObjEventhdlr(setcovering->scip(), EVENTHDLR_NAME, EVENTHDLR_DESC)
{
   setcovering_ = setcovering;
}



EVENT_SOL::~EVENT_SOL()
{
}


/** execution method of event handler
 *
 *  Processes the event. The method is called every time an event occurs, for which the event handler
 *  is responsible. Event handlers may declare themselves resposible for events by calling the
 *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
 *  given event handler and event data.
 */
SCIP_DECL_EVENTEXEC(EVENT_SOL::scip_exec)
{
   assert( scip != 0 );
   assert( eventhdlr != 0 );

   SCIPdebugMessage("Event handler newsol...\n");

   if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
      setcovering_->foundBetterSolution();

   return SCIP_OKAY;
}



//! Initialize event
SCIP_DECL_EVENTINIT(EVENT_SOL::scip_init)
{
   assert( scip != 0);
   assert( eventhdlr != 0 );

   // add event
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, 0, 0) );

   return SCIP_OKAY;
}
