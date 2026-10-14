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

/**@file   branch_intersection.h
 * @brief  Class for constraint branching
 * @author Marc Pfetsch
 */

#ifndef BRANCH_INTERSECTION_H
#define BRANCH_INTERSECTION_H


#include <objscip/objbranchrule.h>


class SETCOVERING;



/**@class BRANCH_INTERSECTION
 * @brief Class for constraint branching
 *
 * Derived from SCIP class wrapper ObjBranchrule.
 */
class BRANCH_INTERSECTION : public scip::ObjBranchrule
{
public:

   //! @name Constructor/Destructor
   //@{
   BRANCH_INTERSECTION(SETCOVERING* setcovering);
   virtual ~BRANCH_INTERSECTION();
   //@}

   //! @name Virtual functions of SCIP to be redefined
   //@{
   //! Generate a branching (call the corresponding function in MORSE_MASTER)
   virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);
   //@}

private:

   //! @name Data
   //@{
   SETCOVERING* setcovering_;    //!< pointer to main class
   //@}
};

#endif
