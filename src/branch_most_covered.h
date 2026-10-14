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

/**@file   branch_most_covered.h
 * @brief  Class for branching on most covered variable
 * @author Marc Pfetsch
 */

#ifndef BRANCH_MOST_COVERED_H
#define BRANCH_MOST_COVERED_H


#include <objscip/objbranchrule.h>


class SETCOVERING;



/**@class BRANCH_MOST_COVERED
 * @brief Class for branching on most covered variable
 *
 * Derived from SCIP class wrapper ObjBranchrule.
 */
class BRANCH_MOST_COVERED : public scip::ObjBranchrule
{
public:

   //! @name Constructor/Destructor
   //@{
   BRANCH_MOST_COVERED(SETCOVERING* setcovering);
   virtual ~BRANCH_MOST_COVERED();
   //@}


   //! @name Virtual functions of SCIP to be redefined
   //@{
   //! Generate a branching (call the corresponding function in MORSE_MASTER)
   virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);
   //@}


private:

   //! @name Data
   //@{
   SETCOVERING* setcovering_;
   int nVars_;
   //@}
};

#endif
