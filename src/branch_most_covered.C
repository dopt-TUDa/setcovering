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

/**@file   branch_most_covered.C
 * @brief  Implementation of class for branching on most covered variable
 * @author Marc Pfetsch
 */


#include "branch_most_covered.h"

#include "setCovering.h"
#include "key_cmp.h"

#include <vector>
#include <algorithm>

#define BRANCHRULE_NAME          "BranchMostCovered"
#define BRANCHRULE_DESC          "Branching on most covered variable"
#define BRANCHRULE_PRIORITY      0
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0



BRANCH_MOST_COVERED::BRANCH_MOST_COVERED(SETCOVERING* setcovering) :
   ObjBranchrule(setcovering->scip(), BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST)
{
   setcovering_ = setcovering;
   nVars_ = setcovering->nVars();
}



BRANCH_MOST_COVERED::~BRANCH_MOST_COVERED()
{
}




//! Generate a branching
SCIP_DECL_BRANCHEXECLP(BRANCH_MOST_COVERED::scip_execlp)
{
   assert( branchrule != 0 );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( scip != 0 );
   assert( result != 0 );

   *result = SCIP_DIDNOTRUN;

   // estimate the number of strong branchings to perform
   int depth = SCIPgetDepth(scip);

   int max_number = nVars_;
   if ( depth > 0 )
   {
      double estimate = double(max_number) / pow(1.5, depth);
      max_number = int(SCIPfloor(scip, estimate));
   }
   SCIPdebugMessage("\nBranching on most covered variable ... [max %d strong branchings]\n", max_number);

   // get branching candidates
   int ncand = 0;
   std::vector<int> indices;
   for (int j = 0; j < nVars_; ++j)
   {
      SCIP_Real x = SCIPgetSolVal(scip, 0, setcovering_->vars(j));
      if ( ! SCIPisIntegral(scip, x) )
      {
	 indices.push_back(j);
	 ++ncand;
      }
   }
   assert( ncand > 0 );
   if (ncand < max_number)
      max_number = ncand;

   // Sort the variables accoring to their number of constraints they are contained in
   int* incident = setcovering_->incidentVars();

   KeyCmp<int> cmp(incident);
   std::sort(indices.begin(), indices.end(), cmp);   // note: keys are not reordered!

   // get current LP objective bound of the local sub problem and global cutoff bound
   SCIP_Real objval = SCIPgetLPObjval(scip);

   // Take some of them and perform strong branching
   SCIP_Real bestscore = -1.0;
   SCIP_VAR* bestvar = setcovering_->vars(indices[0]);    // in doubt: take first variable (w.r.t. to sorting)

   for ( int j = 0; j < max_number; ++j )
   {
      SCIP_VAR* var = setcovering_->vars(indices[j]);

      // take average number of lp iterations as an estimate
      int itlim = 2 * SCIPgetNNodeInitLPIterations(scip)/(SCIPgetNNodeInitLPs(scip) + 1);
      itlim = MAX(itlim, 10);
      SCIP_Real down;
      SCIP_Real up;
      SCIP_Bool lperror;
      SCIP_Bool downinf;
      SCIP_Bool upinf;

      SCIP_CALL( SCIPgetVarStrongbranchFrac(scip, var, itlim, FALSE, &down, &up, 0, 0, &downinf, &upinf, 0, 0, &lperror) );

      if (lperror)
      {
	 SCIPerrorMessage("selectBranchingVariable: unknown error in strong branching!\n");
	 //continue;
	 abort();
      }

      // check if there are infeasible roundings
      if( downinf || upinf )
      {
	 SCIP_Real x = SCIPgetSolVal(scip, 0, var);
	 if( downinf && upinf )
	 {
	    // both roundings are infeasible -> node is infeasible
	    *result = SCIP_CUTOFF;
	    SCIPdebugMessage("Variable %s is infeasible in both directions\n", SCIPvarGetName(var));
	    break;   // terminate initialization loop, because node is infeasible
	 }
	 else if( downinf )
	 {
	    // downwards rounding is infeasible -> change lower bound of variable to upward rounding
	    SCIP_CALL( SCIPchgVarLb(scip, var, SCIPceil(scip, x) ) );
	    SCIPdebugMessage("Variable %s is infeasible in downward branch\n", SCIPvarGetName(var));
	    *result = SCIP_REDUCEDDOM;
	    break;   // terminate initialization loop, because LP was changed
	 }
	 else
	 {
	    // upwards rounding is infeasible -> change upper bound of variable to downward rounding
	    assert(upinf);
	    SCIP_CALL( SCIPchgVarUb(scip, var, SCIPfloor(scip, x)) );
	    SCIPdebugMessage("Variable %s is infeasible in upward branch\n", SCIPvarGetName(var));
	    *result = SCIP_REDUCEDDOM;
	    break;   // terminate initialization loop, because LP was changed
	 }
      }

      down = MAX(down, objval);       // ensure that down >= objval (should be true with exact arithmetic)
      up   = MAX(up, objval);         // ensure that down >= objval (should be true with exact arithmetic)
      SCIP_Real downgain = down - objval;
      SCIP_Real upgain = up - objval;

      // create score
      SCIP_Real score = SCIPgetBranchScore(scip, var, downgain, upgain);
      if( score > bestscore )
      {
	 bestscore = score;
	 bestvar = var;
      }
   }

   // If we could not cut of some child node
   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
   {
      assert( bestvar != 0 );
      SCIPdebugMessage("Choosing variable: %s\n", SCIPvarGetName(bestvar));
      SCIP_CALL( SCIPbranchVar(scip, bestvar, 0, 0, 0) );
      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}
