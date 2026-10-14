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

/**@file   branch_intersection.C
 * @brief  Class for constraint branching
 * @author Marc Pfetsch
 */


#include "branch_intersection.h"

#include "setCovering.h"


extern "C" {
#include <scip/cons_logicor.h>
}

#include <sstream>



#define BRANCHRULE_NAME          "BranchIntersection"
#define BRANCHRULE_DESC          "Branching on intersection"
#define BRANCHRULE_PRIORITY      0
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0



BRANCH_INTERSECTION::BRANCH_INTERSECTION(SETCOVERING* setcovering) :
   ObjBranchrule(setcovering->scip(), BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST)
{
   setcovering_ = setcovering;
}



BRANCH_INTERSECTION::~BRANCH_INTERSECTION()
{
}



/**@brief Generate a branching on a constraint
 *
 * The branching constraint is selected as follows:
 *
 * For all pairs of active IIS-constraints compute their intersection S.
 * Choose the pair that maximizes the sum of the variables in the set S.
 * Generate the following two branching constraints:
 *
 * First  child node: $\\sum_{i \\in S} x_i = 0$.
 *
 * Second child node: $\\sum_{i \\in S} x_i \\geq 1$.
 *
 * This is kind of a Ryan/Foster branching rule.
 */
SCIP_DECL_BRANCHEXECLP(BRANCH_INTERSECTION::scip_execlp)
{
   assert( branchrule != 0 );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( scip != 0 );
   assert( result != 0 );

   SCIPdebugMessage("Branching on the intersection of set covering constraints ...\n");
   *result = SCIP_DIDNOTRUN;

   // get logicor constraint handler
   SCIP_CONSHDLR* logicor = SCIPfindConshdlr(scip, "logicor");
   assert( logicor != 0 );

   // get number of constraints
   int n = SCIPconshdlrGetNConss(logicor);
   SCIP_CONS** cons = SCIPconshdlrGetConss(logicor);
   assert( cons != 0 );

   SCIPdebugMessage("Have %d active set covering constraints.\n", n);

   // We compute the intersections of every two set covering constraints and
   // search for the maximal fractional value of sum over this intersection.
   SCIP_Real max_sum = 0.0;
   SCIP_VAR** max_vars_i = 0;
   SCIP_VAR** max_vars_j = 0;
   int max_ni = 0;
   int max_nj = 0;
   int max_size = 0;

   // loop through all pairs of logicor constraints
   for (int i = 0; i < n; ++i)
   {
      SCIP_CONS* cons_i = cons[i];
      assert( cons_i != 0 );
      SCIP_VAR** vars_i = SCIPgetVarsLogicor(scip, cons_i);
      assert( vars_i != 0 );
      int ni = SCIPgetNVarsLogicor(scip, cons_i);

      for (int j = i+1; j < n; ++j)
      {
	 SCIP_CONS* cons_j = cons[j];
	 assert( cons_j != 0 );
	 SCIP_VAR** vars_j = SCIPgetVarsLogicor(scip, cons_j);
	 assert( vars_j != 0 );
	 int nj = SCIPgetNVarsLogicor(scip, cons_j);

	 // compute sum in intersection
	 SCIP_Real sum = 0.0;
	 int size = 0;
	 for (int k = 0; k < ni; ++k)
	 {
	    for (int l = 0; l < nj; ++l)
	    {
	       if (vars_i[k] == vars_j[l])
	       {
		  sum += SCIPgetSolVal(scip, 0, vars_i[k]);
		  ++size;
		  break;
	       }
	    }
	 }

	 // take only intersections that have a sum in (0,1) - the other cases are not interesting for the branching
	 // (when sum >= 1  one branch is automatically satisfied)
	 if ( SCIPisSumPositive(scip, sum) && SCIPisSumLT(scip, sum, 1.0) )
	 {
	    SCIP_Real dx = sum;
	    if (sum > 0.5)
	       dx = 1.0 - dx;

	    if (dx > max_sum)
	    {
	       max_sum = dx;
	       max_vars_i = vars_i;
	       max_vars_j = vars_j;
	       max_ni = ni;
	       max_nj = nj;
	       max_size = size;
	       if ( fabs(0.5 - max_sum) < 0.1 )
		  break;
	    }
	 }
      }
      if ( fabs(0.5 - max_sum) < 0.1 )
	 break;
   }

   // If no usefull branching was found
   if (max_size == 0)
   {
      SCIPdebugMessage("Could not find branching constraint.\n");
      return SCIP_OKAY;
   }

   // Create childs:
   SCIP_Node* node1;     // Node 1   x(S) == 0
   SCIP_Node* node2;     // Node 2   x(S) >= 1
   SCIP_CONS* c;         // branching constraint for node 2

   SCIP_CALL( SCIPcreateChild(scip, &node1, 0.0, setcovering_->bestCoverValue()) );
   SCIP_CALL( SCIPcreateChild(scip, &node2, 0.0, setcovering_->bestCoverValue()) );

#ifdef SCIP_DEBUG
   std::stringstream str;
#endif
   int cnt = 0;
   SCIP_VAR** vars;
   int ntotalvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, ntotalvars) );
   for (int k = 0; k < max_ni; ++k)
   {
      SCIP_VAR* v = max_vars_i[k];
      int l;
      for (l = 0; l < max_nj; ++l)
      {
	 if (v == max_vars_j[l])
	 {
#ifdef SCIP_DEBUG
	    str << " " << SCIPvarGetName(v);
#endif
	    // set upper bound of variables to 0 for node 1
	    assert( SCIPvarGetLbLocal(v) < 0.5 );
	    if ( SCIPvarGetUbLocal(v) > 0.5 )
	       SCIP_CALL( SCIPchgVarUbNode(scip, node1, v, 0.0) );
	    vars[cnt++] = v;
	    break;
	 }
      }
   }
   assert( cnt < ntotalvars );
#ifdef SCIP_DEBUG
   SCIPdebugMessage("Branching on: %s with size %d.\n", str.str().c_str(), cnt);
#endif

   SCIP_CALL( SCIPcreateConsLogicor(scip, &c, "cb", cnt, vars, FALSE, TRUE, TRUE, FALSE,
         TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddConsNode(scip, node2, c, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &c) );

   SCIPfreeBufferArray(scip, &vars);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}
