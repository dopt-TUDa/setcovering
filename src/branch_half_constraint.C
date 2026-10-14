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

/**@file   branch_half_constraint.C
 * @brief  Class for constraint branching
 * @author Marc Pfetsch
 */


#include "branch_half_constraint.h"
#include "setCovering.h"


extern "C" {
#include <scip/cons_logicor.h>
}

#include <sstream>


#define BRANCHRULE_NAME          "BranchHalfConstraint"
#define BRANCHRULE_DESC          "Branching on halves of constraints"
#define BRANCHRULE_PRIORITY      0
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0



BRANCH_HALF_CONSTRAINT::BRANCH_HALF_CONSTRAINT(SETCOVERING* setcovering) :
   ObjBranchrule(setcovering->scip(), BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST)
{
   setcovering_ = setcovering;
}



BRANCH_HALF_CONSTRAINT::~BRANCH_HALF_CONSTRAINT()
{
}




/** Generate a branching
 *
 * Select a branching constraint as follows:
 *
 * We loop through all active logicor constraints and choose the
 * constraint with a maximum number of total number of incident
 * variables. The we sum up its fractional variables until we reach
 * 0.5. Let S be the set of variables obtained this way. Then we
 * generate the following cut:
 *
 * - First  child node: \f$ \sum_{i \in S} x_i = 0 \f$
 * - Second child node: \f$ \sum_{i \in S} x_i \geq 1 \f$
 *
 */
SCIP_DECL_BRANCHEXECLP(BRANCH_HALF_CONSTRAINT::scip_execlp)
{
   assert( branchrule != 0 );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( scip != 0 );
   assert( result != 0 );

   SCIPdebugMsg(scip, "Half-constraint branching ...\n");
   *result = SCIP_DIDNOTRUN;

   // get logicor constraint handler
   SCIP_CONSHDLR* logicor = SCIPfindConshdlr(scip, "logicor");
   assert( logicor != 0 );

   // get number of constraints
   int n = SCIPconshdlrGetNConss(logicor);
   SCIP_CONS** cons = SCIPconshdlrGetConss(logicor);
   assert( cons != 0 );

   // get pointer to array of IIS-containment counter
   int* incident_IISs = setcovering_->incidentVars();

   int max_cnt = -1;
   SCIP_VAR** max_vars = 0;
   int max_n = 0;
   int nactive = 0;
#ifdef SCIP_DEBUG
   SCIP_Real max_sum = 0.0;
   int max_size = 0;
   SCIP_CONS* max_cons = 0;
#endif
   for (int i = 0; i < n; ++i)
   {
      SCIP_CONS* cons_i = cons[i];
      assert( cons_i != 0 );

      // only consider active logicor constraints
      SCIP_Real d = SCIPgetDualsolLogicor(scip, cons[i]);
      if ( SCIPisPositive(scip, d) )
      {
	 // get the variables of the constraint
	 SCIP_VAR** vars_i = SCIPgetVarsLogicor(scip, cons_i);
	 assert( vars_i != 0 );
	 int ni = SCIPgetNVarsLogicor(scip, cons_i);

	 // sum current LP-solution values until we reach 0.5
	 SCIP_Real sum = 0.0;
	 int cnt = 0;
	 int size = 0;
	 int k;
	 for (k = 0; k < ni; ++k)
	 {
	    SCIP_VAR* var = vars_i[k];
	    int index = SCIPvarGetProbindex(var);
	    assert( 0 <= index && index <= setcovering_->nVars() );
	    cnt += incident_IISs[index];

	    sum += SCIPgetSolVal(scip, 0, var);
	    ++size;

	    if ( SCIPisGE(scip, sum, 0.5) )
	       break;
	 }

	 // sum up incident variables for the second half
	 for (; k < ni; ++k)
	 {
	    SCIP_VAR* var = vars_i[k];
	    int index = SCIPvarGetProbindex(var);
	    assert( 0 <= index && index <= setcovering_->nVars() );
	    cnt += incident_IISs[index];
	 }
	 assert( SCIPisPositive(scip, sum) );

	 // if we have at least 2 variables, 2 variables in the first half and the first half leads to a violated cut
	 if ((size >= 2) && (ni-size >= 2) && SCIPisLT(scip, sum, 1.0) )
	 {
	    // if the number of incident IISs is maximum
	    if ( cnt > max_cnt )
	    {
	       max_cnt = cnt;
	       max_vars = vars_i;
	       max_n = size;
#ifdef SCIP_DEBUG
	       max_sum = sum;
	       max_size = ni;
	       max_cons = cons_i;
#endif
	    }
	 }
	 ++nactive;
      }
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Have %d active set covering constraints.\n", nactive);
   SCIPdebugMsg(scip, "max_cnt: %d   max_sum: %f    max_size: %d\n", max_cnt, max_sum, max_size);
   SCIP_CALL( SCIPprintCons(scip, max_cons, 0) );
#endif

   // If no usefull branching was found
   if (max_cnt == -1)
   {
      SCIPdebugMsg(scip, "Could not find branching constraint.\n");
      return SCIP_OKAY;
   }

   // Create children:
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
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, max_n) );
   for (int k = 0; k < max_n; ++k)
   {
      SCIP_VAR* v = max_vars[k];

#ifdef SCIP_DEBUG
      str << " " << SCIPvarGetName(v);
#endif
      // set upper bound of variables to 0 for node 1
      assert( SCIPvarGetLbLocal(v) < 0.5 );
      if ( SCIPvarGetUbLocal(v) > 0.5 )
	 SCIP_CALL( SCIPchgVarUbNode(scip, node1, v, 0.0) );
      vars[cnt++] = v;
   }
#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Branching on: %s with size %d.\n", str.str().c_str(), cnt);
#endif

   // create branching constraint for second node
   SCIP_CALL( SCIPcreateConsLogicor(scip, &c, "cb", cnt, vars, FALSE, TRUE, TRUE, FALSE,
         TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddConsNode(scip, node2, c, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &c) );

   SCIPfreeBufferArray(scip, &vars);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}
