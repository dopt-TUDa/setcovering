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

/**@file   branch_orbital.cpp
 * @brief  Implementation of orbital branching
 * @author Marc Pfetsch
 */

#include "branch_orbital.h"
#include "symmetry.h"
#include "setCovering.h"
#include <objscip/objprobdata.h>

#include <vector>
#include <unistd.h>
#include <sys/times.h>

/* branching rule properties */
#define BRANCHRULE_NAME            "orbital"
#define BRANCHRULE_DESC            "orbital branching"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0



/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   VecVecUInt*           generators;         /**< returns generators */
   SCIP_Bool             computedsym;        /**< whether we have computed symmetries already */
   int                   nfixedzero;         /**< number of variables fixed to 0 */
   int                   nfixedone;          /**< number of variables fixed to 1 */
   SETCOVERING*          setcovering;        /**< pointer to set covering instance */
};



/*
 * Local methods
 */

/** compute symmetry */
static
SCIP_RETCODE computeSymmetry(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< data of branching rule */
   )
{
   if ( branchruledata->computedsym )
      return SCIP_OKAY;
   branchruledata->computedsym = TRUE;
   assert( branchruledata->generators == nullptr );

   int n = branchruledata->setcovering->n();
   int m = branchruledata->setcovering->m();
   int** A = branchruledata->setcovering->A();
   SCIP_Real* obj = branchruledata->setcovering->obj();

   SC_VARFIXING* fixed;
   SCIP_CALL( SCIPallocBufferArray(scip, &fixed, n) );
   SCIP_VAR** vars;
   vars = SCIPgetOrigVars(scip);
   assert( SCIPgetNOrigVars(scip) == n );
   for (int j = 0; j < n; ++j)
   {
      SCIP_VAR* var = vars[j];
      assert( SCIPisEQ(scip, SCIPvarGetObj(var), obj[j]) );
      if ( SCIPvarGetLbGlobal(var) > 0.5 )
         fixed[j] = SCfixedone;
      else if ( SCIPvarGetUbGlobal(var) < 0.5 )
         fixed[j] = SCfixedzero;
      else
         fixed[j] = SCunfixed;
   }

   SCIPinfoMessage(scip, nullptr, "Computing symmetries for orbital branching ...\n");

   struct tms comptimer_beg;
   struct tms comptimer_end;

   // compute symmetries
   long double groupsize;
   branchruledata->generators = new VecVecUInt();
   (void) times(&comptimer_beg);
   SCIP_CALL( computeAutomorphismsSetCovering(scip, m, n, const_cast<const int**>(A), fixed, obj, FALSE, groupsize, *branchruledata->generators) );
   (void) times(&comptimer_end);

   SCIPinfoMessage(scip, nullptr, "Number of generators of Aut(G):\t\t%lu\n", branchruledata->generators->size());
   SCIPinfoMessage(scip, nullptr, "Size of Aut(G):\t\t\t\t%.1Lg\n", groupsize);
   double t = (comptimer_end.tms_utime - comptimer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
   SCIPinfoMessage(scip, nullptr, "Computing symmetry time:\t\t%4.2f\n", t);

   if ( branchruledata->generators->size() == 0 )
   {
      delete branchruledata->generators;
      branchruledata->generators = nullptr;
   }
   SCIPfreeBufferArray(scip, &fixed);

   return SCIP_OKAY;
}


/** Get binary branching variables on the path to root */
static
SCIP_RETCODE computeBranchingVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   n,                  /**< total number of variables */
   SCIP_Bool*            branchedone,        /**< bool array to mark variables branched to 1 */
   int&                  nbranchedone        /**< retruns number of variables branched to 1 */
   )
{
   /* get curent node */
   SCIP_NODE* node = SCIPgetCurrentNode(scip);
   nbranchedone = 0;

   /* follow path to the root (in the root no domains were changed due to branching) */
   while ( SCIPnodeGetDepth(node) != 0 )
   {
      SCIP_BOUNDCHG* boundchg;
      SCIP_DOMCHG* domchg;
      int nboundchgs;

      /* get domain changes of current node */
      domchg = SCIPnodeGetDomchg(node);
      assert( domchg != 0 );

      /* loop through all bound changes */
      nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
      for (int i = 0; i < nboundchgs; ++i)
      {
         /* get bound change info */
         boundchg = SCIPdomchgGetBoundchg(domchg, i);
         assert( boundchg != 0 );

         /* branching decisions have to be in the beginning of the bound change array */
         if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
            break;

         /* get corresponding branching variable */
         SCIP_VAR* branchvar = SCIPboundchgGetVar(boundchg);
         SCIP_Real scalar = 1.0;
         SCIP_Real constant = 0.0;

         SCIP_CALL( SCIPvarGetOrigvarSum(&branchvar, &scalar, &constant) );

         /* we only consider binary variables */
         if ( SCIPvarIsBinary(branchvar) )
         {
            if ( SCIPvarGetLbLocal(branchvar) > 0.5 )
            {
               int branchvaridx = SCIPvarGetIndex(branchvar);
               assert( 0 <= branchvaridx && branchvaridx < n );
               branchedone[branchvaridx] = TRUE;
               ++nbranchedone;
            }
         }
      }

      node = SCIPnodeGetParent(node);
   }

   return SCIP_OKAY;
}

/** largest orbit branching rule */
static
VecVecInt::const_iterator branchLargestOrbit(
   SCIP*                 scip,               /**< SCIP pointer */
   const SCIP_BRANCHRULEDATA* branchruledata,/**< data of branching rule */
   SCIP_VAR**            vars,               /**< variables */
   int                   n,                  /**< total number of variables */
   const VecVecInt&      orbits,             /**< list of orbits among which we choose zero or one to branch on */
   bool                  branchLP,           /**< whether we branch on an LP solution (vs. pseudo solution) */
   VecInt&               fixzero,            /**< list of variables that should be fixed to 0 in both branches */
   VecInt&               fixone,             /**< list of variables that should be fixed to 1 in both branches */
   bool&                 infeasible          /**< whether the current node is infeasible */
   )
{
   VecVecInt::const_iterator orbitToBranch = orbits.end();
   std::size_t largestSize = 0;

   infeasible = false;

   fixzero.clear();
   fixone.clear();

   for (VecVecInt::const_iterator it = orbits.begin(); it != orbits.end(); ++it)
   {
      /* check whether orbit is large enough */
      std::size_t orbitsize = it->size();
      assert( orbitsize >= 2 );

      /* check all variables in the orbit */
      bool havefixedone = false;
      bool havefixedzero = false;
      bool havefractional = false;
      VecInt::const_iterator oit, oend = it->end();
      for (oit = it->begin(); oit != oend; ++oit)
      {
         int idx = *oit;
         assert( 0 <= idx && idx < n );
         SCIP_VAR* var = vars[idx];
         assert( SCIPvarIsBinary(var) );

         /* if variable is fixed to 1 -> can fix all variables to 1, see paper for explanation (-> orbital fixing) */
         if ( SCIPvarGetLbLocal(var) > 0.5 )
            havefixedone = true;

         /* check for zero-fixed variables */
         if ( SCIPvarGetUbLocal(var) < 0.5 )
            havefixedzero = true;

         if ( branchLP && ! SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, var)) )
            havefractional = true;
      }

      /* check consistency */
      if ( havefixedone && havefixedzero )
      {
         infeasible = true;
         break;
      }

      /* fix all variables to 0 if there is one variable fixed to 0 (-> orbital fixing) */
      if ( havefixedzero )
      {
         assert( ! havefixedone );
         assert( oend == it->end() );
         for (oit = it->begin(); oit != oend; ++oit)
         {
            int idx = *oit;
            assert( idx < n );
            /* only variables that are not yet fixed to 0 */
            if ( SCIPvarGetUbLocal(vars[idx]) > 0.5 )
               fixzero.push_back(idx);
         }
      }
      else
      {
         /* fix all variables to 1 if there is one variable fixed to 1 (-> orbital fixing) */
         if ( havefixedone )
         {
            assert( ! havefixedzero );
            assert( oend == it->end() );
            for (oit = it->begin(); oit != oend; ++oit)
            {
               int idx = *oit;
               assert( idx < n );
               /* only variables that are not yet fixed to 1 */
               if ( SCIPvarGetLbLocal(vars[idx]) < 0.5 )
                  fixone.push_back(idx);
            }
         }
         else
         {
            /* do not chose orbit without variables with fractional values */
            if ( ! branchLP || havefractional )
            {
               /* check size */
               if ( orbitsize > largestSize )
               {
                  orbitToBranch = it;
                  largestSize = orbitsize;
               }
            }
         }
      }
   }

   return orbitToBranch;
}


/** perform branching */
static
SCIP_RETCODE branch(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< data of branching rule */
   bool                  branchLP,           /**< whether we branch on an LP solution (vs. pseudo solution) */
   SCIP_RESULT*          result              /**< result of branching */
   )
{
   assert( scip != nullptr );
   assert( branchruledata != nullptr );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->setcovering == nullptr )
   {
      scip::ObjProbData* objprobdata = SCIPgetObjProbData(scip);
      branchruledata->setcovering = static_cast<SETCOVERING*>(objprobdata);
   }

   // make sure that symmetries have been computed */
   SCIP_CALL( computeSymmetry(scip, branchruledata) );

   // if there are no symmetries -> exit
   if ( branchruledata->generators == nullptr )
      return SCIP_OKAY;

   // get branching variables
   int n = branchruledata->setcovering->n();
   int nbranchedone;
   SCIP_Bool* branchedone;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &branchedone, n) );
   SCIP_CALL( computeBranchingVariables(scip, n, branchedone, nbranchedone) );

   SCIPdebugMessage("Number of variables branched to 1: %d.\n", nbranchedone);

   // compute orbits
   VecVecInt orbits;
   SCIP_CALL( computeOrbitsGensFilter(scip, n, branchedone, *branchruledata->generators, orbits) );
   SCIPdebugMessage("Number of orbits: %zu.\n", orbits.size());

#ifndef NDEBUG
   for (size_t i = 0; i < orbits.size(); ++i)
   {
      for (size_t l = 0; l < orbits[i].size(); ++l)
         assert( ! branchedone[orbits[i][l]] );
   }
#endif

   SCIPfreeBufferArray(scip, &branchedone);

   if ( orbits.size() == 0 || SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* choose orbit to branch on */
   VecInt fixzero;
   VecInt fixone;
   bool infeasible;
   SCIP_VAR** vars = SCIPgetOrigVars(scip);
   assert( SCIPgetNOrigVars(scip) == n );
   VecVecInt::const_iterator orbitToBranch = branchLargestOrbit(scip, branchruledata, vars, n, orbits, branchLP, fixzero, fixone, infeasible);

   /* exit if we detected infeasibility */
   if ( infeasible )
   {
      SCIPdebugMessage("Current node detected to be infeasible.\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if ( orbitToBranch == orbits.end() )
   {
      SCIPdebugMessage("Found no orbit to branch on.\n");
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* prepare output of orbit */
   char name[SCIP_MAXSTRLEN];
   char tmp[SCIP_MAXSTRLEN];
   name[0] = '\0';
#endif

   /* compute estimate for 0-branch: */
   SCIP_Real estimate = SCIPgetLocalLowerbound(scip);
   for (VecInt::const_iterator it = orbitToBranch->begin(); it != orbitToBranch->end(); ++it)
   {
      SCIP_VAR* var = vars[*it];
      estimate = MAX(estimate, SCIPcalcChildEstimate(scip, var, 0.0) );

#ifndef NDEBUG
      (void) strncpy(tmp, name, SCIP_MAXSTRLEN-1);
      if ( branchLP )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN-1, "%s <%s> %u [%1f,%1f,%f]", tmp, SCIPvarGetName(var), *it,
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPgetSolVal(scip, NULL, var));
      }
      else
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN-1, "%s <%s> %u [%1f,%1f]", tmp, SCIPvarGetName(var), *it, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
#endif
   }
#ifndef NDEBUG
   SCIPdebugMessage("Branching on orbit:%s\n", name);
#endif

   /* branch on 0 */
   SCIP_Node* node;
   SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, estimate) );
   for (VecInt::const_iterator it = orbitToBranch->begin(); it != orbitToBranch->end(); ++it)
   {
      SCIP_VAR* var = vars[*it];
      assert( SCIPvarGetUbLocal(var) > 0.5 );
      SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
   }

   /* fix variables to 0 */
   for (VecInt::const_iterator fit = fixzero.begin(); fit != fixzero.end(); ++fit)
   {
      SCIP_VAR* var = vars[*fit];

      assert( SCIPvarGetLbLocal(var) < 0.5 );
      if ( SCIPvarGetUbLocal(var) > 0.5 )
      {
         SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
         ++branchruledata->nfixedzero;
      }
   }

   /* fix variables to 1 */
   for (VecInt::const_iterator fit = fixone.begin(); fit != fixone.end(); ++fit)
   {
      SCIP_VAR* var = vars[*fit];

      assert( SCIPvarGetUbLocal(var) > 0.5 );
      if ( SCIPvarGetLbLocal(var) < 0.5 )
      {
         SCIP_CALL( SCIPchgVarLbNode(scip, node, var, 1.0) );
         ++branchruledata->nfixedone;
      }
   }

   /* take first unfixed variable in orbit as variable to branch to 1 */
   int idx = *(orbitToBranch->begin());
   assert( idx < n );
   SCIP_VAR* branchvar = vars[idx];
   assert( SCIPvarGetLbLocal(branchvar) < 0.5 );    /* variable should not be fixed, by the choice above */

   SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, SCIPcalcChildEstimate(scip, branchvar, 1.0)) );
   SCIP_CALL( SCIPchgVarLbNode(scip, node, branchvar, 1.0) );

   /* fix variables to 0 */
   for (VecInt::const_iterator fit = fixzero.begin(); fit != fixzero.end(); ++fit)
   {
      SCIP_VAR* var = vars[*fit];

      assert( SCIPvarGetLbLocal(var) < 0.5 );
      if ( SCIPvarGetUbLocal(var) > 0.5 )
      {
         SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
      }
   }

   /* fix variables to 1 */
   for (VecInt::const_iterator fit = fixone.begin(); fit != fixone.end(); ++fit)
   {
      SCIP_VAR* var = vars[*fit];

      assert( SCIPvarGetUbLocal(var) > 0.5 );
      if ( SCIPvarGetLbLocal(var) < 0.5 )
      {
         SCIP_CALL( SCIPchgVarLbNode(scip, node, var, 1.0) );
      }
   }

   SCIPdebugMessage("Variables fixed to 0/1 in both branches: [%lu, %lu].\n", (unsigned long) fixzero.size(), (unsigned long) fixone.size());

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/*
 * Callback methods of branching rule
 */


/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeOrbital)
{  /*lint --e{715,818}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   assert( scip != 0 );
   assert( branchrule != 0 );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != 0 );

   if ( branchruledata->generators != nullptr )
      delete branchruledata->generators;
   SCIPfreeBlockMemory(scip, &branchruledata);

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitOrbital)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   assert( scip != 0 );
   assert( branchrule != 0 );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != 0 );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "\nNumber of variables fixed to 0 by orbital branching:\t%6d\n", branchruledata->nfixedzero);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Number of variables fixed to 1 by orbital branching:\t%6d\n\n", branchruledata->nfixedone);

   return SCIP_OKAY;
}


/** Initialization method of branching rule (called after problem was transformed) */
SCIP_DECL_BRANCHINIT(branchInitOrbital)
{
   assert( branchrule != 0 );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( scip != 0 );

   SCIPdebugMessage("Init of orbital branching ...\n");

   // We assume that if the priority is positive, then orbital branching should be used and we have to turn off symmetry
   // handling in SCIP.
   int priority;
   SCIP_CALL( SCIPgetIntParam(scip, "branching/orbital/priority", &priority) );
   if ( priority > 0 )
   {
      SCIPwarningMessage(scip, "Turning off symmetry, because orbital branching might be active.\n");
      SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 0) );
   }

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpOrbital)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   assert( branchrule != 0 );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( scip != 0 );
   assert( result != 0 );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != 0 );

   SCIPdebugMessage("Orbital branching (node: %lld, depth: %d) ...\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), SCIPgetDepth(scip));

   SCIP_CALL( branch(scip, branchruledata, true, result) );

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrbital)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   assert( branchrule != 0 );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( scip != 0 );
   assert( result != 0 );

   SCIPdebugMessage("Orbital branching (node: %lld, depth: %d) ...\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), SCIPgetDepth(scip));

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != 0 );

   SCIP_CALL( branch(scip, branchruledata, false, result) );

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** creates the xyz branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleOrbital(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata = nullptr;

   /* create branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );
   branchruledata->computedsym = FALSE;
   branchruledata->generators = nullptr;
   branchruledata->nfixedzero = 0;
   branchruledata->nfixedone = 0;
   branchruledata->setcovering = nullptr;

   /* include branching rule */
   SCIP_BRANCHRULE* branchrule = nullptr;
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   /* set callbacks */
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitOrbital) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeOrbital) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitOrbital) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpOrbital) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsOrbital) );

   return SCIP_OKAY;
}
