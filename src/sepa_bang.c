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
#pragma ident "@(#) $Id: sepa_bang.c,v 1.2 2006-11-22 11:47:27 bzfpfets Exp $"

/**@file   sepa_bang.c
 * @brief  Balas/Ng cut separator
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
/*#include <scip/cons_setppc.h>*/
#include <scip/cons_logicor.h>


#include "sepa_bang.h"


#define SEPA_NAME              "bang"
#define SEPA_DESC              "Balas/Ng cut separator"
#define SEPA_PRIORITY                -1
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_DELAY                 TRUE
#define SEPA_USESSUBSCIP          FALSE

#define DEFAULT_MAXSIZE             100



/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   int max_size;         /**< maximum number of (randomly chosen) set covering constraints to consider */
   SCIP_Bool strengthen; /**< whether cuts should be strengthend */
};



/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeBaNg)
{
   SCIP_SEPADATA* sepadata;

   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}





/**@brief Try to strengthen a violated cut by considering other rows
 * @param scip
 * @param consIndices	indices of constraints already used
 * @param depth		number of constraints used before
 * @param sum		lhs of previous cut
 * @param cons		constraints
 * @param ncons		number of constraints
 * @param vars		sparse array of variables in cut
 * @param vals		values for variables in cut
 * @param index		indices of variables in cuts
 * @param nindices	number of variables in cut
 * @param value		dense counter for occurence of col. in cut
 * @param nvars		number of variables
 */
static
int strengthenViolatedCut(SCIP* scip, SCIP_SEPA* sepa, int* consIndices, int depth, SCIP_Real sum, SCIP_CONS** cons, int ncons,
   SCIP_VAR** vars, SCIP_Real* vals, int* index, int nindices, int* value, int nvars)
{
   int l, s;
   SCIP_VAR** vars_s;
   int ns;
   int nindices_s;
   SCIP_VAR* var;
   int varindex;
   SCIP_Real sumNew;
   int strengthened;
   SCIP_ROW* row;

   /* Bound on recursion depth */
   /*
   if (depth > 10)
      return 0;
   */

   strengthened = 0;  /* counts the number of generated strengthened cuts */

   /* loop over all constraints */
   for ( s = 0; s < ncons; ++s )
   {
      /* check whether constraints was used before */
      for (l = 0; l < depth; ++l )
      {
	 if (s == consIndices[l])
	    break;
      }
      if (l < depth)
	 continue;

      /* now get constraint description */
      vars_s = SCIPgetVarsLogicor(scip, cons[s]);
      ns = SCIPgetNVarsLogicor(scip, cons[s]);
      assert( vars_s != 0 );
      assert( ns > 0 );

      /* enter coefficients and indices */
      nindices_s = nindices;
      for( l = 0; l < ns; ++l )
      {
	 var = vars_s[l];
	 varindex = SCIPvarGetProbindex(var);
	 assert( 0 <= varindex && varindex < nvars );
	 assert( 0 <= value[varindex] && value[varindex] <= depth);
	 value[varindex] += 1;
	 if (value[varindex] == 1)
	 {
	    index[nindices_s] = varindex;
	    vars[nindices_s++] = var;
	    assert(nindices_s <= nvars);
	 }
      }

      /* Check whether the cut is violated */
      sumNew = 0.0;
      for( l = 0; l < nindices_s; ++l )
      {
	 assert( 0 <= index[l] && index[l] < nvars );
	 assert( 0 < value[index[l]] && value[index[l]] <= (depth+1) );
	 if( value[index[l]] == (depth+1) )
	    sumNew += 2.0 * SCIPgetSolVal(scip, NULL, vars[l]);
	 else
	    sumNew += SCIPgetSolVal(scip, NULL, vars[l]);
      }

      /* if the cut is violated recurse */
      if ( SCIPisLT(scip, sumNew, sum) )
      {
         SCIP_Bool infeasible;
	 int nGen;

	 consIndices[depth] = s;
	 nGen = strengthenViolatedCut(scip, sepa, consIndices, depth+1, sumNew, cons, ncons, vars, vals, index,
				      nindices_s, value, nvars);
	 if ( nGen == 0 )
	 {
	    /* could not strengthen cut - generate own cut */
	    for( l = 0; l < nindices_s; ++l )
	    {
	       if( value[index[l]] == (depth+1) )
		  vals[l] = 2.0;
	       else
		  vals[l] = 1.0;
	    }

	    SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "", 2.0, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
	    SCIP_CALL( SCIPaddVarsToRow(scip, row, nindices_s, vars, vals) );
	    SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            assert( ! infeasible );
	    SCIP_CALL( SCIPaddPoolCut(scip, row) );

	    assert( SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) );

	    /* SCIP_CALL( SCIPprintRow(scip, row, NULL) ); */
	    SCIP_CALL( SCIPreleaseRow(scip, &row) );
	    ++strengthened;
	 }
	 else
	    strengthened += nGen;
      }

      /* correct values array for s */
      for( l = 0; l < ns; ++l )
      {
	 var = vars_s[l];
	 varindex = SCIPvarGetProbindex(var);
	 assert( 0 <= varindex && varindex < nvars );
	 value[varindex] -= 1;
	 assert( 0 <= value[varindex] && value[varindex] <= depth );
      }
   }
   return strengthened;
}








/** separates Balas/Ng cuts and adds them to SCIP */
static
SCIP_DECL_SEPAEXECLP(sepaExecBaNg)
{
   int i, l;     /* iterator variables */
   int ncuts;    /* number of generated cuts */
   int nstr;     /* number of strengthed cuts */

   int nvars;                /* total number of active variables */
   SCIP_CONSHDLR* logicor;   /* logicor constraint handler */
   SCIP_CONS** logicor_cons; /* constraints of logicor cons. hdlr. */
   int logicor_ncons;        /* number of logicor constraints altogether */

   SCIP_CONS** cons;     /* set of choosen set covering constraints */
   int ncons;            /* number of constraints in cons */

   SCIP_VAR** vars_i;    /* variables of the current constraints */
   SCIP_VAR** vars_j;
   SCIP_VAR** vars_k;
   int ni, nj, nk;       /* number of variables in current constraints */
   int nindices_i;       /* number of indices in sparse array of active columns at i-level */
   int nindices_j;       /* number of indices in sparse array of active columns at j-level */
   int nindices_k;       /* number of indices in sparse array of active columns at k-level */

   int* index;           /* sparse array of column indices that occur in one of the three constraints */
   SCIP_VAR** vars;      /* sparse array of variables that occur in one of the three constraints */
   SCIP_Real* vals;      /* sparse array of coefficients for cut */
   int* value;           /* sum of the occurences of column indices among the three constraints (non-sparse) */
   int* consIndices;     /* indices for recursive strengthening */
   SCIP_Real* valsNew;   /* new values for strengthened cuts */

   SCIP_VAR* var;        /* helper variables */
   int varindex;
   SCIP_Real sum;
   SCIP_ROW* row;
   SCIP_SEPADATA* sepadata;
   int max_size;
   SCIP_Bool strengthen;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;
   ncuts = 0;
   nstr = 0;

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );
   max_size = sepadata->max_size;       /* maximal number of set covering constraints to consider */
   strengthen = sepadata->strengthen;   /* whether cuts should be strengthend */

   /* get logicor constraint handler */
   logicor = SCIPfindConshdlr(scip, "logicor");
   assert( logicor != 0 );

   /* get number of constraints */
   logicor_ncons = SCIPconshdlrGetNConss(logicor);
   if ( logicor_ncons == 0 )
      return SCIP_OKAY;
   logicor_cons = SCIPconshdlrGetConss(logicor);
   assert( logicor_cons != 0 );

   /* collect active IIS constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &cons, logicor_ncons) );      /* get space for constraints */
   ncons = 0;
   for ( l = 0; l < logicor_ncons; ++l)
   {
      /*if ( SCIPgetTypeSetppc(scip, setppc_cons[l]) == SCIP_SETPPCTYPE_COVERING &&
	!SCIPconsIsModifiable(setppc_cons[l]) ) */
      if ( !SCIPconsIsModifiable(logicor_cons[l]) )
      {
	 SCIP_Real d = SCIPgetDualsolLogicor(scip, logicor_cons[l]);
	 if ( SCIPisPositive(scip, d) )
	    cons[ncons++] = logicor_cons[l];
      }
   }
   SCIPdebugMessage("BaNg separator: have %d active IIS constraints.\n", ncons);
   fflush(stdout);

   /* if we do not have enough set covering constraints */
   if (ncons <= 2)
   {
      SCIPfreeBufferArray(scip, &cons);
      return SCIP_OKAY;
   }

   /* If # of set covering constraints is too large */
   if ( ncons >= max_size )
   {
      /* take random subset of size max_size. Warning this set is not exactly random ... */
      while ( ncons > max_size )
      {
	 int choose;
	 choose = (int) (((double) ncons) * drand48());
	 assert( choose < ncons );
	 if ( choose < ncons-1 )
	    cons[choose] = cons[ncons-1];
	 --ncons;
      }
      SCIPdebugMessage("BaNg separator: reduced to %d constraints.\n", ncons);
      fflush(stdout);
   }

   /* compute cuts */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &index, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &value, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   BMSclearMemoryArray(value, nvars);

   /* space for strengthening cut */
   SCIP_CALL( SCIPallocBufferArray(scip, &consIndices, ncons) );
   SCIP_CALL( SCIPallocBufferArray(scip, &valsNew, nvars) );

   /* three-fold loop */
   *result = SCIP_DIDNOTFIND;
   for( i = 0; i < ncons; ++i )
   {
      int j;

      vars_i = SCIPgetVarsLogicor(scip, cons[i]);
      ni = SCIPgetNVarsLogicor(scip, cons[i]);
      assert( vars_i != 0 );
      assert( ni > 0 );

      /* enter coefficients and indices */
      nindices_i = 0;
      for( l = 0; l < ni; ++l )
      {
	 var = vars_i[l];
	 varindex = SCIPvarGetProbindex(var);
	 assert( 0 <= varindex && varindex < nvars );
	 assert( value[varindex] == 0 );
	 value[varindex] = 1;
	 index[nindices_i] = varindex;
	 vars[nindices_i++] = var;
      }

      /* second level */
      for( j = i+1; j < ncons; ++j )
      {
	 int k;

	 vars_j = SCIPgetVarsLogicor(scip, cons[j]);
	 nj = SCIPgetNVarsLogicor(scip, cons[j]);
	 assert( vars_j != 0 );
	 assert( nj > 0 );

	 /* enter coefficients and indices */
	 nindices_j = nindices_i;
	 for( l = 0; l < nj; ++l )
	 {
	    var = vars_j[l];
	    varindex = SCIPvarGetProbindex(var);
	    assert( 0 <= varindex && varindex < nvars );
	    assert( value[varindex] == 0 || value[varindex] == 1 );
	    value[varindex] += 1;
	    if (value[varindex] == 1)
	    {
	       index[nindices_j] = varindex;
	       vars[nindices_j++] = var;
	    }
	 }

	 /* third level */
	 for ( k = j+1; k < ncons; ++k )
	 {
	    vars_k = SCIPgetVarsLogicor(scip, cons[k]);
	    nk = SCIPgetNVarsLogicor(scip, cons[k]);
	    assert( vars_k != 0 );
	    assert( nk > 0 );

	    /* enter coefficients and indices */
	    nindices_k = nindices_j;
	    for( l = 0; l < nk; ++l )
	    {
	       var = vars_k[l];
	       varindex = SCIPvarGetProbindex(var);
	       assert( 0 <= varindex && varindex < nvars );
	       assert( 0 <= value[varindex] && value[varindex] <= 2 );
	       value[varindex] += 1;
	       if (value[varindex] == 1)
	       {
		  index[nindices_k] = varindex;
		  vars[nindices_k++] = var;
	       }
	    }

	    /* Check whether the cut is violated */
	    sum = 0.0;
	    for( l = 0; l < nindices_k; ++l )
	    {
	       assert( 0 <= index[l] && index[l] < nvars );
	       assert( 0 < value[index[l]] && value[index[l]] <= 3 );
	       if( value[index[l]] == 3 )
	       {
		  sum += 2.0 * SCIPgetSolVal(scip, NULL, vars[l]);
		  vals[l] = 2.0;
	       }
	       else
	       {
		  sum += SCIPgetSolVal(scip, NULL, vars[l]);
		  vals[l] = 1.0;
	       }
	    }

	    /* if the cut is violated, strengthen cut or generate cut */
	    if ( SCIPisFeasLT(scip, sum, 2.0) )
	    {
               SCIP_Bool infeasible;

	       /* if strengthening is turned on */
	       if ( strengthen )
	       {
		  int nGen;

		  consIndices[0] = i;
		  consIndices[1] = j;
		  consIndices[2] = k;
		  nGen = strengthenViolatedCut(scip, sepa, consIndices, 3, sum, cons, ncons, vars, valsNew, index, nindices_k, value, nvars);
		  if ( nGen == 0 )
		  {
		     /* generate old cut */
		     SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "", 2.0, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
		     SCIP_CALL( SCIPaddVarsToRow(scip, row, nindices_k, vars, vals) );
		     SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
                     assert( ! infeasible );
		     SCIP_CALL( SCIPaddPoolCut(scip, row) );

		     assert( SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) );

		     SCIP_CALL( SCIPreleaseRow(scip, &row) );
		     *result = SCIP_SEPARATED;
		  }
		  ncuts += nGen;
		  nstr += nGen;
	       }
	       else
	       {
		  /* generate cut */
		  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "", 2.0, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
		  SCIP_CALL( SCIPaddVarsToRow(scip, row, nindices_k, vars, vals) );
		  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
                  assert( ! infeasible );
		  SCIP_CALL( SCIPaddPoolCut(scip, row) );
		  /* SCIP_CALL( SCIPprintRow(scip, row, NULL) ); */

		  assert( SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) );

		  SCIP_CALL( SCIPreleaseRow(scip, &row) );
		  *result = SCIP_SEPARATED;
		  ++ncuts;
	       }
	    }

	    /* correct values array for k */
	    for( l = 0; l < nk; ++l )
	    {
	       var = vars_k[l];
	       varindex = SCIPvarGetProbindex(var);
	       assert( 0 <= varindex && varindex < nvars );
	       value[varindex] -= 1;
	       assert( 0 <= value[varindex] && value[varindex] <= 2 );
	    }
	 }
	 /* correct values array for j */
	 for( l = 0; l < nj; ++l )
	 {
	    var = vars_j[l];
	    varindex = SCIPvarGetProbindex(var);
	    assert( 0 <= varindex && varindex < nvars );
	    value[varindex] -= 1;
	    assert( value[varindex] == 0 || value[varindex] == 1 );
	 }
      }
      /* correct values array for i */
      for( l = 0; l < ni; ++l )
      {
	 var = vars_i[l];
	 varindex = SCIPvarGetProbindex(var);
	 assert( 0 <= varindex && varindex < nvars );
	 value[varindex] -= 1;
	 assert( value[varindex] == 0 );
      }
   }
   /* free space for strengthening cuts */
   SCIPfreeBufferArray(scip, &valsNew);
   SCIPfreeBufferArray(scip, &consIndices);

   /* free space */
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &value);
   SCIPfreeBufferArray(scip, &index);
   SCIPfreeBufferArray(scip, &cons);

   if ( strengthen )
   {
      SCIPdebugMessage("BaNg separator: separated %d cuts [%d strengthend].\n", ncuts, nstr);
   }
   else
   {
      SCIPdebugMessage("BaNg separator: separated %d cuts.\n", ncuts);
   }
   fflush(stdout);

   return SCIP_OKAY;
}





/*
 * separator specific interface methods
 */

/** creates the Balas/Ng cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaBaNg(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_USESSUBSCIP, SEPA_DELAY,
         NULL, sepaFreeBaNg, NULL, NULL, NULL, NULL, sepaExecBaNg, NULL, sepadata) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "separating/bang/maxsize",
         "maximal number of considered set covering constraint",
         &sepadata->max_size, FALSE, DEFAULT_MAXSIZE, -1, INT_MAX, NULL, NULL) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/bang/strengthen",
         "whether violated cuts should be strengthend",
         &sepadata->strengthen, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
