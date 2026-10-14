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

/**@file   setCovering.h
 * @brief  Top class for handling optimization data
 * @author Marc Pfetsch
 */


#ifndef SETCOVERING_H
#define SETCOVERING_H

#include "Bitset.h"
#include "type_presolve.h"
#include <objscip/objprobdata.h>

/** @mainpage Solving Set Covering Problems
 *  @author  Marc Pfetsch
 *  @version 0.1
 *  @date 8/2006
 *
 *  The program @c setCovering is a simple implementation of a branch-and-cut method to
 *  solve set covering problems. We use the SCIP framework for this process.
 *
 */


/**@class SETCOVERING
 * @brief Top class for handling optimization data
 */
class SETCOVERING : public scip::ObjProbData
{
public:

   //! Constructor setting up everything and reading input file
   SETCOVERING(const char* problemName, const char* fileName, SCIP* scip);

   //! Destructor
   ~SETCOVERING();


   //! @name General structural functions
   //@{
   //! Read solution file
   void readSolution(const char* filename);

   //! Write solution file
   void writeSolution(const char* filename);

   //! include plugins
   SCIP_RETCODE includePlugins();

   //! Initialize everything needed before optimization, e.g., variables, constraints, pools, etc.
   SCIP_RETCODE initializeOptimization(const char* paramfile);

   //! Cleanup stuff left over from optimization and output statistics
   SCIP_RETCODE terminateOptimization();

   //! Perform symmetry presolving
   SCIP_RETCODE performSymmetryPresolving();
   //@}


   //! @name Access to important SCIP data
   //@{
   //! Return SCIP pointer
   SCIP* scip() { return scip_; }

   //! Returns the array of variables
   SCIP_VAR** vars() { return vars_; }

   //! Returns the i-th variable
   SCIP_VAR* vars(int i) { return vars_[i]; }

   //! Returns the i-th original variable
   SCIP_VAR* origVars(int i) { return origVars_[i]; }

   //! Returns the number of variables
   int nVars() const { return n_; }

   //! Returns pointer to array that stores for each variable the number constraints it is contained in
   int* incidentVars() { return incidentVars_; }

   //! return number of columns
   int n() const { return n_; };

   //! return number of rows
   int m() const { return m_; };

   //! return matrix
   int** A() { return A_; };

   //! return objective
   SCIP_Real* obj() { return obj_; };
   //@}



   //! @name Handling of the primal solutions found
   //@{
   //! Store initial cover (the best cover is stored)
   bool initialCover(const Bitset& S, int size);

   //! Store solution if it is better than the previously found solutions
   bool betterCover(SCIP_Real value, const Bitset& S, int size, SCIP_HEUR* heur);

   //! SCIP found a better solution
   void foundBetterSolution();

   //! Returns the objective function value of the best primal solution found so far.
   SCIP_Real bestCoverValue() const { return bestCoverValue_; }

   //! Returns the best primal solution found so far.
   const Bitset* bestCover() const { return bestCover_; }

   //! Returns the size of the best primal solution found so far.
   int bestCoverSize() const { return bestCoverSize_; }

   //! Tell SCIP about solution
   SCIP_RETCODE submitSolution(SCIP_HEUR* heur);
   //@}


   // ----------------------- functions --------------------------
   //! @name Other functions
   //@{
   //! Read problem file in own format
   SCIP_RETCODE readProblem(const char* filename);

   //! Read problem file in LP format
   SCIP_RETCODE readProblemLP(const char* filename);
   //@}

   //! @name SCIP virtual functions
   //@{
   //! Create data for the transformed problem
   virtual SCIP_RETCODE scip_trans(
      SCIP*         scip,               /**< SCIP data structure */
      ObjProbData** objprobdata,        /**< pointer to store the transformed problem data object */
      SCIP_Bool*    deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
      );

   //! Delete data of the transformed problem
   virtual SCIP_RETCODE scip_deltrans(SCIP* scip);
   //@}

   // ------------------------ data -------------------------------
   //! pointer to SCIP
   SCIP *scip_;

   //! @name Variables
   //@{
   SCIP_VAR** vars_;         //!< SCIP variables (transformed)
   SCIP_VAR** origVars_;     //!< original variables
   //@}


   //! @name best primal solution found so far
   //@{
   Bitset* bestCover_;               //!< best primal solution
   int bestCoverSize_;              //!< size of best primal solution
   SCIP_Real bestCoverValue_;       //!< objective value of best primal solution
   //@}

   //! @name Problem Data
   //@{
   std::string problemName_;   //!< problemname
   int norig_;            // original number of variables
   int morig_;            // original number of rows
   int** Aorig_;          // original constraint matrix
   SCIP_Real* objorig_;   // original objective function

   int n_;                // presolved number of variables
   int m_;                // presolved number of rows
   int** A_;              // presolved constraint matrix
   SCIP_Real* obj_;       // original objective function
   SC_VARFIXING* fixed_;  // fixings applied to original matrix in presolving
   int nfixed_;           // number of fixed variables

   int* incidentVars_;    // incident constraints

   SCIP_Bool presolving_; // perform presolving?
   SCIP_Bool preaddition_;// perform addition presolving?
   SCIP_Bool premerge_;   // perform merge presolving?
   SCIP_Bool prefixation_;//!< Perform fixation presolving?

   SCIP_Real maxfactor_;  // factor for early exit
   int maxpreadd_;        // number of rows to add for each leader/follower combination (-1 = unbounded)
   SCIP_Bool output_optimal_cover_;
   SCIP_Bool writeBestSolution_;
   //@}
};

#endif
