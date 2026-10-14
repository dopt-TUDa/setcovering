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

/**@file   main.C
 * @brief  main file
 * @author Marc Pfetsch
 */


#include "setCovering.h"
#include "macro.h"
#include "get_problem_name.h"
#include "changesoltime.h"
#include <scip/scipdefplugins.h>

#include <iostream>
#include <unistd.h>
#include <sys/times.h>

//! macro to check for a SCIP error and possibly exit
#define SCIP_CALL_ERROR(x)   do                                                                               \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPprintError(_restat_);                                                        \
                             return -1;                                                                       \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

//! read comand line arguments
static
SCIP_Bool readArguments(
   SCIP*                 scip,               //!< SCIP pointer
   int                   argc,               //!< number of shell parameters
   const char**          argv,               //!< array with shell parameters
   const char**          filename,           //!< file name for network from arguments
   const char**          settingsname,       //!< name of settings file
   const char**          solutionfile,       //!< name of solution file
   SCIP_Real&            timelimit,          //!< time limit read from arguments
   SCIP_Real&            memlimit,           //!< memory limit read from arguments
   SCIP_Longint&         nodelimit,          //!< node limit read from arguments
   int&                  dispfreq,           //!< display frequency
   int&                  verblvl,            //!< display verbosity
   SCIP_Real&            cutoffvalue         //!< cutoff value
   )
{
   char usage[5000];
   int i;

   assert( scip != nullptr );
   assert( argv != nullptr );
   assert( filename != nullptr );
   assert( settingsname != nullptr );
   assert( solutionfile != nullptr );

   // init usage text
   (void) snprintf(usage, (unsigned long)5000,
      "usage: <filename> [-s <setting file>] [-l <solution file>] [-t <time limit>] [-m <mem limit>] [-n <node limit>] [-d <display frequency>] [-v <display verbosity>] [-c cutoff]");

   assert( strlen(usage) < 5000 );

   // init arguments
   *filename  = nullptr;
   *settingsname = nullptr;
   *solutionfile = nullptr;
   timelimit = SCIPinfinity(scip);
   memlimit = SCIPinfinity(scip);
   nodelimit = SCIP_LONGINT_MAX;
   dispfreq = -1;
   verblvl = -1;
   cutoffvalue = -SCIPinfinity(scip);

   if ( argc <= 1 )
   {
      SCIPinfoMessage(scip, nullptr, "%s\n", usage);
      return FALSE;
   }

   if ( *(argv[1]) == '-' )
   {
      SCIPinfoMessage(scip, nullptr, "Error: Filename should not start with '-'.\n\n");
      SCIPinfoMessage(scip, nullptr, "%s\n", usage);
      return FALSE;
   }

   *filename = argv[1];

   // check all remaining arguments
   for (i = 2; i < argc; ++i)
   {
      // check for settings file
      if ( strcmp(argv[i], "-s") == 0 )
      {
         ++i;  /*lint !e850*/
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "Error: No setting file name supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         *settingsname = argv[i];
      }
      // check for solution file
      else if ( strcmp(argv[i], "-l") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "Error: No solution file name supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         *solutionfile = argv[i];
      }
      // check for time limit
      else if ( strcmp(argv[i], "-t") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No time limit supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         timelimit = atof(argv[i]);
      }
      // check for memory limit
      else if ( strcmp(argv[i], "-m") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No memory limit supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         memlimit = atof(argv[i]);
      }
      // check for node limit
      else if ( strcmp(argv[i], "-n") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No node limit supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         nodelimit = atol(argv[i]);
      }
      // check for display frequency
      else if ( strcmp(argv[i], "-d") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No display frequency supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         dispfreq = atoi(argv[i]);
      }
      // check for display verbosity
      else if ( strcmp(argv[i], "-v") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No verbosity level supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         verblvl = atoi(argv[i]);
      }
      // check for cutoff value
      else if ( strcmp(argv[i], "-c") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No cutoff value supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         cutoffvalue = atof(argv[i]);
      }
      else
      {
	SCIPinfoMessage(scip, nullptr, "Unknown argument %s.\n\n", argv[i]);
	SCIPinfoMessage(scip, nullptr, "%s\n", usage);
        return FALSE;
      }
   }

   return TRUE;
}


//! main function for solving the stable set problem
int main(
   int                   argc,               //!< number of command line arguments
   const char**          argv                //!< array of command line arguments
   )
{
   // initialize SCIP
   SCIP* scip;
   SCIP_CALL_ERROR( SCIPcreate(&scip) );
   SCIP_CALL_ERROR( SCIPincludeDefaultPlugins(scip) );

   // output set covering banner
   SCIPinfoMessage(scip, nullptr, "SC: set covering branch-and-cut\n\n");
   const char* filename;
   const char* settingsname;
   const char* solutionfile;
   SCIP_Real timelimit;
   SCIP_Real memlimit;
   SCIP_Longint nodelimit;
   int dispfreq;
   int verblvl;
   SCIP_Real cutoffvalue;

   if ( ! readArguments(scip, argc, argv, &filename, &settingsname, &solutionfile, timelimit, memlimit, nodelimit, dispfreq, verblvl, cutoffvalue) )
   {
      SCIP_CALL_ERROR( SCIPfree(&scip) );
      BMScheckEmptyMemory();
      return -1;
   }

   if ( ! SCIPfileExists(filename) )
   {
      SCIPerrorMessage("Problem file <%s> does not exist.\n", filename);
      return SCIP_NOFILE;
   }

   // output SCIP banner
   SCIPprintVersion(scip, nullptr);
   SCIPinfoMessage(scip, nullptr, "\n");
   SCIPprintExternalCodes(scip, nullptr);
   SCIPinfoMessage(scip, nullptr, "\n");

   // initialize problem class
   std::string problem_name = get_problem_name(filename);

   SETCOVERING* setcovering = new SETCOVERING(problem_name.c_str(), filename, scip);
   SCIP_CALL_ERROR( SCIPcreateObjProb(scip, problem_name.c_str(), setcovering, TRUE) );
   SCIP_CALL_ERROR( setcovering->includePlugins() );

#ifdef DEBUGSOL
   SCIPenableDebugSol(scip);
#endif

   if ( argc == 3 )
   {
      setcovering->readSolution(argv[2]);
   }

   // check for parameters
   if ( settingsname != nullptr )
   {
      if ( SCIPfileExists(settingsname) )
      {
         SCIPinfoMessage(scip, nullptr, "Reading parameter file <%s> ...\n\n", settingsname);
         SCIP_CALL_ERROR( SCIPreadParams(scip, settingsname) );
      }
      else
      {
         SCIPerrorMessage("Parameter file <%s> not found.\n", settingsname);
         return SCIP_NOFILE;
      }
   }

   // perform graph presolving
   struct tms timer_beg;
   struct tms timer_end;
   (void) times(&timer_beg);
   SCIP_CALL_ERROR( setcovering->performSymmetryPresolving() );
   (void) times(&timer_end);
   double presolvingtime = (timer_end.tms_utime - timer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Presolving time: %.2f\n", presolvingtime);

   // solve the problem ...
   SCIP_CALL_ERROR( setcovering->initializeOptimization(nullptr) );

   // SCIP_CALL_ERROR( SCIPwriteOrigProblem(scip, "debug.cip", "cip", FALSE) );

   // set time, node, and memory limit
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/time", timelimit) );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Adjusting the solving time to: %.2f\n", presolvingtime);
   changeSolvingTime(scip, presolvingtime);

   if ( ! SCIPisInfinity(scip, memlimit) )
   {
      SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/memory", memlimit) );
   }

   if ( nodelimit < SCIP_LONGINT_MAX )
   {
      SCIP_CALL_ERROR( SCIPsetLongintParam(scip, "limits/nodes", nodelimit) );
   }

   if ( dispfreq >= 0 )
   {
      SCIP_CALL_ERROR( SCIPsetIntParam(scip, "display/freq", dispfreq) );
   }

   if ( verblvl >= 0 && verblvl <= 5 )
   {
      SCIP_CALL_ERROR( SCIPsetIntParam(scip, "display/verblevel", verblvl) );
   }

   // set cutoff and turn heuristics off
   if ( ! SCIPisInfinity(scip, -cutoffvalue) )
   {
      SCIPinfoMessage(scip, nullptr, "set limits objective %g\n", cutoffvalue);
      SCIP_CALL_ERROR( SCIPsetObjlimit(scip, cutoffvalue) );
      SCIPinfoMessage(scip, nullptr, "Turn all heuristics off!\n");
      SCIP_CALL_ERROR( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   }

   // solve problem
   SCIP_CALL_ERROR( SCIPsolve(scip) );
   SCIP_CALL_ERROR( SCIPprintBestSol(scip, NULL, FALSE) );

   SCIP_CALL_ERROR( setcovering->terminateOptimization() );

   /*
   SCIPprintMemoryDiagnostic(scip);
   */

   SCIP_CALL_ERROR( SCIPfree(&scip) );
}
