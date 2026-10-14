# setcovering
Solving set covering problems using symmetry presolving

This project contains code to solve set covering instances, based on
SCIP (https://github.com/scipopt/scip). To compile, define the
enivornment variable SCIP_PATH to the SCIP directory (SCIP needs to be
compiled using make, not cmake). The the code can be compiled in the
same way as SCIP (see the installation instructions there).

The code is then placed in the "bin/" directory. It is run by passing
an LP- oder MPS-file that contains a set covering instances. Settings
can be passed with "-s <setting file>". For example,

bin/setCovering <file> -s settings/presolve.set

This enables symmetry presolving as explained in https://optimization-online.org/?p=32590.
