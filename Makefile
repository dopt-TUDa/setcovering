# -------------------------------------------------------------------------------
# Makefile for setCovering
# -------------------------------------------------------------------------------

#----------------------------------------------------------------------------
# paths
#----------------------------------------------------------------------------

# SCIP_PATH is an enivironment variable that should be set outside
# this makefile and point to the main directory of SCIP.

SCIPDIR         =       $(SCIP_PATH)
SCIPREALPATH	=	$(realpath $(SCIPDIR))

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project
-include $(SCIPDIR)/make/local/make.$(HOSTNAME)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP).$(OPT)

#-----------------------------------------------------------------------------
# basic definitions
#-----------------------------------------------------------------------------

VERSION		=	1.0

# overwrite flags for dependencies
DFLAGS          =       -MMD

# detect with symmetry code has been used by SCIP
SCIPVERSION	:=	$(shell $(SCIPDIR)/bin/scip.$(BASE).$(LPS).$(TPI)$(EXEEXTENSION) -v | sed -e 's/$$/@/')
override SYM:=	$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')

# detect symmetry handling code
ifeq ($(SYM),nauty)
USRFLAGS        +=      -DSC_WITH_NAUTY
USRCXXFLAGS     +=      $(CXX14FLAG)
ifeq ($(NAUTYEXTERNAL),true)
USRLDFLAGS      +=      $(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)nauty.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif
endif

# use pure nauty here
ifeq ($(SYM),snauty)
USRFLAGS        +=      -DSC_WITH_NAUTY
USRCXXFLAGS     +=      $(CXX14FLAG)
ifeq ($(NAUTYEXTERNAL),true)
USRLDFLAGS      +=      $(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)nauty.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif
endif

ifeq ($(SYM),bliss)
USRFLAGS        +=      -DSC_WITH_BLISS -I$(SCIPDIR)/$(LIBDIR)/include/
USRCXXFLAGS     +=      $(CXX14FLAG)
USRLDFLAGS      +=      $(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)bliss.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif

ifeq ($(SYM),sbliss)
USRFLAGS        +=      -DSC_WITH_SBLISS -I$(SCIPDIR)/$(LIBDIR)/include/
USRCXXFLAGS     +=      $(CXX14FLAG)
USRLDFLAGS      +=      $(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)bliss.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif


#-----------------------------------------------------------------------------
# main program
#-----------------------------------------------------------------------------

MAINNAME	=	setCovering

MAINCXXOBJ 	= 	setCovering.o get_problem_name.o branch_most_covered.o branch_orbital.o \
			branch_intersection.o branch_half_constraint.o Bitset.o event_sol.o symmetry.o main.o
MAINCOBJ 	=	sepa_bang.o changesoltime.o

MAINSRC 	=       $(addprefix $(SRCDIR)/,$(MAINCXXOBJ:.o=.C))
MAINSRC 	+=	$(addprefix $(SRCDIR)/,$(MAINCOBJ:.o=.c))

MAIN            =       $(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE        =       $(BINDIR)/$(MAIN)
MAINSHORTLINK   =       $(BINDIR)/$(MAINNAME)
MAINOBJFILES    =       $(addprefix $(OBJDIR)/,$(MAINCXXOBJ)) $(addprefix $(OBJDIR)/,$(MAINCOBJ))

#-----------------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------------

# running the makefile without arguments will build everything
.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

# create directory for objective files
$(OBJDIR):
		@-mkdir -p $(OBJDIR)

# create directory for binaries
$(BINDIR):
		@-mkdir -p $(BINDIR)

# possibly compile SCIP
.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

# call pclint code checker on source files
.PHONY: pclint
pclint:		$(MAINSRC)
		-rm -f pclint.out
ifeq ($(FILES),)
		@$(SHELL) -ec 'echo "-> running pclint ..."; \
			for i in $^; \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
else
		@$(SHELL) -ec  'echo "-> running pclint on specified files ..."; \
			for i in $(FILES); \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
endif

# clean all objective files
.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rm -f $(OBJDIR)/*.d && rmdir $(OBJDIR));
		@echo "-> remove objective files"
endif

# create tags
.PHONY: tags
tags:
		rm -f TAGS; etags src/*.C src/*.h src/*.c $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h --output=TAGS; sed 's!\#undef .*!!g' TAGS > tags; mv tags TAGS

# run tests
.PHONY: test
test:		$(MAINFILE)
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) \
		$(MEM) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(USEBESTSOL) $(SETCUTOFF);

# run tests on cluster
.PHONY: testcluster
testcluster:	$(MAINFILE)
		cd check; \
		$(SHELL) ./check_cluster.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) \
		$(MEM) $(DISPFREQ) $(VERSION) $(LPS) $(SETCUTOFF) $(QUEUE) $(PPN) $(CLIENTTMPDIR) $(EXCLUSIVE);

#-----------------------------------------------------------------------------
# project rules
#-----------------------------------------------------------------------------

# do not attempt to include .d files if there will definitely be any (empty DFLAGS), because it slows down the build on Windows considerably
ifneq ($(DFLAGS),)
-include $(MAINOBJFILES:.o=.d)
else
ifeq ($(VERBOSE),true)
$(info No compilation dependencies. If changing header files, do a make clean before building.)
endif
endif

# turn off output unless VERBOSE is true
ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

# main targets:
$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBBASEFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> using SYM=$(SYM)"
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) $(LINKCXXSCIPALL) $(LINKCXX_o)$@

# create short link
$(MAINSHORTLINK): $(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

generateRandomSC: $(BINDIR)
		$(CXX) $(CCFLAGS) -o $(BINDIR)/generateRandomSC $(SRCDIR)/generateRandomSC.C


# rules to build objective files from source code
$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(DFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.C
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) $(DFLAGS) -c $< $(CXX_o)$@
