# $Id: Makefile,v 2.1 2004/05/25 17:11:42 jhoey Exp $
#
#    Cudd - DD package
#---------------------------

# Heavily modified by Don Murray, Sept 2000

# set the debug level desired
#DDDEBUG	= -DDD_DEBUG -DDD_CACHE_PROFILE -DDD_VERBOSE -DDD_UNIQUE_PROFILE
DDDEBUG	=

CFLAGS	+=  -O4	$(DDDEBUG)

LIBS	= -lcudd -lmtr -lst -lutil -lterminal

#MNEMLIB =
#LDFLAGS	=

# files for the package
PACKAGE	= cudd
PSRC	= cuddAPI.c cuddAddAbs.c cuddAddApply.c cuddAddFind.c cuddAddIte.c \
	  cuddAddInv.c cuddAddNeg.c cuddAddWalsh.c cuddAndAbs.c \
	  cuddAnneal.c cuddApa.c cuddApprox.c cuddBddAbs.c cuddBddCorr.c \
	  cuddBddIte.c cuddBridge.c cuddCache.c cuddCheck.c cuddClip.c \
	  cuddCof.c cuddCompose.c cuddDecomp.c cuddEssent.c \
	  cuddExact.c cuddExport.c cuddGenCof.c cuddGenetic.c \
	  cuddGroup.c cuddHarwell.c cuddInit.c cuddInteract.c \
	  cuddLCache.c cuddLevelQ.c \
	  cuddLinear.c cuddLiteral.c cuddMatMult.c cuddPriority.c \
	  cuddRead.c cuddRef.c cuddReorder.c cuddSat.c cuddSign.c \
	  cuddSolve.c cuddSplit.c cuddSubsetHB.c cuddSubsetSP.c cuddSymmetry.c \
	  cuddTable.c cuddUtil.c cuddWindow.c cuddZddCount.c cuddZddFuncs.c \
	  cuddZddGroup.c cuddZddIsop.c cuddZddLin.c cuddZddMisc.c \
	  cuddZddPort.c cuddZddReord.c cuddZddSetop.c cuddZddSymm.c \
	  cuddZddUtil.c
PHDR    = cudd.h cuddInt.h
POBJ	= $(PSRC:.c=.o)

LIBRARY	= lib$(PACKAGE).a

# files for the test program
EXEC	= test$(PACKAGE)
EXECSRC	= $(EXEC).c
EXECOBJ	= $(EXECSRC:.c=.o)

#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all: $(LIBRARY) 

install: installinclude installlib installbin

# copy include files to the include directory
installinclude: 
	$(INSTALL) $(INSTALLFLAGS) $(PHDR) $(INCLUDEDIR)

# build library and copy to lib directory
installlib: $(LIBRARY)
	$(INSTALL) $(INSTALLFLAGS) $(LIBRARY) $(LIBDIR)

# build executables and copy to the bin directory
installbin: $(EXEC)
	$(INSTALL) $(INSTALLFLAGS) $(EXEC) $(BINDIR)

$(LIBRARY): $(POBJ)
	ar rv $@ $^
	$(RANLIB) $@

# file suffixes to pay attention to
.SUFFIXES: .o .c

# rules to build .o's from .c's
# note: use CXX since now includes .hh files
%.o: %.c $(PHDR)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXEC): $(EXECSRC) $(EXECOBJ) $(PHDR) $(MNEMLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EXECOBJ) $(LIBS) $(MNEMLIB) -lm

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(BINDIR)/$(EXEC) $(LIBDIR)/$(LIBRARY)
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))

# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)

