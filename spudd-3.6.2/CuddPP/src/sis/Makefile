


NOT FINISHED, IS THIS PACKAGE NECESSARY??

# $Id: Makefile,v 2.0 2003/02/07 00:13:55 staubin Exp $
#
#    Cudd - sis package
#---------------------------

# Heavily modified by Don Murray, Sept 2000

# set the debug level desired
#DDDEBUG	= -DDD_DEBUG -DDD_CACHE_PROFILE -DDD_VERBOSE -DDD_UNIQUE_PROFILE
DDDEBUG	=

CFLAGS	+=	$(DDDEBUG)

LIBS	= -lcudd -lmtr -lst -lutil

#MNEMLIB =
#LDFLAGS	=

# files for the package
PACKAGE	= bdd
PSRC	= cuddBdd.h  cuddBddPort.c  cuddPwPt.c 
PHDR    = cuddBdd.h 
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
%.o: %.c $(PHDR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

# if the header files change, recompile
# I don't think this is necessary
#$(POBJ): $(PHDR)
#$(OBJ): $(PHDR)

$(EXEC): $(EXECSRC) $(EXECOBJ) $(PHDR) $(MNEMLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EXECOBJ) $(LIBS) $(MNEMLIB) -lm

# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)



CAD	= /projects/octtools/octtools/$(MACHINE)
SIS	= ..
LINTCREATEFLAG = -C

# files for the package
P	= bdd
PSRC	= cuddAPI.c cuddAddAbs.c cuddAddApply.c cuddAddFind.c cuddAddIte.c \
	  cuddAddInv.c cuddAddNeg.c cuddAddWalsh.c cuddAndAbs.c \
	  cuddAnneal.c cuddApa.c cuddApprox.c cuddBddAbs.c cuddBddCorr.c \
	  cuddBddIte.c cuddBddPort.c cuddBridge.c cuddCache.c cuddCheck.c \
	  cuddClip.c cuddCof.c cuddCompose.c cuddDecomp.c cuddEssent.c \
	  cuddExact.c cuddExport.c cuddGenCof.c \
	  cuddGenetic.c cuddGroup.c cuddHarwell.c cuddInit.c cuddInteract.c \
	  cuddLCache.c cuddLevelQ.c cuddLinear.c cuddLiteral.c \
	  cuddMatMult.c cuddPriority.c cuddPwPt.c \
	  cuddRead.c cuddRef.c cuddReorder.c cuddSat.c cuddSign.c \
	  cuddSolve.c cuddSplit.c cuddSubsetHB.c cuddSubsetSP.c cuddSymmetry.c \
	  cuddTable.c cuddUtil.c cuddWindow.c cuddZddCount.c cuddZddFuncs.c \
	  cuddZddGroup.c cuddZddIsop.c cuddZddLin.c cuddZddMisc.c \
	  cuddZddPort.c cuddZddReord.c cuddZddSetop.c cuddZddSymm.c \
	  cuddZddUtil.c
POBJ	= $(PSRC:.c=.o)
PHDR    = cudd.h cuddInt.h cuddBdd.h

# files for the test program
TARGET	= testcudd
SRC	= testcudd.c
OBJ	= $(SRC:.c=.o)
HDR	=

LIBS	= ../util/libutil.a ../st/libst.a
LINTLIBS= ../util/llib-lutil.ln ../st/llib-lst.ln
INCLUDE = -I$(CAD)/include -I$(SIS)/include

CFLAGS	= $(INCLUDE)
LDFLAGS = 
LINTFLAGS = $(INCLUDE) ${LINTEXTRAS}

#------------------------------------------------------

$(TARGET): $(PHDR) $(OBJ) $(POBJ) $(LIBS)
	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJ) $(POBJ) $(LIBS)

lint: $(PSRC) $(PHDR) $(SRC) $(HDR)
	lint $(LINTFLAGS) $(SRC) $(PSRC) $(LINTLIBS)

install: lib$(P).a llib-l$(P).ln

lib$(P).a: $(POBJ)
	ar cr $@ $?
	$(RANLIB) $@

unpack: lib$(P).a
	@for i in $(POBJ); do				\
	    ln -s $(SIS)/$(P)/$$i $(SIS)/unpack;	\
	done

llib-l$(P).ln: $(PSRC) $(PHDR)
	lint $(LINTFLAGS) ${LINTCREATEFLAG}$(P) -n $(PSRC)

clean:
	rm -f $(TARGET) *.a *.ln *.o \
	    [Tt]ags [Mm]ake.out lint malloc.out gmon.out __.SYMDEF

tags: _force
	@for i in $(PSRC) $(PHDR); do			\
	    cwd=`pwd`; ctags -a $$cwd/$$i;		
	done;

strip_depend:
	sed '/^#--DO NOT CHANGE ANYTHING AFTER THIS LINE/,$$d' Makefile >mktemp
	mv mktemp Makefile

depend:
	sed '/^#--DO NOT CHANGE ANYTHING AFTER THIS LINE/,$$d' Makefile >mktemp
	echo '#--DO NOT CHANGE ANYTHING AFTER THIS LINE' >>mktemp
	$(CAD)/bin/cc-M $(CFLAGS) $(PSRC) | \
	    sed 's|$(CAD)|$$(CAD)|g' | \
	    grep -v '/usr/include' >>mktemp
	mv mktemp Makefile

#-------------------------- IBM 3090 support -----------------
IBMHOST = opua
IBMDIST = /users2/sis
ibmdist: $(PSRC) $(PHDR)
	rdist -Richw $(PSRC) $(PHDR) $(IBMHOST):$(IBMDIST)
#-------------------------------------------------------------
_force:

