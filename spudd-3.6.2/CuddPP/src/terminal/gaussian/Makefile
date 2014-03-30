# $Id: Makefile,v 1.2 2005/03/07 22:24:36 jhoey Exp $
#
# gaussian:
#---------------------------------------------------------------------------

# Heavily modified by Don Murray, Sept 2000

# set the debug level desired
#DDDEBUG	= -DDD_DEBUG -DDD_CACHE_PROFILE -DDD_VERBOSE -DDD_UNIQUE_PROFILE
DDDEBUG	=

CXXFLAGS	+=	$(DDDEBUG)
CXXFLAGS	+=	-DDD_STATS

LIBS	= -lcudd -lmtr -lst -lutil

#MNEMLIB =
#LDFLAGS	=

# files for the package

PACKAGE	= mixgauss
PSRC	= gaussian.cpp mixgauss.cpp nrutil.cpp
PHDR    = gaussian.h mixgauss.h nrutil.h
POBJ	= $(PSRC:.cpp=.o)

LIBRARY	= ../lib$(PACKAGE).a

# files for the test program
EXEC	= testmixgauss
EXECSRC	= $(EXEC).cpp  $(PSRC)
EXECOBJ	= $(EXECSRC:.cpp=.o) 

# mathematica stuff (if it exists)
MMA_VERSION=5.1
MMA_MLINKDIR = /h/23/jhoey/Mathematica/${MMA_VERSION}/AddOns/MathLink/DeveloperKit
MMA_SYS = Linux
MMA_CADDSDIR = ${MMA_MLINKDIR}/${MMA_SYS}/CompilerAdditions
MMA_INCDIR = ${MMA_CADDSDIR}
MMA_LIBDIR = ${MMA_CADDSDIR}

MMA_IFLAGS =  -I$(MMA_INCDIR)
MMA_LFLAGS = -L$(MMA_LIBDIR)
#MMA_LIBS = -lML -lrt

#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all: $(POBJ) $(EXEC)

# not including installbin since this package has not executable
install: installinclude installlib 

# copy include files to the include directory
installinclude: 
	$(INSTALL) $(INSTALLFLAGS) $(PHDR) ../$(INCLUDEDIR)

# build library and copy to lib directory
installlib: $(POBJ)
#	$(INSTALL) $(INSTALLFLAGS) $(LIBRARY) $(LIBDIR)
# no action since its done in terminal Makefile

# build executables and copy to the bin directory
installbin: $(EXEC)
#	$(INSTALL) $(INSTALLFLAGS) $(EXEC) $(BINDIR)
# no action for installbin since this package has not executable

$(LIBRARY): $(POBJ)
#	ar rv $@ $^
#	$(RANLIB) $@

# file suffixes to pay attention to
.SUFFIXES: .o .cpp

# rules to build .o's from .c's
%.o: %.cpp $(PHDR) $(DEPENDHDRS)
	$(CXX) -I../../../include/ $(MMA_IFLAGS) -DUNIX -c -o $@ $<

$(EXEC): $(EXECOBJ) $(PHDR) 
	$(CC) $(CFLAGS) -L ./ $(MMA_LFLAGS) -o $@ $(EXECOBJ) $(MMA_LIBS)

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(LIBDIR)/$(LIBRARY)
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))

# include the general rules files to customize this makefile
include ../../Makefile.rules
include ../../Makefile.rules.$(BUILD)



