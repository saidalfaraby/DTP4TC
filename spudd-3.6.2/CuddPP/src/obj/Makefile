# $Id: Makefile,v 2.0 2003/02/07 00:13:51 staubin Exp $
#
# obj: CUDD's C++ interface
#---------------------------------------------------------------------------

# Heavily modified by Don Murray, Sept 2000

# set the debug level desired
#DDDEBUG	= -DDD_DEBUG -DDD_CACHE_PROFILE -DDD_VERBOSE -DDD_UNIQUE_PROFILE
DDDEBUG	=

CFLAGS	+=	$(DDDEBUG)

LIBS	= -lobj -lcudd -lmtr -lst -lutil -lterminal

#MNEMLIB =
#LDFLAGS	=

# files for the package
PACKAGE	= obj
PSRC	= cuddObj.cc
PHDR	= cuddObj.hh 
POBJ	= $(PSRC:.cc=.o)
DEPENDHDRS	= $(INCDIR)/cudd.h

LIBRARY	= lib$(PACKAGE).a

# files for the test program
EXEC	= test$(PACKAGE)
# exec file seems to be broken
EXECSRC	= $(EXEC).cc
EXECOBJ	= $(EXECSRC:.cc=.o)

#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all: $(LIBRARY) $(EXEC)

install: installinclude installlib installbin

# copy include files to the include directory
installinclude: 
	$(INSTALL) $(INSTALLFLAGS) $(PHDR) $(INCLUDEDIR)

# build library and copy to lib directory
installlib: $(LIBRARY)
	$(INSTALL) $(INSTALLFLAGS) $(LIBRARY) $(LIBDIR)

# build executables and copy to the bin directory
installbin: $(EXEC)
#	$(INSTALL) $(INSTALLFLAGS) $(EXEC) $(BINDIR)
# exec file seems to be broken

$(LIBRARY): $(POBJ)
	ar rv $@ $^
	$(RANLIB) $@

# file suffixes to pay attention to
.SUFFIXES: .o .cc

# rules to build .o's from .c's
%.o: %.cc $(PHDR) $(DEPENDHDRS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXEC): $(EXECSRC) $(EXECOBJ) $(PHDR) $(MNEMLIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(EXECOBJ) $(LIBS) $(MNEMLIB) -lm

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
#	@\rm -f $(BINDIR)/$(EXEC) 
	@\rm -f $(LIBDIR)/$(LIBRARY)
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))


# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)


