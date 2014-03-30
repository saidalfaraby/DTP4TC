# $Id: Makefile,v 2.2 2003/09/03 16:50:37 jhoey Exp $
#
# pair: 
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

PACKAGE	= terminal
PSRC	= terminal.cc 
PHDR    = terminal.hh 
POBJ	= $(PSRC:.cc=.o)

LIBRARY	= ../lib$(PACKAGE).a

# files for the test program
EXEC	= 
EXECSRC	= $(EXEC).cc
EXECOBJ	= $(EXECSRC:.cc=.o)

#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all: $(POBJ) 

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
.SUFFIXES: .o .cc

# rules to build .o's from .c's
%.o: %.cc $(PHDR) $(DEPENDHDRS)
	$(CXX) -I../../../include/ -DUNIX -c -o $@ $<

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(LIBDIR)/$(LIBRARY)
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))

# include the general rules files to customize this makefile
include ../../Makefile.rules
include ../../Makefile.rules.$(BUILD)



