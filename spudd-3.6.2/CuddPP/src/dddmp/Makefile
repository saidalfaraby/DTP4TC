# $Id: Makefile,v 2.0 2003/02/07 00:13:13 staubin Exp $
#
# dddmp: Decision Diagram storage and retrieval
#---------------------------------------------------------------------------

# Heavily modified by Don Murray, Sept 2000

# set the debug level desired
#DDDEBUG	= -DDD_DEBUG -DDD_CACHE_PROFILE -DDD_VERBOSE -DDD_UNIQUE_PROFILE
DDDEBUG	=

CFLAGS	+=	$(DDDEBUG)

LIBS	= -ldddmp -lcudd -lmtr -lst -lutil -lterminal

#MNEMLIB =
#LDFLAGS =

# files for the package
PACKAGE	= dddmp
PSRC    = dddmpDump.c dddmpBinary.c dddmpConvert.c dddmpDbg.c \
          dddmpDdNode.c
PHDR    = dddmp.h dddmpInt.h 
# note: depend headers not checked now
DEPENDHDR= $(INCLUDEDIR)/cudd.h $(INCLUDEDIR)/cuddInit.h
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
# use CXX since this package is including .hh files
%.o: %.c $(PHDR) 
	@echo "building $@"
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

