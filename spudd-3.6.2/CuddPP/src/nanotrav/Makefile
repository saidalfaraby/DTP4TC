# $Id: Makefile,v 2.0 2003/02/07 00:13:43 staubin Exp $
#
# nanotrav: simple FSM traversal program
#---------------------------------------------------------------------------

# Heavily modified by Don Murray, Sept 2000

# set the debug level desired
#DDDEBUG	= -DDD_DEBUG -DDD_CACHE_PROFILE -DDD_VERBOSE -DDD_UNIQUE_PROFILE
DDDEBUG	=

CFLAGS	+=	$(DDDEBUG)

LIBS	= -ldddmp -lcudd -lmtr -lst -lutil -lterminal

MNEMLIB =
#MNEMLIB	= -lmnem

# files for the package
PACKAGE	= nanotrav
PSRC	= main.c bnet.c ntr.c ntrHeap.c ntrBddTest.c ntrMflow.c ntrZddTest.c \
	  chkMterm.c ucbqsort.c 
# no install includes for this package
PHDR    = 
POBJ	= 
DEPENDHDR	= bnet.h		\
		  ntr.h			\
		  $(INCDIR)/cuddInt.h 	\
		  $(INCDIR)/cudd.h	\
	 	  $(INCDIR)/dddmp.h

# no library for this package
LIBRARY	= 

# files for executable
EXEC	= $(PACKAGE)
EXECSRC	= $(PSRC)
EXECOBJ	= $(EXECSRC:.c=.o)

#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all: $(EXEC) 

install: installinclude installbin
# installlib - don't include as this package doesn't have a library

# copy include files to the include directory
installinclude: 
#	$(INSTALL) $(INSTALLFLAGS) $(PHDR) $(INCLUDEDIR)

# build library and copy to lib directory
installlib: $(LIBRARY)
# no library for this one
#	$(INSTALL) $(INSTALLFLAGS) $(LIBRARY) $(LIBDIR)

# build executables and copy to the bin directory
installbin: $(EXEC)
	$(INSTALL) $(INSTALLFLAGS) $(EXEC) $(BINDIR)

$(LIBRARY): $(POBJ)
	ar rv $@ $^
	$(RANLIB) $@

# file suffixes to pay attention to
.SUFFIXES: .o .c

# rules to build .o's from .c's
%.o: %.c $(PHDR) $(DEPENDHDRS)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXEC): $(EXECSRC) $(EXECOBJ) $(PHDR) $(MNEMLIB)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $(EXECOBJ) $(LIBS) $(MNEMLIB) -lm

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(BINDIR)/$(EXEC) 
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))


# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)





