# $Id: Makefile,v 2.0 2003/02/07 00:14:10 staubin Exp $
#
# util -- miscellaneous utility routines
#---------------------------------------------------------------------------

# Heavily modified by Don Murray, Sept 2000

# set the debug level desired
#DDDEBUG	= -DDD_DEBUG -DDD_CACHE_PROFILE -DDD_VERBOSE -DDD_UNIQUE_PROFILE
DDDEBUG	=

CFLAGS	+=	$(DDDEBUG)

LIBS	= -lcudd -lmtr -lst -lutil

#MNEMLIB =
#LDFLAGS	=

# files for the package
PACKAGE	= util
PSRC	= cpu_time.c cpu_stats.c getopt.c safe_mem.c strsav.c texpand.c \
	  ptime.c prtime.c pipefork.c pathsearch.c stub.c \
	  restart.c tmpfile.c datalimit.c
PHDR	= util.h
POBJ	= $(PSRC:.c=.o)


LIBRARY	= lib$(PACKAGE).a

# files for the test program
EXEC	= 
EXECSRC	= 
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
#	$(INSTALL) $(INSTALLFLAGS) $(EXEC) $(BINDIR)
# no executable for this package

$(LIBRARY): $(POBJ)
	ar rv $@ $^
	$(RANLIB) $@

# file suffixes to pay attention to
.SUFFIXES: .o .c

# rules to build .o's from .c's
# use CXX since there is a problem with declaring macro
# strchr() on line 39
%.o: %.c $(PHDR)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

# if the header files change, recompile
# I don't think this is necessary
#$(POBJ): $(PHDR)
#$(OBJ): $(PHDR)

$(EXEC): $(EXECSRC) $(EXECOBJ) $(PHDR) $(MNEMLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EXECOBJ) $(LIBS) $(MNEMLIB) -lm

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(LIBDIR)/$(LIBRARY)
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))

# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)





