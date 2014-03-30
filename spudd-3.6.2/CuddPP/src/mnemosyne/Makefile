# $Id: Makefile,v 2.0 2003/02/07 00:13:29 staubin Exp $
#
#	Makefile for the Mnemosyne memory allocation tracker.
#
#	Marcus J. Ranum, 1990
#

# Heavily modified by Don Murray, Sept 2000

#Options:
# define MALLOC_IS_VOIDSTAR if your system's malloc is declared as a (void *)
# otherwise, it is assumed to be a (char *). a "mall_t" is typedeffed in
# mnemconf.h and mnemosyne.h to implement this.
#OPTNS	=	-DMALLOC_IS_VOIDSTAR

CPPFLAGS	+= -DMALLOC_IS_VOIDSTAR

LIBS	= -lmnem 

# files for the package
PACKAGE	= mnem
PSRC	= mnemosyne.c
PHDR	= mnemosyne.h mnemconf.h
POBJ	= $(PSRC:.c=.o)

LIBRARY	= lib$(PACKAGE).a

# files for the test program
EXEC1		= mnemalyse
EXEC1SRC	= $(EXEC1).c
EXEC1OBJ	= $(EXEC1SRC:.c=.o)
EXEC2		= mtest
EXEC2SRC	= $(EXEC2).c
EXEC2OBJ	= $(EXEC2SRC:.c=.o)

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
installbin: $(EXEC1) $(EXEC2)
	$(INSTALL) $(INSTALLFLAGS) $(EXEC1) $(BINDIR)
	$(INSTALL) $(INSTALLFLAGS) $(EXEC2) $(BINDIR)

$(LIBRARY): $(POBJ)
	ar rv $@ $^
	$(RANLIB) $@

# file suffixes to pay attention to
.SUFFIXES: .o .c

# rules to build .o's from .c's
%.o: %.c $(PHDR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXEC1): $(EXEC1SRC) $(EXEC1OBJ) $(PHDR) $(MNEMLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EXEC1OBJ) $(LIBS) $(MNEMLIB) -lm

$(EXEC2): $(EXEC2SRC) $(EXEC2OBJ) $(PHDR) $(MNEMLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EXEC2OBJ) $(LIBS) $(MNEMLIB) -lm

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(BINDIR)/$(EXEC) $(LIBDIR)/$(LIBRARY)
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))


# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)

