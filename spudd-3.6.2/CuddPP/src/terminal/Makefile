
# $Id: Makefile,v 2.1 2004/12/14 13:51:06 jhoey Exp $
#
# Makefile for the CUDD distribution kit 
#---------------------------------------------------------------------------


#
# Heavily modified from Fabio's original distribution by 
# Don Murray, Sept, 2000
#

PACKAGE	= terminal
LIBRARY	= lib$(PACKAGE).a

# set default build to linux or solaris
#BUILD	= linux
#BUILD	= solaris

# define where to find the include and library files
TOP		= ../..
INCLUDEDIR	= $(TOP)/include
LIBDIR		= $(TOP)/lib/$(BUILD)
BINDIR		= $(TOP)/bin/$(BUILD)

RANLIB	= ranlib

INSTALLFLAGS 	+= -p --mode=444
INSTALL		= install


# packages to compile
PKGS	 	= terminal \
		    pair    \
	             gaussian
#------------------------------------------------------------------------

all: installinclude installbin installlib


installinclude installbin clean distclean depend:
	for pkg in $(PKGS) ; do \
		( echo $$pkg:$@ ; cd $$pkg ; $(MAKE) $(MFLAGS) $@ ) ;\
	done

installlib:
	for pkg in $(PKGS) ; do \
		( echo $$pkg:$@ ; cd $$pkg ; $(MAKE) $(MFLAGS) $@ ) ;\
	done
	for pkg in $(PKGS) ; do \
		( ar rv $(LIBRARY) $$pkg/*.o ) ;\
	done
	$(INSTALL) $(INSTALLFLAGS) $(LIBRARY) $(LIBDIR)

# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)
