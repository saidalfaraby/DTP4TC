
# $Id: Makefile,v 1.12 2005/02/17 01:34:58 jhoey Exp $
#
# Makefile for the CUDD distribution kit 
#---------------------------------------------------------------------------


#
# Heavily modified from Fabio's original distribution by 
# Don Murray, Sept, 2000
#

# packages to compile

PKGS	 	= MDP \
		MVPSpudd
#------------------------------------------------------------------------

all: installinclude installlib installbin

installinclude installlib installbin clean distclean depend:
	for pkg in $(PKGS) ; do \
		( echo $$pkg:$@ ; cd $$pkg ; $(MAKE) $(MFLAGS) $@ ) ;\
	done


