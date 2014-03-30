
# $Id: Makefile,v 2.0 2003/02/07 00:11:54 staubin Exp $
#
# Makefile for the CUDD distribution kit 
#---------------------------------------------------------------------------


#
# Heavily modified from Fabio's original distribution by 
# Don Murray, Sept, 2000
#

# packages to compile
#PKGS	 	= cudd	
PKGS	 	= cudd	\
                  dddmp-2.0 \
		  mtr	\
		  st	\
		  util	\
		  terminal \
		  nanotrav \
		  epd 

#------------------------------------------------------------------------

all: installinclude installlib installbin

installinclude installlib installbin clean distclean depend:
	for pkg in $(PKGS) ; do \
		( echo $$pkg:$@ ; cd $$pkg ; $(MAKE) $(MFLAGS) $@ ) ;\
	done


