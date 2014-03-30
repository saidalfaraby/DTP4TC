# $Id: Makefile,v 2.9 2004/04/15 13:16:34 jhoey Exp $
#
#    Spudd - Spudd package
#---------------------------

# Heavily modified by Don Murray, Sept 2000

LIBS	= -lterminal -lcudd -lmtr -lst -lutil -ldddmp -lm

# files for the package
PACKAGE	= MVPSpudd
PCPPSRC	= JAPRICODD.cc pspuddproc.cc PSpudd.cc
PCSRC	= dumpdot.c mvdumpdot.c parser.tab.c parser.c 
PHDR    = pspudd.hh
POBJ	= $(PCSRC:.c=.o) $(PCPPSRC:.cc=.o) PSpuddMain.o 
P1OBJ   = $(PCSRC:.c=.o) $(PCPPSRC:.cc=.o) policyQuery.o policyQueryMain.o
P2OBJ   = $(PCSRC:.c=.o) $(PCPPSRC:.cc=.o) policyServerMain.o policyServer.o

LIBRARY	= 

# files for the program
EXEC	= $(PACKAGE)
EXECSRC	= 
EXECOBJ	= $(POBJ)

EXEC1	= pquery
EXEC1SRC = 
EXEC1OBJ = $(P1OBJ)

EXEC2	= pserve
EXEC2SRC = 
EXEC2OBJ = $(P2OBJ)


FLTKDIR = /usr/local/src/fltk-1.1.4
#FLTKDIR = /lci/project/little/jhoey/software/fltk-1.1.2

IFLAGS1 =  -I$(FLTKDIR)
IFLAGS2 =  -I$(FLTKDIR) -I../../../../homer/include/ 
LFLAGS1 = -L/usr/X11R6/lib  -L$(FLTKDIR)/lib/ 
LFLAGS2 = -L/usr/X11R6/lib  -L$(FLTKDIR)/lib/ -L../../../../homer/lib/ 
LIBS1 = -lfltk -lXext -lX11 -lXmu  -lXt -lXi 
LIBS2 = -lfltk -lXext -lX11 -lXmu  -lXt -lXi 

#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all:  $(EXEC) $(EXEC1)
install: installinclude installlib installbin

# copy include files to the include directory
installinclude: 
#	install -p --mode=444 $(PHDR) $(INCLUDEDIR)
# no includes

# build library and copy to lib directory
installlib: $(LIBRARY)
#	install -p --mode=444 $(LIBRARY) $(LIBDIR)
# no library for this package

# build executables and copy to the bin directory
installbin: $(EXEC) $(EXEC1)
	install -p --mode=555 $(EXEC) $(EXEC1) $(BINDIR)

$(LIBRARY): $(POBJ)
	ar rv $@ $^
	$(RANLIB) $@

# file suffixes to pay attention to
#.SUFFIXES: .o .c

# rules to build .o's from .c's
# note: use CXX since now includes .hh files
%.o: %.c $(PHDR)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXEC): $(EXECOBJ) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS)  -o $@ $(EXECOBJ) $(LIBS) $(MNEMLIB) -lm

$(EXEC1): $(EXEC1OBJ) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS1) -o $@ $(EXEC1OBJ) $(LIBS) $(LIBS1) $(MNEMLIB) -lm

$(EXEC2): $(EXEC2OBJ) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS2) -o $@ $(EXEC2OBJ) $(LIBS) $(LIBS2) $(MNEMLIB) -lm



distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(BINDIR)/$(EXEC) 
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))

# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)

# specify indivisual targets
parser.tab.c:parser.y pspudd.hh
	bison -d -v parser.y
parser.c:
	echo
parser.o:parser.c 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c parser.c
lex.yy.c:scanner.h scanner.l parser.tab.h
	flex scanner.l
parser.tab.o:parser.tab.c lex.yy.c pspudd.hh
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c parser.tab.c
JAPRICODD.o:JAPRICODD.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c JAPRICODD.cc
pspuddproc.o:pspuddproc.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c pspuddproc.cc
aspuddproc.o:aspuddproc.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c aspuddproc.cc
PSpudd.o:PSpudd.cc pspudd.hh
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c PSpudd.cc
PSpuddMain.o:PSpuddMain.cc pspudd.hh
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c PSpuddMain.cc
lex.yy.o:lex.yy.c scanner.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c lex.yy.c
mvdumpdot.o: mvdumpdot.c
	 $(CXX) $(CFLAGS) $(CPPFLAGS) -c mvdumpdot.c
dumpdot.o:dumpdot.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c dumpdot.c
policyQuery.o: policyQuery.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c policyQuery.cpp
policyQueryMain.o: policyQueryMain.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c policyQueryMain.cpp
policyServer.o: policyServer.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS2) -c policyServer.cpp
policyServerMain.o: policyServerMain.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS2) -c policyServerMain.cpp















