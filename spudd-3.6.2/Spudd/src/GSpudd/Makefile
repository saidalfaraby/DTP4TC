# $Id: Makefile,v 1.2 2005/03/07 22:24:20 jhoey Exp $
#
#    Spudd - Spudd package
#---------------------------

# Heavily modified by Don Murray, Sept 2000

LIBS	= -lterminal -lobj -lcudd -lmtr -lst -lutil -ldddmp -lm

# files for the package
PACKAGE	= gspudd
PCPPSRC	= gspudd.cpp MDP.cpp POMDP.cpp
PCSRC	= mdpparser.tab.c mvdumpdot.c
PHDR    = 
POBJ	= $(PCPPSRC:.cpp=.o) $(PCSRC:.c=.o)  

LIBRARY	= 

# files for the program
EXEC	= $(PACKAGE)
EXECSRC	= 
EXECOBJ	= $(POBJ)

# mathematica stuff (if it exists)
MMA_VERSION=5.1
MMA_MLINKDIR = /h/23/jhoey/Mathematica/${MMA_VERSION}/AddOns/MathLink/DeveloperKit
MMA_SYS = Linux
MMA_CADDSDIR = ${MMA_MLINKDIR}/${MMA_SYS}/CompilerAdditions
MMA_INCDIR = ${MMA_CADDSDIR}
MMA_LIBDIR = ${MMA_CADDSDIR}

FLTKDIR = /usr/local/src/fltk-1.1.4
#FLTKDIR = /lci/project/little/jhoey/software/fltk-1.1.2

IFLAGS1 =  -I$(FLTKDIR) -I$(MMA_INCDIR)
LFLAGS1 = -L/usr/X11R6/lib  -L$(FLTKDIR)/lib/
LIBS1 = -lfltk -lXext -lX11 -lXmu  -lXt -lXi 
LIBS2 = -lfltk -lXext -lX11 -lXmu  -lXt -lXi 
MMA_LFLAGS = -L$(MMA_LIBDIR)
MMA_LIBS = -lML -lrt
#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all:  $(EXEC)
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
installbin: $(EXEC) 
	install -p --mode=555 $(EXEC) $(BINDIR)

$(LIBRARY): $(POBJ)
	ar rv $@ $^
	$(RANLIB) $@

# file suffixes to pay attention to
#.SUFFIXES: .o .cpp

# rules to build .o's from .c's
# note: use CXX since now includes .hh files
%.o: %.c $(PHDR)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXEC): $(EXECOBJ) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(MMA_LFLAGS) -o $@ $(EXECOBJ) $(LIBS) $(MMA_LIBS) $(MNEMLIB) -lm




distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx
	@\rm -f $(BINDIR)/$(EXEC) 
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))

# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)

# specify individual targets
mdpparser.tab.c:mdpparser.y 
	bison -d -v mdpparser.y
lex.yy.c:scanner.h scanner.l mdpparser.tab.h
	flex scanner.l
lex.yy.o:lex.yy.c scanner.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c lex.yy.c
mdpparser.tab.o:mdpparser.tab.c lex.yy.c 
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c mdpparser.tab.c
gspudd.o:gspudd.cpp 
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c gspudd.cpp
MDP.o: MDP.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c MDP.cpp
POMDP.o: POMDP.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c POMDP.cpp
mvdumpdot.o: mvdumpdot.c
	 $(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c mvdumpdot.c














