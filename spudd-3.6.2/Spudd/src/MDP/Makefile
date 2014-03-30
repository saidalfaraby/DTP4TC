# $Id: Makefile,v 1.31 2005/02/14 14:43:18 jhoey Exp $
#
#    Spudd - Spudd package
#---------------------------

# Heavily modified by Don Murray, Sept 2000

LIBS	= -lterminal -lobj -lcudd -lmtr -lst -lutil -ldddmp -lm

# files for the package
LIBRARY	= libspudd.a

PHDR = MDP.h POMDP.h net2.h networkSettings.h

LIBSRC	=  MDP.cc POMDP.cc FMDP.cc net2.cc
LIBCSRC = dumpdot.c mvdumpdot.c mdpparser.tab.c mdp_parser.c 
LIBOBJ 	= $(LIBSRC:%.cc=%.o) $(LIBCSRC:%.c=%.o) 

# files for the program
EXEC3	= Spudd
EXEC3SRC = 
EXEC3OBJ = testmdp.o

EXEC4	= testmdpserver
EXEC4SRC = 
EXEC4OBJ = testmdpserver.o

EXEC5	= pquery
EXEC5OBJ = policyQuery.o policyQueryMain.o

EXEC6	= mdpquery
EXEC6OBJ = mdpquery.o

EXEC7	= pcomp
EXEC7OBJ = policyCompare.o

EXEC8	= hweval
EXEC8OBJ = testhw.o


#FLTKDIR = /usr/local/lib/fltk-1.1.2
FLTKDIR = /lci/project/little/jhoey/software/fltk-1.1.2

IFLAGS1 =  -I./
IFLAGS2 =   -I./ -I$(FLTKDIR)
IFLAGS3 =   -I./ -I$(FLTKDIR)
LFLAGS1 = -L/usr/X11R6/lib 
LFLAGS2 = -L/usr/X11R6/lib  -L$(FLTKDIR)/lib/
LIBS1 = -lXext -lX11 -lXmu  -lXt -lXi 
LIBS2 = -lfltk -lXext -lX11 -lXmu  -lXt -lXi 

#------------------------------------------------------

# don't forget required targets 'installinclude' 'installlib' 'installbin'

all:   $(LIBRARY) $(EXEC3) $(EXEC4) $(EXEC7) $(EXEC8)
install: installinclude installlib installbin

# copy include files to the include directory
installinclude: 
	install -p --mode=444 $(PHDR) $(INCLUDEDIR)
# no includes

# build library and copy to lib directory
installlib: $(LIBRARY)
	install -p --mode=444 $(LIBRARY) $(LIBDIR)
# no library for this package

# build executables and copy to the bin directory
installbin: $(EXEC3) $(EXEC4) $(EXEC6)  $(EXEC7) $(EXEC8)
	install -p --mode=555 $(EXEC3) $(BINDIR)

# include the general rules files to customize this makefile
include ../Makefile.rules
include ../Makefile.rules.$(BUILD)

CUDDLIBS = $(CUDDLIBDIR)/lib*.a 


$(LIBRARY): $(LIBOBJ)  unpacklibs
	rm -f $(LIBRARY)
	ar r $(LIBRARY) $(LIBOBJ) tmp/*.o; 
	rm -fr tmp;

unpacklibs: 
	rm -fr tmp
	mkdir tmp
	for lib in $(CUDDLIBS) ; do \
	( cd tmp; \
	  ar x ../$$lib; \
          touch ../$$lib ); \
	done



# file suffixes to pay attention to
#.SUFFIXES: .o .c

# rules to build .o's from .c's
# note: use CXX since now includes .hh files
%.o: %.c $(PHDR)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXEC3): $(EXEC3OBJ)  $(LIBRARY) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS1) -L ./ -o $@ $(EXEC3OBJ) -lspudd
$(EXEC8): $(EXEC8OBJ)  $(LIBRARY) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS1) -L ./ -o $@ $(EXEC8OBJ) -lspudd
$(EXEC7): $(EXEC7OBJ)  $(LIBRARY) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS1) -L ./ -o $@ $(EXEC7OBJ) -lspudd
$(EXEC4): $(EXEC4OBJ) $(LIBRARY) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS1) -L ./ -o $@ $(EXEC4OBJ) -lspudd
$(EXEC5): $(EXEC5OBJ) $(LIBRARY) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS2) -L ./ -o $@ $(EXEC5OBJ) -lspudd $(LIBS2) 
$(EXEC6): $(EXEC6OBJ)  $(LIBRARY) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS1) -L ./ -o $@ $(EXEC6OBJ) -lspudd

$(EXECt): $(EXECtOBJ)  $(LIBRARY) $(PHDR) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LFLAGS1) -L ./ -o $@ $(EXECtOBJ) -lspudd

distclean: clean
	@echo "dist cleaning $(PACKAGE)"
	@\rm -f $(EXEC) $(LIBRARY) *.bak *~ .gdb_history *.qv *.qx 
	@\rm -f $(BINDIR)/$(EXEC) 
	@(cd $(INCLUDEDIR); \rm -f $(PHDR))

# specify individual targets
mdpparser.tab.c:mdpparser.y mdp-pspudd.h
	bison -d -v mdpparser.y
mdp_parser.o:mdp_parser.c 
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c mdp_parser.c
lex.yy.c:scanner.h scanner.l mdpparser.tab.h
	flex scanner.l
lex.yy.o:lex.yy.c scanner.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c lex.yy.c
mdpparser.tab.o:mdpparser.tab.c lex.yy.c mdp-pspudd.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c mdpparser.tab.c
testmdp.o:testmdp.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c testmdp.cc
testhw.o:testhw.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c testhw.cc
mdpquery.o:mdpquery.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c mdpquery.cpp
testmdpserver.o:testmdpserver.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c testmdpserver.cpp
mvdumpdot.o: mvdumpdot.c
	 $(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c mvdumpdot.c
dumpdot.o:dumpdot.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c dumpdot.c
MDP.o: MDP.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c MDP.cc
POMDP.o: POMDP.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c POMDP.cc
FMDP.o: FMDP.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c FMDP.cc
policyQuery.o: policyQuery.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c policyQuery.cpp
policyQueryMain.o: policyQueryMain.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c policyQueryMain.cpp
net2.o: net2.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c net2.cc
policyCompare.o: policyCompare.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(IFLAGS1) -c policyCompare.cc















