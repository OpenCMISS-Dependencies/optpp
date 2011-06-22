##############################################################################
#
#	Makefile for opt++ Examples
#
##############################################################################
#
# What Compiler do we want to use
#
#
OPT_ROOT	= ../../..
OPT_ARCH	= SGI
#OPT_ARCH	= HPPA
#OPT_ARCH	= SUN4
#OPT_ARCH	= PGON

include $(OPT_ROOT)/util/$(OPT_ARCH)opt++.defs

SRC_DIR		= $(OPT_ROOT)

CFLAGS          = -c -O2 -mips2
CINCLUDES	= -I../

C++		= $(OPT_COMP)
C++FLAGS	= $(OPT_FLAGS) -g
C++DEFINES      = $(OPT_DEFINES) -D$(OPT_ARCH)
C++INCLUDES	= $(OPT_INCLUDES) -I$(OPT_ROOT)/libopt -I$(OPT_ROOT)/libfcn -I$(OPT_ROOT)/newmat7 -I$(OPT_ROOT)/libopt/GA

LIBOPT		= $(OPT_ROOT)/libopt/gen/$(OPT_ARCH)/libopt.a 
LIBFCN		= $(OPT_ROOT)/libfcn/gen/$(OPT_ARCH)/libfcn.a 
LIBNEWMAT	= $(OPT_ROOT)/newmat7/gen/$(OPT_ARCH)/libnewmat7.a 
LIBPDS 		= $(OPT_ROOT)/libpds/gen/$(OPT_ARCH)/libpds.a 

OPT_LIBS	= $(LIBOPT) $(LIBFCN) $(LIBNEWMAT)

LIBS		= $(OPT_LIBS) -lfastm -lm -lfpe

#  End of user defined flags
##############################################################################

PGMS		= tstga 

all:		$(PGMS)
	@echo "==== > OPT++ Examples made"

GA_OBJS		= eval.o polymer.o dnumrec.o vector.o

verify:	$(PGMS)
	make all
	@verify
#
#  Main test routine
#
tstga:	tstga.o rosen.o $(GA_OBJS)
	$(C++) -o tstga tstga.o rosen.o $(GA_OBJS) $(LIBS) 

tstga.o: tstga.C
	$(C++) $(C++FLAGS) $(C++DEFINES) $(C++INCLUDES) tstga.C 

eval.o :  $(GA_HEADERS) eval.c
	$(CC) $(CFLAGS) $(CCDEFINES) $(CINCLUDES) eval.c

polymer.o :  $(GA_HEADERS) polymer.c
	$(CC) $(CFLAGS) $(CCDEFINES) $(CINCLUDES) polymer.c

vector.o :  $(GA_HEADERS) vector.c
	$(CC) $(CFLAGS) $(CCDEFINES) $(CINCLUDES) vector.c

dnumrec.o :  $(GA_HEADERS) dnumrec.c
	$(CC) $(CFLAGS) $(CCDEFINES) $(CINCLUDES) dnumrec.c

rosen.o: rosen.C
	$(C++) $(C++FLAGS) $(C++DEFINES) $(C++INCLUDES) rosen.C 
clean:
	/bin/rm -f core *.o *.a *.d GAout.dat GApop*.dat GAbest*.dat GAarray*.dat

clobber:
	/bin/rm -f core *.o *.a *.d *.out $(PGMS) SCHEME
