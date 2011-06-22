################################################################
#
#	SERIAL SGI Makefile for GA.
#
################################################################

SRC_DIR		= $(OPT_ROOT)/libopt/GA

include $(OPT_ROOT)/util/$(OPT_ARCH)opt++.defs

CFLAGS          = -c

C++		= $(OPT_COMP)
C++FLAGS        = $(OPT_FLAGS)
C++DEFINES      = $(OPT_DEFINES) -D$(OPT_ARCH) -D$(PLL_TARGET)
C++INCLUDES     = $(OPT_INCLUDES) -I$(OPT_ROOT)/libopt -I$(OPT_ROOT)/libfcn \
		  -I$(OPT_ROOT)/newmat7
C++HEADERS	= $(OPT_HEADERS)

LIBOPT		= $(OPT_ROOT)/libopt/gen/$(OPT_ARCH)/libopt.a

GA_HEADERS	= $(SRC_DIR)/pga.h $(SRC_DIR)/optga.h
GA_OBJS     	= ga.o pga_util.o

all:		$(LIBOPT)
	@echo "==== > SERIAL GA make complete"

$(LIBOPT):	$(GA_OBJS)
	ar rcv $(LIBOPT) $(GA_OBJS) 

clean:
	rm -f core *.o
	@-rm -f *~

clobber:	clean

delete: 	clobber

ga.o:		$(C++HEADERS) $(GA_HEADERS) $(SRC_DIR)/ga.C
	$(C++) $(C++FLAGS) $(C++DEFINES) $(C++INCLUDES) $(SRC_DIR)/ga.C

pga_util.o : $(GA_HEADERS) $(SRC_DIR)/pga_util.c
	$(CC) $(CFLAGS) $(CDEFINES) $(CINCLUDES) $(SRC_DIR)/pga_util.c


