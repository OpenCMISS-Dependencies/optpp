#######################################################################
#
#	Top-level makefile for PDS.
#
#######################################################################

SHELL	= /bin/sh

all:
	cd gen; aimk all

install:
	cd gen; aimk install

clean:
	@-rm -f *~ gen/*~
	cd gen; aimk clean 

clobber:
	cd gen; aimk clobber 

delete:	
	@(cd gen; \
	    for DIR in * ;\
	    do \
		    if [ -d $$DIR ] ;\
		    then \
			    env OPT_ARCH=$$DIR aimk delete; \
		    fi; \
	    done \
	);
