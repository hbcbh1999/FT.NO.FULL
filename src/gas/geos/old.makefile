#
#		Makefile for the geos:
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#

ifndef GAS
 ifdef absdirpath
   GAS = $(absdirpath)/gas
 else #ifdef absdirpath
   GAS = ..
 endif #ifdef absdirpath
endif #ifndef GAS

lib = gas/geos

sources = eosplot.c     \
	  generic-eos.c \
	  geosvec.c   \
	  giniteos.c    \
	  gphriem.c     \
	  gseshyp.c     \
	  gsesintrp.c   \
	  gsesphase.c   \
	  gsesspline.c  \
	  gsesprint.c   \
	  gsesinout.c   \
	  mpoly-eos.c   \
	  poly-eos.c    \
	  sesame-eos.c  \
	  spoly-eos.c   \
	  gpertsub.F    \
	  gsestoc.F     \
	  sesadd.F      \
	  sesame.F      \
	  sesinv.c      \
	  sesspln.c     \
	  sesstate.c    \
	  jwl-eos.c     \
	  mg-eos.c      \
	  gentest-eos.c \
	  s2phase-eos.c

includes = geosdecs.h   \
	   geosprotos.h \
	   gsesprotos.h \
	   mpoly.h      \
	   poly.h       \
	   sesame.h     \
	   spoly.h      \
	   jwl.h        \
	   mg.h         \
	   gentest.h	\
	   s2phase.h

morefiles = mkseslib.c testseslib.c

include $(GAS)/make-gas.defs

ifeq "$(findstring -DTWOD,$(CPPFLAGS))" "-DTWOD"
 CPPFLAGS += -DSESAME_CODE $(MKSESLIB)
 ifeq "$(findstring -DPHASE_CODE,$(CPPFLAGS))" "-DPHASE_CODE"
  ifneq "$(findstring -DSESAME_CODE,$(CPPFLAGS))" "-DSESAME_CODE"
   CPPFLAGS += -DSESAME_CODE
  endif #ifneq "$(findstring -DSESAME_CODE,$(CPPFLAGS))" "-DSESAME_CODE"
 endif #ifeq "$(findstring -DPHASE_CODE,$(CPPFLAGS))" "-DPHASE_CODE"
endif #ifeq "$(findstring -DTWOD,$(CPPFLAGS))" "-DTWOD"

ifeq "$(MKSESLIB)" "-DMAKE_SES_LIB"
 mkseslib: sesame.o gsestoc.o mkseslib.o
	-$(LD) -o mkseslib $^ $(LDFLAGS) $(LOADLIBES)
	-@$(RM) *.o
else #ifeq "$(MKSESLIB)" "-DMAKE_SES_LIB"
 mkseslib: sesame.F gsestoc.F mkseslib.c
	$(MAKE) "MKSESLIB=-DMAKE_SES_LIB" mkseslib
endif #ifeq "$(MKSESLIB)" "-DMAKE_SES_LIB"

testseslib: sesame.o gsestoc.o testseslib.o
	-$(LD) -o testseslib $^ $(LDFLAGS) $(LOADLIBES)
	-@$(RM) *.o
