#
#		Makefile for program ghyp:
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

lib = gas/ghyp

sources = gflw.F        \
          ggodunov.c    \
	  ghyperbolic.c \
          ghypvec.c     \
          glf.c         \
          glw.c         \
          gplm.c        \
          gmuscl.c      \
          gcgrsolve.c   \
          grsolve.c     \
	  ghypprt.c     \
	  ghypsub.c     \
	  gmoc.c	\
	  gvisc.c	

includes = ghyp.h       \
           ghypprotos.h

include $(GAS)/make-gas.defs
