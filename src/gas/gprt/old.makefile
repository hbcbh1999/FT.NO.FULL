#
#		Makefile for program gprt:
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

lib = gas/gprt

sources =                \
	  gdiagnostic.c  \
	  gintext.c	 \
	  glayeravg.c	 \
          g2dprint.c     \
          gdriverstat.c  \
          gintstat.c     \
          glpdiff.c      \
          gprcur.c       \
          gprint.c       \
          gprstate.c     \
          grectstat.c    \
          grmdata.c

includes =               \
	   glayer.h      \
	   gprtprotos.h

include $(GAS)/make-gas.defs
