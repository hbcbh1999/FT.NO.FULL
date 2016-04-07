#
#		Makefile for the gprop:
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

lib = gas/gprop

sources = griemann.c        \
          griemcombst.c     \
          gpolar.c          \
          gipolar.c         \
          gprop.c           \
          gwspeed.c

includes = gpropprotos.h

include $(GAS)/make-gas.defs
