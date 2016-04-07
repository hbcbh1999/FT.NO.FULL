#
#		Makefile for program libgstate:
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

lib = gas/gstate

sources = gintrp.c         \
          gstate.c         \
          gbstate.c        \
          grstate.c        \
          gstglobs.c       \
          gparams.c        \
          goverinterp.c

includes = gstateprotos.h

include $(GAS)/make-gas.defs
