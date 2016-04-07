#
#		Makefile for program libgintfc:
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

lib = gas/gintfc

sources = gcheck3d.c       \
	  gintfcpert.c     \
          gtop.c           \
          gtypes.c         \
          guserhooks.c     \
          guserintfc.c

includes = gintfcprotos.h

include $(GAS)/make-gas.defs
