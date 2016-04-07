#
#		Makefile for program gbifur:
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

lib = gas/gbifur

sources = gbifur.c    \
	  gcapture.c  \
          grefl.c     \
          grp.c       \
          gsc1.c      \
          gsc2.c      \
          gsc3.c      \
          guntan1d.c  \
          guntan2d.c  \
          gvecuntan.c

includes = gbifurprotos.h

include $(GAS)/make-gas.defs
