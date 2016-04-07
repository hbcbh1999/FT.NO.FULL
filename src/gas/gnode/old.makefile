#
#		Makefile for program libgnode:
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

lib = gas/gnode

sources =             \
          gRP.c       \
          ganom.c     \
          gbnode.c    \
          gccnode.c   \
	  gcurve.c    \
          gmdnode.c   \
	  gnode.c     \
          gnodesub.c  \
          gpcnode.c   \
          gscnode.c   \
          gsndnode.c  \
          gscnsts.c   \
          gssnode.c   \
          gssnsts.c

includes = gnodeprotos.h

include $(GAS)/make-gas.defs
