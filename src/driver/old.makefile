#
#		Makefile for the DRIVER package:
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#

ifndef HDF
 HDF = yes
 export HDF
endif #ifndef HDF
lib = driver
ifndef UTIL
 UTIL := $(shell pwd | sed -e "s/\/$(lib).*//")/util
endif #ifndef UTIL

sources = dmain.c    \
          dinit.c    \
          diprt.c    \
          dprint.c   \
          dsub.c     \
          dinout.c   \
          dstat.c    \
	  dscatter.c \
          dpatchmesh.c \
	  doverturepatch.c \
          doverturepatch2.c

includes = damr.h    \
           ddecs.h   \
           dprt.h    \
           dprotos.h

include $(UTIL)/make.defs
