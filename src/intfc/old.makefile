#		Makefile for the INTFC library:
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#

lib = intfc
liblibs = $(util)
ifndef UTIL
 UTIL := $(shell pwd | sed -e "s/\/$(lib).*//")/util
endif #ifndef UTIL

sources =                 \
	  arrayutils.c    \
	  iblkb.c	  \
	  iblkc2.c        \
	  imkcurve.c      \
	  imksurf.c 	  \
	  igview.c        \
	  ixgraph.c       \
          comp.c          \
          comp1d.c        \
          comp2d.c        \
          comp3d.c        \
          cross2d.c       \
          geomutils.c     \
          icheck3d.c      \
	  idiagnostic.c   \
          iecomps.c       \
          ifourier.c      \
          igrid.c         \
          int3d.c         \
          intfc.c         \
          iprt3d.c        \
          irefl.c         \
          iredist.c       \
          iscatter.c      \
          isect2d.c       \
          isect3d.c       \
          isub.c          \
          iuserintfc.c    \
          map.c           \
          ppcopy.c        \
          setb1d.c        \
          setb2d.c        \
          setb3d.c        \
          shift.c         \
          top.c           \
	  triangle.c      \
          trisurf.c       \
          userhooks.c	  \
          zoom.c          \
          intfc_amr.c

includes =                \
	   array.h        \
           geom.h         \
	   iloc.h         \
           ilocprotos.h   \
           int.h          \
           iprotos.h      \
           table.h        \
	   triangledefs.h \
           userint.h      \
           int_amr.h

#
# Instructions for making testintfc
#

ifndef PRECISION
 ifneq "$(unix)" "UNICOS"
  PRECISION = double
 else #ifneq "$(unix)" "UNICOS"
  PRECISION = single
 endif #ifneq "$(unix)" "UNICOS"
 export PRECISION
endif #ifndef PRECISION

ifndef DIMENSION
 DIMENSION = -DONED -DTWOD -DTHREED
 export DIMENSION
endif #ifndef DIMENSION

include $(UTIL)/make.defs

LOADLIBES := -L$(liblibdir) -$(LNK)intfc$(CCtag) $(LOADLIBES)

testintfc.o testtriangle.o: libutil libintfc$(CCtag)

libutil:
	-@$(CD) ../util; $(MAKE) libutil

ifeq "$(shared)" "yes"
 libintfc$(CCtag): libintfc$(CCtag).so
else #ifeq "$(shared)" "yes"
 libintfc$(CCtag): libintfc$(CCtag).a
endif #ifeq "$(shared)" "yes"

.PHONY: libutil libintfc
