#
#		Makefile for the TRI library:
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#

lib = tri
liblibs = $(util)
ifndef UTIL
 UTIL := $(shell pwd | sed -e "s/\/$(lib).*//")/util
endif #ifndef UTIL


sources =                 \
          tri1d.c         \
          tri2d.c         \
          tri3d.c         \
	  tri3dutils.c    \
          tricpy.c        \
          tricrx.c        \
          triel1.c        \
          triel2.c        \
          triel3.c        \
	  trigrid1.c      \
          trigrid2.c      \
          triloc.c        \
          tripcs.c        \
          triprint.c      \
	  trisurgery.c    \
	  triuserintfc.c  \
          overture_trigrid1.c 

includes =                \
	   trigrid.h      \
	   trilocaldecs.h \
	   tri3ddefs.h    \
           triprotos.h

include $(UTIL)/make.defs

ifeq "$(shared)" "yes"
 libtri$(CCtag): libtri$(CCtag).so
else #ifeq "$(shared)" "yes"
 libtri$(CCtag): libtri$(CCtag).a
endif #ifeq "$(shared)" "yes"
