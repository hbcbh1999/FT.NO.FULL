#
#		Makefile for the HYP library:
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#

lib = hyp
ifndef UTIL
 UTIL := $(shell pwd | sed -e "s/\/$(lib).*//")/util
endif #ifndef UTIL

sources = hdriver.c  \
          hinit.c    \
          hsub.c     \
          hsrc.c     \
          hscatter.c \
          hwave.c    \
          hnpt.c     \
          hvec.c     \
          hprint.c   \
	  hpseudo.c  \
          hsoln.c    \
          hoverture_driver.c

includes = hdecs.h   \
           hprotos.h

include $(UTIL)/make.defs
