#
#		Makefile for the FRONT library:
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#

lib = front
ifndef UTIL
 UTIL := $(shell pwd | sed -e "s/\/$(lib).*//")/util
endif #ifndef UTIL


sources = fadv.c       \
          fbdry1.c     \
          fbdry2.c     \
          fbdry3.c     \
          fbdry4.c     \
	  fcheck3d.c   \
          fcorrspnd.c  \
          fcrosscur.c  \
          fcrossext.c  \
          fcrstatus.c  \
	  fdiagnostic.c\
          finit.c      \
          fint.c       \
          fnode.c      \
          fnodesub.c   \
          fprint.c     \
          fprop2d.c    \
          fprop3d.c    \
          fredist.c    \
          fredist3d.c  \
          frp1.c       \
          frp2.c       \
          fscatter.c   \
          fredist1d.c  \
          fredist2d.c  \
          fscat1d.c    \
          fscat2d.c    \
          fscat3d1.c   \
          fscat3d2.c   \
          fstate2d.c   \
          fstate.c     \
          ftop.c       \
          fsub.c       \
          funtan2d.c   \
          funtan3d.c   \
          fuserintfc.c \
          fuserhooks.c \
          foverture_patch.c \
          foverture_adv.c

includes = fdecs.h     \
           fprotos.h   \
           frp.h       \
           fuserint.h


#
# Instructions for making testfront (excutable)
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

LOADLIBES := -L$(liblibdir) -$(LNK)front$(CCtag) -$(LNK)intfc$(CCtag) -$(LNK)tri$(CCtag) $(LOADLIBES) -lm -llinpak

testfront.o : libutil libintfc$(CCtag) libtri$(CCtag) libfront$(CCtag)

libutil:
	-@$(CD) ../util; $(MAKE) libutil

libintfc$(CCtag): 
	-@$(CD) ../intfc; $(MAKE) libintfc$(CCtag)

libfront$(CCtag): libfront$(CCtag).a

libtri$(CCtag): 
	-@$(CD) ../tri; $(MAKE) libtri$(CCtag)


.PHONY: libutil libintfc libfront
