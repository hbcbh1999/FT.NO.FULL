#
#		Makefile for the directory ginit:
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

lib = gas/ginit

sources =                \
	  gibifur.c      \
	  gijet.c        \
	  giphysprompt.c \
          gicc.c         \
          gielrp.c       \
          giglobs.c      \
          gihypinit.c    \
          gilayer.c      \
          gictype.c      \
          gimkcur.c      \
          gimksurf.c     \
          ginitintfc.c   \
          ginitphys.c    \
          giparams.c     \
          gipert.c       \
          girt.c         \
          girpregion.c   \
          gipppert.c     \
          giprt.c        \
          gireadstate.c  \
          girefl.c       \
          girestrt.c     \
          girm_linear.c  \
          girstate.c     \
          gisc.c         \
          gistate.c      \
	  gitabregion.c  \
	  giboone.c	 \
          spolars.c	 \
	  testsolver.c

includes =               \
           gicomptype.h  \
           ginit.h       \
	   ginitprotos.h

include $(GAS)/make-gas.defs
