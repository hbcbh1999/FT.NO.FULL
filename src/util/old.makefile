# -*-Makefile-*-
#		Makefile for the UTIL library:
#		Requires gnu make version 3.70 or higher
#
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#


#
#               Specification of variable storage precision
# 
#       Default equals single precision
# 


ifndef UTIL
 UTIL = .
endif #ifndef UTIL
lib = util

ifeq "$(PRECISION)" "double"
 ads = d$(T_ARCH)
endif #ifeq "$(PRECISION)" "double"

sources =		\
          cleanup.c     \
          debug.c       \
          error.c       \
          fgetstrin.c   \
          machine.c	\
          matrix.c      \
	  navigate.c	\
          other.c       \
	  fsort.c       \
          output.c      \
          ppsub.c       \
          quad.c        \
          roots.c       \
          runga.c       \
	  sphhar.c	\
          screen.c      \
          sigplot.c     \
          simpleio.c    \
          times.c       \
          uinit.c       \
          vectormalloc.c

includes =		\
	   cdecs.h      \
           fnamedebug.h \
	   navdecs.h	\
           plotdecs.h   \
           vmalloc.h    \
           uprotos.h

libutil:

all: libutil libutil_p

morefiles = make.defs make.master $(wildcard test*.c)

include $(UTIL)/make.defs

ifeq "$(TRACE_BACK)" "yes"
 CFLAGS += -D__TRACE_BACK__
endif #ifeq "$(TRACE_BACK)" "yes"

libutil:
	-$(MAKE) $(JFLAG) "useCC=no" "ads=d$(T_ARCH)" "PRECISION=double" libdutil$(SS).$(LIB_SUFFIX)
	-$(MAKE) $(JFLAG) "useCC=no" "PRECISION=" libutil$(SS).$(LIB_SUFFIX)
ifeq "$(PP)" ""
 ifneq "$(findstring IRIX,$(unix))" ""
	-$(MAKE) $(JFLAG) "useCC=no" "sgiABIopts=-o32" "PRECISION=" libutil$(SS).$(LIB_SUFFIX)
 endif #ifneq "$(findstring IRIX,$(unix))" ""
 ifneq "$(findstring v9,$(sunABIopts))" ""
	-$(MAKE) $(JFLAG) "useCC=no" "sunABIopts=" "PRECISION=" libutil$(SS).$(LIB_SUFFIX)
 endif #ifneq "$(findstring v9,$(sunABIopts))" ""
endif #ifeq "$(PP)" ""
	-$(MAKE) $(JFLAG) "useCC=yes" "ads=d$(T_ARCH)-CC" "PRECISION=double" libdutil$(SS)-CC.$(LIB_SUFFIX)
	-$(MAKE) $(JFLAG) "useCC=yes" "PRECISION=" libutil$(SS)-CC.$(LIB_SUFFIX)
ifndef NOMPI
ifdef PPVERSION
	-$(MAKE) $(JFLAG) "useCC=no" "ads=d$(T_ARCH)" "PRECISION=double" "PP=$(PPVERSION)" libdutil$(SS).$(LIB_SUFFIX)
	-$(MAKE) $(JFLAG) "useCC=no" "PRECISION=" "PP=$(PPVERSION)" libutil$(SS).$(LIB_SUFFIX)
	-$(MAKE) $(JFLAG) "useCC=yes" "ads=d$(T_ARCH)-CC" "PRECISION=double" "PP=$(PPVERSION)" libdutil$(SS)-CC.$(LIB_SUFFIX)
	-$(MAKE) $(JFLAG) "useCC=yes" "PRECISION=" "PP=$(PPVERSION)" libutil$(SS)-CC.$(LIB_SUFFIX)
endif #ifdef PPVERSION
endif #ifndef NOMPI

libutil_p:
ifndef noprofile
	-$(MAKE) "PROF=-p" libutil
endif #ifndef noprofile

ifeq "$(MAKELEVEL)" "0"
dirpath = ..
VPATH = .
endif #ifeq "$(MAKELEVEL)" "0"

USEROPTIONS =

testfiles = testdebug testoutpu testsigplot testscree testtimes \
	    testvmall testfgets testcleanup    testpp   testnav \
	    testquad  testroots

USERREALCLEAN = $(testfiles)

test: $(testfiles)

lintall:
	-$(MAKE) "PRECISION=" llib-lutil.ln
	-$(MAKE) "PRECISION=double" llib-ldutil.ln

interfaces = vmalloc_i.c

interfaces: $(foreach var,$(interfaces),../$(hosttype)-insight/$(var:.c=.tqs))

../$(hosttype)-insight/%.tqs: %.c
	iic -I$(UTIL) $< -o $@

clean: clean-util

clean-util:
	@$(RM) */*.d

testmach: testmach.o fepsilon.o

.PHONY: libutil libutil_p
