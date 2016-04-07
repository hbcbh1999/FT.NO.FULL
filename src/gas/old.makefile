#
#		Makefile for the program GAS:
#
#	Copyright 1999 by The University at Stony Brook, All rights reserved.
#
#	Remark about undefined functions:
#	
#	Due to a bug in most unix loaders,  it is possible (and likely)
#	that the linker will falsely identify several functions
#	as undefined.  The -u option to ld is designed to get around
#	this problem by forcing the loader to look for the undefined
#	functions.  Since the functions that are not found can include
#	system functions (such as in libF77.a),  this list may be
#	machine dependent.  For each executable produced by this
#	makefile,  there is a an executable file in this directory 
#	called undef/[executable name].  This file simply echos
#	the list of undefined functions.  This list can be updated
#	as needed.
#

ifndef GAS
 GAS = .
endif #ifndef GAS

lib = gas

sources = $(lib)-main.c

include $(GAS)/make-gas.libs
include $(GAS)/make-gas.defs

#
#	Additional Executable Programs
#

spolars eosplot testsolver: gas
	cd $(STORE); if [ -h $@ -o -f $@ ] ; then $(RM) $@; fi
	cd $(STORE); ln -s gas $@

install.spolars install.eosplot install.testsolver: install
	cd $(FronTier)/bin/$(hosttype); \
	if [ -h $@ -o -f $(subst .,,$(suffix $@)) ] ; \
		then $(RM) $(subst .,,$(suffix $@)); \
	fi; \
	ln -s gas $(subst .,,$(suffix $@))


rcp_shock:
	make "RHOST=shock@$(hostname)"			\
		"RPATH=/nfs/n/shock/shock/libs" rcp_remote
	-make -i "DIFFER=$(DIFFER)" "DIFFNAME=shock" \
		"DIFFPATH=/nfs/n/shock/shock/libs" diff_local

rcp_remote:
	-make "host=$(RHOST)" "sitepath=$(RPATH)" rcp_all
	-cd ../util; make "host=$(RHOST)" "sitepath=$(RPATH)" rcp

diff_remote:
	-@make -i "DIFFER=$(DIFFER)" "DIFFNAME=$(RHOST)"  "DIFFPATH=$(RHOST):$(RPATH)" differ_all > DIFF-$(RHOST).log 2>&1
	-@cd ../util; make -i "DIFFER=$(DIFFER)" "DIFFNAME=$(RHOST)"  "DIFFPATH=$(RHOST):$(RPATH)" differ > DIFF-$(RHOST).log 2>&1

diff_local:
	-@make -i "DIFFER=$(DIFFER)" "DIFFNAME=$(DIFFNAME)" \
		"DIFFPATH=$(DIFFPATH)" differ_all > DIFF-$(DIFFNAME).log 2>&1
	-@cd ../util; make -i "DIFFER=$(DIFFER)" "DIFFNAME=$(DIFFNAME)" \
		"DIFFPATH=$(DIFFPATH)" differ > DIFF-$(DIFFNAME).log 2>&1

.PHONY: $(foreach var, $(rhosts), update.$(var) xupdate.$(var))

libintfc: $(mylibdir)/libintfc$(CCtag).$(LIB_SUFFIX)
libfront: $(mylibdir)/libfront$(CCtag).$(LIB_SUFFIX)
libtri: $(mylibdir)/libtri$(CCtag).$(LIB_SUFFIX)
libhyp: $(mylibdir)/libhyp$(CCtag).$(LIB_SUFFIX)
libdriver: $(mylibdir)/libdriver$(CCtag).$(LIB_SUFFIX)
libginit: $(mylibdir)/gas/libginit$(CCtag).$(LIB_SUFFIX)
libgprt: $(mylibdir)/gas/libgprt$(CCtag).$(LIB_SUFFIX)
libgbifur: $(mylibdir)/gas/libgbifur$(CCtag).$(LIB_SUFFIX)
libgnode: $(mylibdir)/gas/libgnode$(CCtag).$(LIB_SUFFIX)
libghyp: $(mylibdir)/gas/libghyp$(CCtag).$(LIB_SUFFIX)
libgprop: $(mylibdir)/gas/libgprop$(CCtag).$(LIB_SUFFIX)
libgstate: $(mylibdir)/gas/libgstate$(CCtag).$(LIB_SUFFIX)
libgeos: $(mylibdir)/gas/libgeos$(CCtag).$(LIB_SUFFIX)
libgintfc: $(mylibdir)/gas/libgintfc$(CCtag).$(LIB_SUFFIX)

.PHONY: libintfc libfront libtri libhyp libdriver
.PHONY: ginit gprt gbifur gnode ghyp gprop gstate geos gintfc

#ifeq "$(useCC)" "no"
#dummy := $(shell /bin/rm ../intfc/$(T_ARCH)/*.d)
#endif #ifeq "$(useCC)" "no"
