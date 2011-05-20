# $Id: Makefile.unit_test_mol_wt.app 194095 2010-06-10 13:47:39Z dicuccio $

APP = unit_test_mol_wt
SRC = unit_test_mol_wt

CPPFLAGS = $(ORIG_CPPFLAGS) $(BOOST_INCLUDE)

LIB = test_boost xobjutil $(SOBJMGR_LIBS)
LIBS = $(NETWORK_LIBS) $(DL_LIBS) $(ORIG_LIBS)

REQUIRES = Boost.Test.Included

CHECK_CMD =

WATCHERS = ucko dicuccio
