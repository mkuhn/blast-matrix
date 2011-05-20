# $Id: Makefile.unit_test_seq_translator.app 194094 2010-06-10 13:46:59Z dicuccio $

APP = unit_test_seq_translator
SRC = unit_test_seq_translator

CPPFLAGS = $(ORIG_CPPFLAGS) $(BOOST_INCLUDE)

LIB = test_boost xobjutil $(SOBJMGR_LIBS)
LIBS = $(NETWORK_LIBS) $(DL_LIBS) $(ORIG_LIBS)

REQUIRES = Boost.Test.Included

CHECK_CMD =

WATCHERS = vasilche dicuccio bollin
