# $Id: Makefile.blastdb_format_unit_test.app 186851 2010-03-24 18:48:36Z ucko $

APP = blastdb_format_unit_test
SRC = seq_writer_unit_test

CPPFLAGS = $(ORIG_CPPFLAGS) $(BOOST_INCLUDE)
CXXFLAGS = $(FAST_CXXFLAGS)
LDFLAGS = $(FAST_LDFLAGS)

LIB_ = test_boost blastdb_format xobjutil seqdb blastdb $(SOBJMGR_LIBS)

LIB = $(LIB_:%=%$(STATIC))
LIBS = $(DL_LIBS) $(ORIG_LIBS)

CHECK_REQUIRES = in-house-resources
CHECK_CMD = blastdb_format_unit_test
CHECK_COPY = data

REQUIRES = Boost.Test.Included


WATCHERS = blastsoft
