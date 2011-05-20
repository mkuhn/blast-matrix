# $Id: Makefile.blastdbcheck.app 208174 2010-10-14 16:52:22Z serovav $

WATCHERS = camacho maning 

APP = blastdbcheck
SRC = blastdbcheck blastdb_aux

LIB_ = $(BLAST_INPUT_LIBS) $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = $(LIB_:%=%$(STATIC))

CFLAGS   = $(FAST_CFLAGS)
CXXFLAGS = $(FAST_CXXFLAGS)
LDFLAGS  = $(FAST_LDFLAGS)

CPPFLAGS = $(ORIG_CPPFLAGS)
LIBS = $(CMPRS_LIBS) $(DL_LIBS) $(NETWORK_LIBS) $(ORIG_LIBS)

REQUIRES = objects -Cygwin
