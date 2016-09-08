#
# $Id$
#
# $Log$
#
AUTOMAKE_OPTIONS = 1.7
ACLOCAL_AMFLAGS  = -I ../../system

INCLUDES = -I${srcdir}/../../include -I${srcdir}/../../triangulation_server/include -I/usr/local/include

#LDFLAGS = 
#
#
#

bin_PROGRAMS = quasi test

quasi_SOURCES = main.cc stamp.cc
quasi_LDADD = ../../src/libquasi.la ../../triangulation_server/src/libtriang.la -L/usr/local/lib -lCGAL -lCGAL_Core -lgmp -lboost_thread -lmpfr @DMALLOC_C_LIBS@

test_SOURCES = test.cc stamp.cc
test_LDADD = ../../src/libquasi.la ../../triangulation_server/src/libtriang.la -L/usr/local/lib -lCGAL -lCGAL_Core -lgmp -lboost_thread -lmpfr @DMALLOC_C_LIBS@


all-local: triang_server

triang_server:
	$(RM) triang_server
	$(LN_S) ../../triangulation_server/src/triang_server .

clean-local:
	$(RM) triang_server
