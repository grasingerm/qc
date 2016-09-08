## Makefile for triangle server and client

ACLOCAL_AMFLAGS  = -I ../../system

INCLUDES = -I${srcdir}/../include -I${srcdir}/../../include


## Client triangulation library

noinst_LTLIBRARIES = libtriang.la
libtriang_la_SOURCES = triangulation_client.c util.c


## Server program

bin_PROGRAMS = triang_server
triang_server_SOURCES = triangulation_server.c error.f
triang_server_LDADD = ../geompack/deltr3d/libdeltr3d.la \
../geompack/basic3d/libbasic3d.la ../geompack/general/libgeneral.la \
libtriang.la

