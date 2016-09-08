## Makefile.am for main QC3D code
##
## Notes:
## Build all binaries _HERE_ for now, and fix later
## Do not build itpack, ONLY triangulation_server and quasi
## Currently relies on higher-level makefile to build libtriang.a first


## Need all includes
AUTOMAKE_OPTIONS = 1.7
ACLOCAL_AMFLAGS  = -I ../system

INCLUDES = -I${srcdir}/../include -I${srcdir}/../triangulation_server/include -I$/usr/local/include


## main program: compile qc3d.c and link with all libraries

#bin_PROGRAMS = quasi
#quasi_SOURCES = qc3d.c stamp.cc
#quasi_LDADD = libquasi.la ../triangulation_server/src/libtriang.la  @DMALLOC_C_LIBS@


## Build library of main QC3D code

noinst_LTLIBRARIES = libquasi.la
libquasi_la_SOURCES = aluminum.f monitor.c numa.c range_search_generic.c read_pipe.c \
site_element_map.c site_element_map_cache.c threads.c \
SolverFunction.cc NonLinearSolver.cc Quasicontinua.cc Quasicontinuum.cc \
QuasicontinuaForceFunction.cc C_Interface.cc CGNonLinearSolver.cc \ CrossNeighborList.cc Electrostatics.cc Element.cc \
CreateMesh.cc ForceEnergyCalculation.cc Indent.cc Input.cc Lattice.cc \
MiscFunctions.cc Node.cc Output.cc PairPotentials.cc QuadraturePoints.cc \
Shape.cc Void.cc RunData.cc






