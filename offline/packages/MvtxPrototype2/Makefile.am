AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmvtxprototype2.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

libmvtxprototype2_la_LIBADD = \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lmvtx

pkginclude_HEADERS = \
  MvtxPrototype2Align.h \
  MvtxPrototype2Geom.h \
  MvtxPrototype2UnpackPRDF.h \
  MvtxPrototype2Clusterizer.h \
  MvtxStandaloneTracking.h

libmvtxprototype2_la_SOURCES = \
  MvtxPrototype2Align.cc \
  MvtxPrototype2Geom.cc \
  MvtxPrototype2UnpackPRDF.cc \
  MvtxPrototype2Clusterizer.cc \
  MvtxStandaloneTracking.cc

################################################
# linking tests

BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_tracker

testexternals_tracker_SOURCES = testexternals.C
testexternals_tracker_LDADD = libmvtxprototype2.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f $(BUILT_SOURCES)
