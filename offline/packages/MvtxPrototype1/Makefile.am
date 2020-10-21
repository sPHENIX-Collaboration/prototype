AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libMvtxPrototype1.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

libMvtxPrototype1_la_LIBADD = \
  -lCLHEP \
  -lphool \
  -lfun4all \
  -lphparameter \
  -ltrack_io \
  -ltrack_util \
  -lmvtx

pkginclude_HEADERS = \
  MvtxAlign.h \
  MvtxApplyHotDead.h \
  MvtxClusterizer.h \
  MvtxRunInfoUnpackPRDF.h \
  MvtxStandaloneTracking.h \
  MvtxUnpackPRDF.h \
  MvtxQAHisto.h 

libMvtxPrototype1_la_SOURCES = \
  MvtxAlign.C \
  MvtxApplyHotDead.C \
  MvtxClusterizer.C \
  MvtxRunInfoUnpackPRDF.C \
  MvtxStandaloneTracking.C \
  MvtxUnpackPRDF.C \
  MvtxQAHisto.C

################################################
# linking tests

BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_tracker

testexternals_tracker_SOURCES = testexternals.C
testexternals_tracker_LDADD = libMvtxPrototype1.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f $(BUILT_SOURCES)
