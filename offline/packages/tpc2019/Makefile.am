##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libtpc2019_io.la \
  libtpc2019.la

AM_CXXFLAGS = \
 -DRAVE -DRaveDllExport=

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  TpcPrototypeTrack.h \
  TpcPrototypeCluster.h \
  TpcPrototypeClusterizer.h \
  TpcPrototypeGenFitTrkFinder.h \
  TpcPrototypeGenFitTrkFitter.h \
  TpcPrototypeDefs.h \
  TpcPrototypeUnpacker.h

ROOTDICTS = \
  TpcPrototypeTrack_Dict.cc \
  TpcPrototypeCluster_Dict.cc \
  TpcPrototypeDefs_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA =  \
  TpcPrototypeTrack_Dict_rdict.pcm \
  TpcPrototypeCluster_Dict_rdict.pcm \
  TpcPrototypeDefs_Dict_rdict.pcm \
  TpcPrototypeUnpacker_Dict_rdict.pcm

libtpc2019_io_la_SOURCES = \
  $(ROOTDICTS) \
  TpcPrototypeTrack.cc \
  TpcPrototypeCluster.cc \
  TpcPrototypeDefs.cc

libtpc2019_la_SOURCES = \
  TpcPrototypeClusterizer.cc \
  TpcPrototypeGenFitTrkFinder.cc  \
  TpcPrototypeGenFitTrkFitter.cc  \
  TpcPrototypeUnpacker.cc  \
  TpcPrototypeUnpacker_Dict.cc

libtpc2019_io_la_LIBADD = \
  -ltpc_io \
  -ltrack_reco_io \
  -lphool

libtpc2019_la_LIBADD = \
  libtpc2019_io.la \
  -lgsl \
  -lgslcblas \
  -lgenfit2 \
  -lgenfit2exp \
  -lPHGenFit \
  -ltrackbase_historic_io \
  -lphparameter \
  -lphgeom \
  -ltpc

# Rule for generating CINT dictionaries from class headers.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


# link check
BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = \
  testexternals.cc

testexternals_io_LDADD = \
  libtpc2019_io.la

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libtpc2019.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
