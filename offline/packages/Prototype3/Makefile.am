##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libPrototype3_io.la \
  libPrototype3.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

pkginclude_HEADERS = \
  CaloCalibration.h \
  CaloUnpackPRDF.h \
  EventInfoSummary.h \
  GenericUnpackPRDF.h \
  PROTOTYPE3_FEM.h \
  Prototype3DSTReader.h \
  RawTower_Prototype3.h \
  RawTower_Temperature.h \
  RunInfoUnpackPRDF.h \
  TempInfoUnpackPRDF.h

ROOTDICTS = \
  RawTower_Prototype3_Dict.cc \
  RawTower_Temperature_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  RawTower_Prototype3_Dict_rdict.pcm \
  RawTower_Temperature_Dict_rdict.pcm

libPrototype3_io_la_SOURCES = \
  $(ROOTDICTS) \
  PROTOTYPE3_FEM.cc \
  RawTower_Prototype3.cc \
  RawTower_Temperature.cc

libPrototype3_la_SOURCES = \
  Prototype3DSTReader.cc \
  CaloUnpackPRDF.cc \
  TempInfoUnpackPRDF.cc \
  RunInfoUnpackPRDF.cc \
  EventInfoSummary.cc \
  CaloCalibration.cc \
  GenericUnpackPRDF.cc

libPrototype3_io_la_LIBADD = \
  -lcalo_io \
  -lphool

libPrototype3_la_LIBADD = \
  libPrototype3_io.la \
  -lSubsysReco \
  -lfun4all \
  -lphparameter

# Rule for generating CINT dictionaries from class headers.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = \
  testexternals.cc

testexternals_io_LDADD = \
  libPrototype3_io.la

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libPrototype3.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
