##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libPrototype2_io.la \
  libPrototype2.la

AM_CXXFLAGS = \
 -Werror

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`


pkginclude_HEADERS = \
  CaloCalibration.h \
  CaloUnpackPRDF.h \
  GenericUnpackPRDF.h \
  PROTOTYPE2_FEM.h \
  Prototype2DSTReader.h \
  RawTower_Prototype2.h \
  RawTower_Temperature.h \
  RunInfoUnpackPRDF.h \
  TempInfoUnpackPRDF.h

ROOTDICTS = \
  RawTower_Prototype2_Dict.cc \
  RawTower_Temperature_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    RawTower_Prototype2_Dict_rdict.pcm \
    RawTower_Temperature_Dict_rdict.pcm
else
  ROOT5_IO_DICTS = \
    PROTOTYPE2_FEM_Dict.cc
  ROOT5_DICTS = \
    CaloCalibration_Dict.cc \
    CaloUnpackPRDF_Dict.cc \
    GenericUnpackPRDF_Dict.cc \
    Prototype2DSTReader_Dict.cc \
    RunInfoUnpackPRDF_Dict.cc \
    TempInfoUnpackPRDF_Dict.cc
endif

libPrototype2_io_la_SOURCES = \
  $(ROOTDICTS) \
  $(ROOT5_IO_DICTS) \
  PROTOTYPE2_FEM.cc \
  RawTower_Prototype2.cc \
  RawTower_Temperature.cc

libPrototype2_la_SOURCES = \
  $(ROOT5_DICTS) \
  CaloCalibration.cc \
  CaloUnpackPRDF.cc \
  GenericUnpackPRDF.cc \
  TempInfoUnpackPRDF.cc \
  Prototype2DSTReader.cc \
  RunInfoUnpackPRDF.cc

libPrototype2_io_la_LIBADD = \
  -lcalo_io \
  -lphool

libPrototype2_la_LIBADD = \
  libPrototype2_io.la \
  -lSubsysReco \
  -lphool \
  -lfun4all \
  -lcalo_util \
  -lphparameter

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = \
  testexternals.cc

testexternals_io_LDADD = \
  libPrototype2_io.la

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libPrototype2.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
