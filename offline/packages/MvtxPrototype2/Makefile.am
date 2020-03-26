AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmvtxprototype2.la

AM_CPPFLAGS = 								\
  -I$(includedir) 						\
  -I$(OFFLINE_MAIN)/include		\
  -I$(ROOTSYS)/include				\
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS =							\
  -L$(libdir) 						\
  -L$(ROOTSYS)/lib 				\
  -L$(OFFLINE_MAIN)/lib

libmvtxprototype2_la_LIBADD = \
  -lCLHEP 										\
  -lphool 										\
  -lfun4all 									\
  -lphparameter 							\
	-ltrack_io 									\
	-ltrackbase_historic_io 		\
  -lmvtx

if ! MAKEROOT6
ROOT5_DICTS = \
  MvtxPrototype2Align_Dict.C 				\
	MvtxPrototype2Geom_Dict.C 				\
	MvtxPrototype2UnpackPRDF_Dict.C 	\
	MvtxPrototype2Clusterizer_Dict.C 	\
  MvtxStandaloneTracking_Dict.C
endif

pkginclude_HEADERS = 					\
  MvtxPrototype2Align.h 			\
  MvtxPrototype2Geom.h 				\
  MvtxPrototype2UnpackPRDF.h 	\
  MvtxPrototype2Clusterizer.h \
  MvtxStandaloneTracking.h

libmvtxprototype2_la_SOURCES = \
	$(ROOT5_DICTS) 							 \
  MvtxPrototype2Align.cc 	     \
  MvtxPrototype2Geom.cc 	     \
  MvtxPrototype2UnpackPRDF.cc  \
  MvtxPrototype2Clusterizer.cc \
  MvtxStandaloneTracking.cc

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
	rm -f *Dict* *.pcm
