AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4caloprototype_io.la \
    libg4caloprototype.la

AM_CPPFLAGS = \
    -DCGAL_DISABLE_ROUNDING_MATH_CHECK \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I$(G4_MAIN)/include \
    -I$(OPT_SPHENIX)/include

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4caloprototype_io_la_LIBADD = \
  -lphool \
  -lphparameter \
  -lg4detectors_io 


libg4caloprototype_la_LIBADD = \
  libg4caloprototype_io.la \
  -lphg4gdml \
  -lphool  \
  -lCGAL \
  -lSubsysReco \
  -lg4testbench \
  -lg4detectors

pkginclude_HEADERS = \
  Prototype2RawTowerBuilder.h \
  PHG4PrototypeHcalDefs.h \
  PHG4Prototype2HcalCellReco.h \
  PHG4Prototype2OuterHcalSubsystem.h \
  PHG4Prototype2InnerHcalSubsystem.h \
  PHG4Prototype3InnerHcalSubsystem.h \
  PHG4SpacalPrototypeSubsystem.h \
  PHG4SpacalPrototype4Subsystem.h 

# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOT_IO_DICTS = 
# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# MAKEROOT6 is set in the configure.ac, true for root6
if MAKEROOT6
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
pcmdir = $(libdir)
nobase_dist_pcm_DATA = 
else
  ROOT5_IO_DICTS = 
  ROOT5_DICTS = \
	  Prototype2RawTowerBuilder_Dict.cc \
	  PHG4Prototype2InnerHcalSubsystem_Dict.cc \
	  PHG4Prototype2HcalCellReco_Dict.cc \
	  PHG4Prototype2OuterHcalSubsystem_Dict.cc \
	  PHG4Prototype3InnerHcalSubsystem_Dict.cc \
	  PHG4SpacalPrototypeSubsystem_Dict.cc \
	  PHG4SpacalPrototype4Subsystem_Dict.cc

endif

libg4caloprototype_io_la_SOURCES = \
  $(ROOT_IO_DICTS) \
  $(ROOT5_IO_DICTS) 

libg4caloprototype_la_SOURCES = \
  $(ROOT5_DICTS) \
  Prototype2RawTowerBuilder.cc \
  PHG4Prototype2InnerHcalDetector.cc \
  PHG4Prototype2InnerHcalSteppingAction.cc \
  PHG4Prototype2InnerHcalSubsystem.cc \
  PHG4Prototype2HcalCellReco.cc \
  PHG4Prototype2OuterHcalDetector.cc \
  PHG4Prototype2OuterHcalSteppingAction.cc \
  PHG4Prototype2OuterHcalSubsystem.cc \
  PHG4Prototype3InnerHcalDetector.cc \
  PHG4Prototype3InnerHcalSteppingAction.cc \
  PHG4Prototype3InnerHcalSubsystem.cc \
  PHG4HcalPrototypeDetector.cc \
  PHG4HcalPrototypeDetectorMessenger.cc \
  PHG4HcalPrototypeSteppingAction.cc \
  PHG4HcalPrototypeSubsystem.cc \
  PHG4SpacalPrototypeDetector.cc \
  PHG4SpacalPrototypeSteppingAction.cc \
  PHG4SpacalPrototypeSubsystem.cc \
  PHG4SpacalPrototype4Detector.cc \
  PHG4SpacalPrototype4SteppingAction.cc \
  PHG4SpacalPrototype4Subsystem.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4caloprototype \
  testexternals_g4caloprototype_io

testexternals_g4caloprototype_SOURCES = testexternals.cc
testexternals_g4caloprototype_LDADD = libg4caloprototype.la

testexternals_g4caloprototype_io_SOURCES = testexternals.cc
testexternals_g4caloprototype_io_LDADD = libg4caloprototype_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
