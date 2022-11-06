AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4caloprototype.la

AM_CXXFLAGS = `geant4-config --cflags | sed s/-Woverloaded-virtual// | sed s/-W\ // | sed s/-Wshadow// `

AM_CPPFLAGS = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -isystem$(ROOTSYS)/include \
    -I$(G4_MAIN)/include

# set in configure.in to check if gcc version >= 4.8
#if GCC_GE_48
#  AM_CXXFLAGS = -std=c++11
#endif

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4caloprototype_la_LIBADD = \
  -lcalo_io \
  -lg4testbench \
  -lg4detectors \
  -lphool  \
  -lSubsysReco


pkginclude_HEADERS = \
  Prototype2RawTowerBuilder.h \
  PHG4PrototypeHcalDefs.h \
  PHG4Prototype2HcalCellReco.h \
  PHG4Prototype2OuterHcalSubsystem.h \
  PHG4Prototype2InnerHcalSubsystem.h \
  PHG4Prototype3InnerHcalSubsystem.h \
  PHG4SpacalPrototypeSubsystem.h \
  PHG4SpacalPrototype4Subsystem.h 

libg4caloprototype_la_SOURCES = \
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

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4caloprototype

testexternals_g4caloprototype_SOURCES = testexternals.cc
testexternals_g4caloprototype_LDADD = libg4caloprototype.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
