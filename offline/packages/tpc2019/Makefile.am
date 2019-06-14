##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libtpc2019_io.la \
  libtpc2019.la

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
  TpcPrototypeClusterizer.h

ROOTDICTS = 
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA =  
else
  ROOT5_IO_DICTS = 
  ROOT5_DICTS = \
    TpcPrototypeClusterizer_Dict.cc 
endif

libtpc2019_io_la_SOURCES = \
  $(ROOTDICTS) \
  $(ROOT5_IO_DICTS) 

libtpc2019_la_SOURCES = \
  $(ROOT5_DICTS) \
  TpcPrototypeClusterizer.cc 

libtpc2019_io_la_LIBADD = \
  -ltpc_io \
  -lphool

libtpc2019_la_LIBADD = \
  libtpc2019_io.la \
  -ltpc

# Rule for generating CINT dictionaries from class headers.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
