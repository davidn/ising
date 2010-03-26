## Process this file with automake to produce Makefile.in

## Created by Anjuta

AM_CPPFLAGS = \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir)"\" \
	-DPACKAGE_DATA_DIR=\""$(datadir)"\" 

AM_CFLAGS =\
	 -Wall\
	 -g

bin_PROGRAMS = ising \
	runs

runs_SOURCES = \
	runs.cc\
	../gnuplot-iostream/gnuplot-iostream.cc \
	../gnuplot-iostream/gnuplot-iostream.h

runs_LDADD = \
	-lboost_iostreams \
	-lutil

ising_SOURCES = \
	ising.cc\
	lattice.h \
	lattice.cc \
	../gnuplot-iostream/gnuplot-iostream.cc \
	../gnuplot-iostream/gnuplot-iostream.h

ising_LDADD = \
	-lboost_iostreams \
	-lutil

ising_LDFLAGS = 

