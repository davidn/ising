## Process this file with automake to produce Makefile.in

## Created by Anjuta

AM_CPPFLAGS = \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir)"\" \
	-DPACKAGE_DATA_DIR=\""$(datadir)"\" 

AM_CFLAGS =\
	 -Wall\
	 -g

bin_PROGRAMS = ising

ising_SOURCES = \
	main.cc\
	../gnuplot-iostream/gnuplot-iostream.cc \
	../gnuplot-iostream/gnuplot-iostream.h \
	lattice.h \
	lattice.cc

ising_LDFLAGS = 

ising_LDADD = \
	-lboost_iostreams\
	-lutil

