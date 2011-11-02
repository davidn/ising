## Process this file with automake to produce Makefile.in
## Created by Anjuta

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src \
	tests

isingdocdir = ${prefix}/doc/ising
isingdoc_DATA = \
	README\
	COPYING\
	AUTHORS\
	ChangeLog\
	INSTALL\
	NEWS


EXTRA_DIST = $(isingdoc_DATA)


# Remove doc directory on uninstall
uninstall-local:
	-rm -r $(isingdocdir)
