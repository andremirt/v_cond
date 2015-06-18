AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@
include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = libbindings.la
libbindings_la_SOURCES = fortran.cpp

if LIBQC_WITH_GAMESS
lib_LTLIBRARIES = libqc_gamess.la
libqc_gamess_la_SOURCES = gamess.cpp
endif
