include $(top_srcdir)/config/flymake.mk

lib_LTLIBRARIES = libqc.la
libqc_la_SOURCES = mp2.cpp
libqc_la_LIBADD = core/libcore.la util/libutil.la bindings/libbindings.la
libqc_la_LIBADD += hf/libhf.la
libqc_la_LIBADD += matrix/libmatrix.la

SUBDIRS = core hf util bindings matrix adapter

# libqc_la_LIBADD += contrib/libcontrib.la
# SUBDIRS += contrib
