#include $(top_srcdir)/config/cheetah.mk
include $(top_srcdir)/config/flymake.am

AM_CPPFLAGS =  -I$(top_srcdir)/src -I$(CUDA_ROOT)/include @AM_CPPFLAGS@

# BUILT_SOURCES = specialization.hpp
CLEANFILES = $(BUILT_SOURCES)

noinst_LTLIBRARIES = libkernel.la

libkernel_la_SOURCES = link.cpp
# libkernel_la_SOURCES += recurrence.cu transfer.cu
# libkernel_la_SOURCES += ssss.cu fock.cu
# libkernel_la_SOURCES += psss.cu
# libkernel_la_SOURCES += dsss.cu
# libkernel_la_SOURCES += fsss.cu
libkernel_la_SOURCES += eri.cu
libkernel_la_SOURCES += fock.cu
libkernel_la_SOURCES += reduce.cu

headers = device.hpp quartet.hpp

# eri.o: eri.cu quartet.hpp
# eri.lo: eri.cu quartet.hpp

# .cu.o:
# 	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(AM_CPPFLAGS)
# #	source='$<' object='$@' libtool=no \
# #	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
# #	$(CXXCOMPILE) -c -o $@ $<

# .cu.lo: $(headers)
# 	$(top_srcdir)/config/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $(AM_CPPFLAGS) $<