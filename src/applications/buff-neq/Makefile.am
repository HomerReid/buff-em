bin_PROGRAMS = buff-neq GetTemperatureDependence

buff_neq_SOURCES =	\
 buff-neq.cc		\
 Quadrature.cc		\
 CreateBNEQData.cc	\
 GetFlux.cc		\
 buff-neq.h

buff_neq_LDADD = $(top_builddir)/src/libs/libbuff/libbuff.la

GetTemperatureDependence_SOURCES = GetTemperatureDependence.cc

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libbuff
