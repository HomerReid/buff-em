EXTRA_DIST = 					\
 E10Sphere_533.buffgeo				\
 Sphere_533.msh    				\
 EPFile.XAxis

LIBBUFF = $(top_builddir)/src/libs/libbuff/libbuff.la
AM_CPPFLAGS = -DBUFF -I$(top_srcdir)/src/libs/libbuff 

noinst_PROGRAMS = 		\
 unit-test-LFField       	\
 unit-test-FIBBICache 

check_PROGRAMS = 		\
 unit-test-LFField		\
 unit-test-FIBBICache

TESTS = 			\
 unit-test-LFField		\
 unit-test-FIBBICache

unit_test_LFField_SOURCES = unit-test-LFField.cc
unit_test_LFField_LDADD   = $(LIBBUFF)

unit_test_FIBBICache_SOURCES = unit-test-FIBBICache.cc
unit_test_FIBBICache_LDADD   = $(LIBBUFF)
