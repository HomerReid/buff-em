EXTRA_DIST = 					\
 E10Sphere_533.buffgeo				\
 Sphere_533.msh    				\
 EPFile.XAxis               

LIBBUFF = $(top_builddir)/src/libs/libbuff/libbuff.la
AM_CPPFLAGS = -DBUFF -I$(top_srcdir)/src/libs/libbuff 

noinst_PROGRAMS = 		\
 unit-test-LFField       	

check_PROGRAMS = 		\
 unit-test-LFField

TESTS = 			\
 unit-test-LFField

unit_test_LFField_SOURCES = unit-test-LFField.cc
unit_test_LFField_LDADD   = $(LIBBUFF)
