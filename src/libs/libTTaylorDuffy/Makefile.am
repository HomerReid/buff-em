noinst_LTLIBRARIES = libTTaylorDuffy.la
pkginclude_HEADERS = TTaylorDuffy.h
libTTaylorDuffy_la_SOURCES = 	\
 TTaylorDuffy.cc      	\
 TTaylorDuffy.h

# the UpsilonFiles/Upsilon_*.cc files are #included by 
# TTaylorDuffy.cc, so I don't want to include them in 
# SOURCES above, but I do need them to be packaged together 
# with the distribution
SUBDIRS = UpsilonFiles
