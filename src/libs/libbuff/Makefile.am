lib_LTLIBRARIES = libbuff.la
pkginclude_HEADERS = libbuff.h IHAIMatProp.h FIBBICache.h
libbuff_la_SOURCES = 	\
 GetFields.cc    	\
 GMatrixElements.cc	\
 Cubature.cc     	\
 FIBBICache.cc   	\
 IHAIMatProp.cc  	\
 InitFaceList.cc 	\
 ReadGMSHFile.cc 	\
 RHSVector.cc    	\
 SWGGeometry.cc  	\
 SWGVolume.cc    	\
 TaskTiming.cc		\
 TetCR.cc     		\
 VIEMatrix.cc		\
 Visualize.cc		\
 OPFT.cc         	\
 JDEPFT.cc         	\
 DSIPFT.cc       	\
 GetPFT.cc          	\
 IHAIMatProp.h		\
 FIBBICache.h		\
 rwlock.h    		\
 libbuff.h

# combine all auxiliary libraries into a single library 
libbuff_la_LIBADD = ../libTTaylorDuffy/libTTaylorDuffy.la

# set the shared-library version number (DIFFERENT from human version #)
libbuff_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

AM_CPPFLAGS = -I../libTTaylorDuffy
