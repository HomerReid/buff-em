lib_LTLIBRARIES = libbuff.la
pkginclude_HEADERS = libbuff.h SVTensor.h FIBBICache.h
libbuff_la_SOURCES = 	\
 GetFields.cc    	\
 GMatrixElements.cc	\
 Cubature.cc     	\
 FIBBICache.cc   	\
 SVTensor.cc     	\
 InitFaceList.cc 	\
 ReadGMSHFile.cc 	\
 RHSVector.cc    	\
 SWGGeometry.cc  	\
 SWGVolume.cc    	\
 TetCR.cc     		\
 VIEMatrix.cc		\
 Visualize.cc		\
 OverlapIntegrals.cc	\
 OPFT.cc         	\
 JDEPFT.cc         	\
 DSIPFT.cc       	\
 MomentPFT.cc     	\
 GetPFT.cc          	\
 SVTensor.h    		\
 FIBBICache.h		\
 libbuff.h

# combine all auxiliary libraries into a single library 
libbuff_la_LIBADD = ../libTTaylorDuffy/libTTaylorDuffy.la

# set the shared-library version number (DIFFERENT from human version #)
libbuff_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

AM_CPPFLAGS = -I../libTTaylorDuffy
