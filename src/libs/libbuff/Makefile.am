lib_LTLIBRARIES = libbuff.la
pkginclude_HEADERS = libbuff.h IHAIMatProp.h 
libbuff_la_SOURCES = \
 GetFields.o    \
 Cubature.o     \
 IHAIMatProp.o  \
 InitFaceList.o \
 ReadGMSHFile.o \
 RHSVector.o    \
 SWGGeometry.o  \
 SWGVolume.o    \
 TetCR.o     	\
 VIEMatrix.o	\
 Visualize.o

# set the shared-library version number (DIFFERENT from human version #)
libbuff_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@
