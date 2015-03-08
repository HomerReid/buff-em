lib_LTLIBRARIES = libbuff.la
pkginclude_HEADERS = libbuff.h IHAIMatProp.h TTaylorDuffy.h
libbuff_la_SOURCES = 	\
 GetFields.cc    	\
 GMatrixElements.cc	\
 Cubature.cc     	\
 IHAIMatProp.cc  	\
 InitFaceList.cc 	\
 PFT.cc          	\
 ReadGMSHFile.cc 	\
 RHSVector.cc    	\
 SWGGeometry.cc  	\
 SWGVolume.cc    	\
 TaskTiming.cc		\
 TetCR.cc     		\
 TTaylorDuffy.cc      	\
 VIEMatrix.cc		\
 Visualize.cc		\
 TTaylorDuffy.h		\
 IHAIMatProp.h		\
 libbuff.h

# set the shared-library version number (DIFFERENT from human version #)
libbuff_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

# the UpsilonFiles/Upsilon_*.cc files are #included by 
# TTaylorDuffy.cc, so I don't want to include them in 
# SOURCES above, but I do need them to be packaged together 
# with the distribution
SUBDIRS = UpsilonFiles
