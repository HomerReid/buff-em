##################################################
##################################################
##################################################
SCUFFEM  = /home/homer/work/scuff-em-installation
#SCUFFEM  = /home/homer/work/scuff-em-debug
LIBDIR   = $(SCUFFEM)/lib

BUFFLIBDIR = ../../libs/libbuff
BUFFINCDIR = ../../libs/libbuff

##################################################
##################################################
##################################################
CPPFLAGS += -I$(SCUFFEM)/include/scuff-em -I$(HOME)/include -I. -fopenmp
CPPFLAGS += -I/home/homer/work/scuff-em/src/libs/libscuff
CPPFLAGS += -I$(BUFFINCDIR)
LDFLAGS += -L$(LIBDIR) -Wl,-rpath,$(LIBDIR)
LDFLAGS += -L$(HOME)/codes/lib -L$(HOME)/lib -L$(BUFFLIBDIR) -L.
CXXFLAGS += -O3 
#CXXFLAGS = -ggdb -O0

HR_LIBS=-lbuff -lscuff

RDL_LIBS=-lreadline -lncurses
HDF5_LIBS= -lhdf5_hl -lhdf5 
LB_LIBS=-llapack -lopenblas -lgomp -lgfortran 
OTHER_LIBS = $(RDL_LIBS) $(HDF5_LIBS) $(LB_LIBS)

LIBS = $(HR_LIBS) $(OTHER_LIBS)

##################################################
##################################################
##################################################
OBJS =              \
 buff-cas3D.o       \
 CasimirIntegrand.o \
 CreateBC3Data.o    \
 SumsIntegrals.o

buff-cas3D:  	$(OBJS)
		$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		/bin/rm *.o *.a
