
INITIAL  = testbinary
HYDRO    = euler
BOUNDARY = polar
OUTPUT   = ascii
RESTART  = no_restart
PLANETS  = sun

UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /home/install/app/hdf5
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/local
endif

CC = /usr/local/openmpi-1.4.3/bin/mpicc
CPP = /usr/local/openmpi-1.4.3/bin/mpic++
FOR = gfortran
FLAGS = -O3 -Wall -g

INC = 
LIB =  -lm -lgfortran
#-lhdf5 

OBJ = main.o readpar.o timestep.o onestep.o riemann.o mpisetup.o gridsetup.o domain.o misc.o geometry.o faces.o exchange.o plm.o report.o profiler.o planet.o omega.o analysis.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o $(BOUNDARY).o $(RESTART).o $(PLANETS).o helmholtz.o inverthelmeos.o BD.o interface.o iso7.o screen5.o sg_alloc.o sg_bndry.o sg_fft.o sg_fft.o sg_force.o sg_interp.o sg_routine.o

default: disco

%.o: %.c paul.h helm.h
	$(CC) $(FLAGS) $(INC) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(PLANETS).o : Planets/$(PLANETS).c paul.h
	$(CC) $(FLAGS) $(INC) -c Planets/$(PLANETS).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c paul.h
	$(CC) $(FLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

$(RESTART).o : Restart/$(RESTART).c paul.h
	$(CC) $(FLAGS) $(INC) -c Restart/$(RESTART).c

helmholtz.o : helmholtz/helmholtz.f90 
	$(FOR) -c helmholtz/helmholtz.f90

inverthelmeos.o : helmholtz/inverthelmeos.c helm.h	
	$(CC) -c helmholtz/inverthelmeos.c
	
BD.o : burning/BD.cpp nuclear.h
	$(CPP) $(FLAGS) $(INC) -c burning/BD.cpp
screen5.o : burning/screen5.cpp 
	$(CPP) $(FLAGS) $(INC) -c burning/screen5.cpp
interface.o : burning/interface.cpp nuclear.h paul.h
	$(CPP) $(FLAGS) $(INC) -c burning/interface.cpp
iso7.o : burning/iso7.cpp nuclear.h
	$(CPP) $(FLAGS) $(INC) -c burning/iso7.cpp	

%.o: SG/%.c SG.h paul.h
	$(CC) $(FLAGS) $(INC) -c $<
disco: $(OBJ) paul.h
	$(CPP) $(FLAGS)  -o disco $(OBJ) $(LIB)

clean:
	rm -f *.o disco
