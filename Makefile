COMP = mpic++
OPTIMFLAGS = -O1 
OPENMPFLAG = -fopenmp
FFTFLAG = -lfftw3 -lfftw3_omp

OBJS = main.o functions.o global.o Events.o kmc.o


kmc.ex : $(OBJS)   
	$(COMP) -o kmc.ex main.o kmc.o functions.o global.o Events.o $(FFTFLAG) $(OPENMPFLAG) $(OPTIMFLAGS)

main.o : main.cpp kmc.h
	$(COMP) -c main.cpp  $(FFTFLAG) $(OPENMPFLAG) $(OPTIMFLAGS)

functions.o : functions.cpp functions.h
	$(COMP) -c functions.cpp $(FFTFLAG) $(OPENMPFLAG) $(OPTIMFLAGS)

global.o : global.h global.cpp functions.h
	$(COMP) -c global.cpp $(FFTFLAG) $(OPENMPFLAG) $(OPTIMFLAGS) 

Events.o : Events.h Events.cpp functions.h
	$(COMP) -c Events.cpp $(FFTFLAG) $(OPENMPFLAG) $(OPTIMFLAGS)

kmc.o : Events.h Adatom.h Island.h global.h kmc.cpp
	$(COMP) -c kmc.cpp $(FFTFLAG) $(OPENMPFLAG) $(OPTIMFLAGS)

clean:
	rm *.o
