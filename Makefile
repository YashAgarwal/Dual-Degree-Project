# Makefile for diffusion microstructural evolution

# Compilation and linking options

COMPOPS = -g -Wall
LINKOPS = -lgsl -lgslcblas -lfftw3 -lfftw3_threads -lpthread -lm

# List of object files

objects = diffusion.o evolve.o functions.o

# List of header files

headers = stdio.h stdlib.h math.h complex.h gsl_rng.h fftw3.h \
functions.h

# List of source files

sources = diffusion.c evolve.c functions.c

# Directory paths for the source and header files

vpath %.c src
vpath %.h ./headers
vpath %.h /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/include/
vpath %.h /usr/include/
vpath %.h /usr/local/include/
vpath %.h /usr/local/include/gsl/

# Dependencies and the actions on the source files

diffusion: $(objects) $(headers)
	gcc -o spinodal.out $(objects) $(LINKOPS)
	rm -rf nohup.out

diffusion.o: $(sources) $(headers)
	gcc -o $@ -c ./src/diffusion.c $(COMPOPS)

evolve.o: evolve.c $(headers)
	gcc -o $@ -c ./src/evolve.c $(COMPOPS)

functions.o: functions.c $(headers)
	gcc -o $@ -c ./src/functions.c $(COMPOPS)

clean:
	rm -rf *.o

# End of Makefile
