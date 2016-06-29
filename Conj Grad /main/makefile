IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj

LIBS=-lm

_DEPS = Point.h Cube.h Nbhd.h Energy.h Descent.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

OBJ = main.o Point.o Cube.o Nbhd.o Energy.o Descent.o
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


#$(ODIR)/%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)

#$(ODIR): $(IDIR)/%.c $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: $(IDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

energy: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
	mv *.o ./$(ODIR)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ energy



























#CC=gcc
#CFLAGS=-lm
#DEPS = ../include/Point.h ../include/Cube.h ../include/Nbhd.h ../include/Energy.h ../include/Descent.h

#%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)

#energy: main.o Point.o Cube.o Nbhd.o Energy.o Descent.o
#	$(CC) -o energy main.o Point.o Cube.o Nbhd.o Energy.o Descent.o -lm

#clean: 
#	rm main.o Point.o Cube.o Nbhd.o Energy.o Descent.o energy
