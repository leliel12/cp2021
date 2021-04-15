CC=clang -g 
### CFLAGS= -O3 -march=native
CFLAGS= -O3 -march=native
LDFLAGS=

TARGETS=demo headless headless_peano
SOURCES=$(shell echo *.c)
COMMON_OBJECTS=solver.o
COMMON_OBJECTS2=solver_peano.o

all: $(TARGETS)

demo: demo.o $(COMMON_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lGL -lGLU -lglut

headless: headless.o $(COMMON_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

headless_peano: headless_peano.o $(COMMON_OBJECTS2)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
clean:
	rm -f $(TARGETS) *.o .depend *~

.depend: *.[ch]
	$(CC) -MM $(SOURCES) >.depend

compara: headless headless_peano
	 ./headless_peano 
	 ./headless
	 sh ./reorder_peano.sh

-include .depend

.PHONY: clean all
