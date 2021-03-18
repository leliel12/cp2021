CC=gcc

solucion0: clean
	$(CC) -o $@ solucion0.c
	./$@ || echo -n

solucion1: clean
	$(CC) -O1 -o $@ solucion0.c
	./$@ || echo -n

solucion2: clean
	$(CC) -O2 -o $@ solucion0.c
	./$@ || echo -n

solucion3: clean
	$(CC) -O3 -o $@ solucion0.c
	./$@ || echo -n

solucion4: clean
	$(CC) -O3 -march=native -o $@ solucion0.c
	./$@ || echo -n

solucion5: clean
	$(CC) -O3 -march=native -o $@ solucion5.c
	./$@ || echo -n

solucion6: clean
	$(CC) -O2 -o $@ solucion6.c -lopenblas
	OPENBLAS_NUM_THREADS=1 ./$@ || echo -n

solucion7: clean
	$(CC) -O2 -o $@ solucion6.c -lopenblas
	for i in $(shell seq 2 `grep processor /proc/cpuinfo | wc -l`); do echo $$i hilos; OPENBLAS_NUM_THREADS=$$i ./$@ || echo -n; done

.PHONY=clean

TARGETS=solucion0 solucion1 solucion2 solucion3 solucion4 solucion5 solucion6 solucion7

clean:
	rm -f $(TARGETS) *.o
