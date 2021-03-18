# Compilador
CC=gcc

# Parametros para el compilador
CFLAGS=-Wall -Wextra

# Parametros para el preprocesador
# -Dnombre_a_definir=valor
# -Ipath_donde_encontrar_encabezados
CPPFLAGS=

# Bibliotecas a enlazar
# -Lpath_donde_encontrar_bibliotecas
# -lnombre_de_biblioteca
LDFLAGS=


TARGETS=matmul
SOURCES=matmul.c

# Make ya tiene una regla implicita que compila .c en .o usando
# los nombres de variables estandar de arriba

# Regla para enlazar matmul
matmul: $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


.PHONY=clean

clean:
	rm -f $(TARGETS) *.o

