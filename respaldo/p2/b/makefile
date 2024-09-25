# Definir el compilador y las banderas de compilación
CXX = g++
CXXFLAGS = -Wall -std=c++17 -Og -g

# Nombre del archivo ejecutable
target = miprograma.exe

# Buscar automáticamente todos los archivos fuente con extensión .cxx
SRCS = $(wildcard *.cxx spins.cc)

# Convertir los archivos fuente (.cxx) en archivos objeto (.o)
OBJS = $(patsubst %.cxx, %.o, $(SRCS))

# Regla para compilar y generar el ejecutable
$(target): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Regla para compilar cada archivo fuente .cxx en un archivo objeto .o
%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Limpiar los archivos objeto y otros archivos temporales
clean:
	@rm -f *.o *~

# Borra todos los archivos objeto, ejecutables y temporales
borra:
	@rm -f *.o *.exe *~

# Ejecutar el programa
run: $(target)
	@./$(target)

# Limpiar los archivos objeto y temporales
clear:
	@rm -f *.o *~
