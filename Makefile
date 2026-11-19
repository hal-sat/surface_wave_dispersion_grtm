cc = g++
CFLAGS = -Wall -Wextra -O3 -std=c++20 -mtune=native -march=native -pedantic

INCDIR = -I/usr/include/eigen3
SRCDIR = ./src/
BINDIR = ./bin/

TASKS = main
all: $(TASKS)

main:
	$(cc) $(CFLAGS) -o $(BINDIR)$@ $(INCDIR) $(SRCDIR)surface_wave_dispersion.cpp $(SRCDIR)$@.cpp

clean:
	rm -f *.o $(BINDIR)*
