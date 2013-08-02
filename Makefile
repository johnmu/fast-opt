CC = g++
CFLAGS = -O2 -m64 -msse2 -lm -Wall -pthread -funroll-loops
#CFLAGS = -g3 -m64 -msse2 -lm -Wall -pthread

install: fast-opt

clean:
	rm -f fast-opt *.o

fast-opt: Makefile main.cpp main.h opt_tree.h copt_tree.h dfopt.h map_tree.h stl.h gamma_table.h general_utils.h disopt_tree.h llopt_tree.h lsopt_tree.h opt_utils.h density_store.h
	$(CC) $(CFLAGS) -o fast-opt main.cpp






