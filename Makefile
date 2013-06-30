CC = g++

CFLAGS = -g -O2 -m64 -msse2 -lm -Wall -pthread
#CFLAGS = -g -m64 -msse2 -lm -Wall -pthread

install: opt-fast

clean:
	rm -f opt-fast *.o

opt-fast: Makefile main.cpp main.h opt_tree.h dfopt.h map_tree.h stl.h gamma_table.h general_utils.h disopt_tree.h llopt_tree.h lsopt_tree.h opt_utils.h
	$(CC) $(CFLAGS) -o opt-fast main.cpp






