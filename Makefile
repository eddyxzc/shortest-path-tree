CC=gcc
CFLAGS=-Wall
OFLAGS=-O

all: sssp

sssp: main.c
	$(CC) $(CFLAGS) $(OFLAGS) -o sssp main.c

clean:
	rm sssp