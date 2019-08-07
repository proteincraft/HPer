CFLAGS=-g -Wall -lm

all:globalhepo

globalhepo: protein.o geometry.o mem.o read_structure.o score.o surface.o paxis.o globalhepo4.o wyang.o
	mpicc protein.o geometry.o mem.o read_structure.o score.o surface.o paxis.o globalhepo4.o wyang.o -lm -o globalhepo

globalhepo4.o: globalhepo4.c
	mpicc -lm -c globalhepo4.c

protein.o: protein.c
	gcc -lm -c protein.c

geometry.o: geometry.c
	gcc -lm -c geometry.c

mem.o: mem.c
	gcc -lm -c mem.c

read_structure.o: read_structure.c
	gcc -lm -c read_structure.c

surface.o: surface.c
	gcc -lm -c surface.c
score.o:score.c
	gcc -lm -c score.c
paxis.o:paxis.c
	mpicc -lm -c paxis.c
wyang.o:wyang.c
	gcc -lm -c wyang.c
