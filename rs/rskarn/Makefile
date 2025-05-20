CC=gcc
CFLAGS=-O4

rstest: rstest.o rs.o
	$(CC) -o rstest rstest.o rs.o

rstest.o: rstest.c rs.h
rs.o: rs.c rs.h

