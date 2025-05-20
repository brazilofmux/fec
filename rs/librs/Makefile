all: librs.a

librs.a: rs.o crc.o
	ar rvs librs.a rs.o crc.o

rs.o: rs.cpp rs.h projdefs.h
	g++ -O3 -fomit-frame-pointer -c rs.cpp

crc.o: crc.cpp crc.h
	g++ -O3 -fomit-frame-pointer -c crc.cpp
