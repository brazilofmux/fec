all: rsframe

rsframe: rsframe.o ../librs/librs.a
	g++ -O3 -fomit-frame-pointer -o rsframe rsframe.o ../librs/librs.a

rsframe.o: rsframe.cpp crc.h rsframe.h projdefs.h
	g++ -O3 -fomit-frame-pointer -c rsframe.cpp
