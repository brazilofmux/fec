all: rsframe

.PHONY: all clean

clean:
	rm -f *.o rsframe

# librs headers come from ../librs (single source of truth) — this directory
# no longer carries its own copies.
rsframe: rsframe.o ../librs/librs.a
	g++ -O3 -fomit-frame-pointer -o rsframe rsframe.o ../librs/librs.a

rsframe.o: rsframe.cpp rsframe.h rsframe_processor.h ../librs/rs_factory.h ../librs/crc.h ../librs/projdefs.h
	g++ -O3 -fomit-frame-pointer -I../librs -c rsframe.cpp
