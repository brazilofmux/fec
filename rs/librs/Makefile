all: librs.a

librs.a: crc.o rs_common.o rs_factory.o rs_tables.o rs_decoder_general.o rs_encoder_general.o rs_encoder_general_descending.o rs_decoder_base.o rs_decoder_t1.o rs_encoder_t1.o rs_encoder_t2.o rs_encoder_t4.o rs_encoder_t8.o rs_encoder_t14.o rs_encoder_t16.o rs_encoder_t32.o rs_encoder_t64.o
	ar rvs librs.a crc.o rs_common.o rs_factory.o rs_tables.o rs_decoder_general.o rs_encoder_general.o rs_encoder_general_descending.o rs_decoder_base.o rs_decoder_t1.o rs_encoder_t1.o rs_encoder_t2.o rs_encoder_t4.o rs_encoder_t8.o rs_encoder_t14.o rs_encoder_t16.o rs_encoder_t32.o rs_encoder_t64.o

rs_common.o: rs_common.cpp rs_common.h
	g++ -O3 -fomit-frame-pointer -c rs_common.cpp

rs_factory.o: rs_factory.cpp rs_factory.h rs_tables.h rs_codec.h rs_encoder_t1.h rs_encoder_t2.h rs_encoder_t4.h rs_encoder_t8.h rs_encoder_t14.h rs_encoder_t16.h rs_encoder_t32.h rs_encoder_t64.h rs_encoder_general.h rs_decoder_t1.h rs_decoder_general.h
	g++ -O3 -fomit-frame-pointer -c rs_factory.cpp

rs_tables.o: rs_tables.cpp rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_tables.cpp

rs_decoder_general.o: rs_decoder_general.cpp rs_decoder_base.h rs_decoder_general.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_general.cpp

rs_encoder_general.o: rs_encoder_general.cpp rs_encoder_base.h rs_encoder_general.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_general.cpp

rs_encoder_general_descending.o: rs_encoder_general_descending.cpp rs_encoder_base.h rs_encoder_general_descending.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_general_descending.cpp

rs_decoder_base.o: rs_decoder_base.cpp rs_decoder_base.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_base.cpp

rs_decoder_t1.o: rs_decoder_t1.cpp rs_decoder_t1.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_t1.cpp

rs_encoder_t1.o: rs_encoder_t1.cpp rs_encoder_t1.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t1.cpp

rs_encoder_t2.o: rs_encoder_t2.cpp rs_encoder_t2.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t2.cpp

rs_encoder_t4.o: rs_encoder_t4.cpp rs_encoder_t4.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t4.cpp

rs_encoder_t8.o: rs_encoder_t8.cpp rs_encoder_t8.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t8.cpp

rs_encoder_t14.o: rs_encoder_t14.cpp rs_encoder_t14.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t14.cpp

rs_encoder_t16.o: rs_encoder_t16.cpp rs_encoder_t16.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t16.cpp

rs_encoder_t32.o: rs_encoder_t32.cpp rs_encoder_t32.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t32.cpp

rs_encoder_t64.o: rs_encoder_t64.cpp rs_encoder_t64.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_encoder_t64.cpp

crc.o: crc.cpp crc.h projdefs.h
	g++ -O3 -fomit-frame-pointer -c crc.cpp
