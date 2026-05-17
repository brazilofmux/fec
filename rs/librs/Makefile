all: librs.a

# Architecture detection for optional SIMD compilation (AVX2 + AVX-512)
ifeq ($(shell uname -m),x86_64)
AVX2_FLAG    = -mavx2
AVX512_FLAGS = -mavx512f -mavx512bw
else
AVX2_FLAG    =
AVX512_FLAGS =
endif

librs.a: crc.o rs_factory.o rs_tables.o rs_decoder_general.o rs_decoder_direct.o rs_decoder_direct_neon.o rs_decoder_direct_avx2.o rs_decoder_direct_avx512.o rs_flipped_encoder_general.o rs_standard_encoder_general.o rs_decoder_base.o rs_decoder_t1.o rs_flipped_encoder_t1.o rs_flipped_encoder_t2.o rs_flipped_encoder_t4.o rs_flipped_encoder_t8.o rs_flipped_encoder_t14.o rs_flipped_encoder_t16.o rs_flipped_encoder_t32.o rs_flipped_encoder_t64.o
	ar rvs librs.a crc.o rs_factory.o rs_tables.o rs_decoder_general.o rs_decoder_direct.o rs_decoder_direct_neon.o rs_decoder_direct_avx2.o rs_decoder_direct_avx512.o rs_flipped_encoder_general.o rs_standard_encoder_general.o rs_decoder_base.o rs_decoder_t1.o rs_flipped_encoder_t1.o rs_flipped_encoder_t2.o rs_flipped_encoder_t4.o rs_flipped_encoder_t8.o rs_flipped_encoder_t14.o rs_flipped_encoder_t16.o rs_flipped_encoder_t32.o rs_flipped_encoder_t64.o

rs_factory.o: rs_factory.cpp rs_factory.h rs_tables.h rs_codec.h rs_flipped_encoder_t1.h rs_flipped_encoder_t2.h rs_flipped_encoder_t4.h rs_flipped_encoder_t8.h rs_flipped_encoder_t14.h rs_flipped_encoder_t16.h rs_flipped_encoder_t32.h rs_flipped_encoder_t64.h rs_flipped_encoder_general.h rs_standard_encoder_general.h rs_decoder_t1.h rs_decoder_general.h rs_decoder_direct.h rs_decoder_direct_neon.h rs_decoder_direct_avx2.h rs_decoder_direct_avx512.h
	g++ -O3 -fomit-frame-pointer -c rs_factory.cpp

rs_decoder_direct.o: rs_decoder_direct.cpp rs_decoder_base.h rs_decoder_direct.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_direct.cpp

rs_decoder_direct_neon.o: rs_decoder_direct_neon.cpp rs_decoder_base.h rs_decoder_direct_neon.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_direct_neon.cpp

rs_decoder_direct_avx2.o: rs_decoder_direct_avx2.cpp rs_decoder_base.h rs_decoder_direct_avx2.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer $(AVX2_FLAG) -c rs_decoder_direct_avx2.cpp

rs_decoder_direct_avx512.o: rs_decoder_direct_avx512.cpp rs_decoder_base.h rs_decoder_direct_avx512.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer $(AVX512_FLAGS) -c rs_decoder_direct_avx512.cpp

rs_tables.o: rs_tables.cpp rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_tables.cpp

rs_decoder_general.o: rs_decoder_general.cpp rs_decoder_base.h rs_decoder_general.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_general.cpp

rs_flipped_encoder_general.o: rs_flipped_encoder_general.cpp rs_encoder_base.h rs_flipped_encoder_general.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_general.cpp

rs_standard_encoder_general.o: rs_standard_encoder_general.cpp rs_encoder_base.h rs_standard_encoder_general.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_standard_encoder_general.cpp

rs_decoder_base.o: rs_decoder_base.cpp rs_decoder_base.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_base.cpp

rs_decoder_t1.o: rs_decoder_t1.cpp rs_decoder_t1.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_decoder_t1.cpp

rs_flipped_encoder_t1.o: rs_flipped_encoder_t1.cpp rs_flipped_encoder_t1.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t1.cpp

rs_flipped_encoder_t2.o: rs_flipped_encoder_t2.cpp rs_flipped_encoder_t2.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t2.cpp

rs_flipped_encoder_t4.o: rs_flipped_encoder_t4.cpp rs_flipped_encoder_t4.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t4.cpp

rs_flipped_encoder_t8.o: rs_flipped_encoder_t8.cpp rs_flipped_encoder_t8.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t8.cpp

rs_flipped_encoder_t14.o: rs_flipped_encoder_t14.cpp rs_flipped_encoder_t14.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t14.cpp

rs_flipped_encoder_t16.o: rs_flipped_encoder_t16.cpp rs_flipped_encoder_t16.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t16.cpp

rs_flipped_encoder_t32.o: rs_flipped_encoder_t32.cpp rs_flipped_encoder_t32.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t32.cpp

rs_flipped_encoder_t64.o: rs_flipped_encoder_t64.cpp rs_flipped_encoder_t64.h rs_common.h rs_tables.h
	g++ -O3 -fomit-frame-pointer -c rs_flipped_encoder_t64.cpp

crc.o: crc.cpp crc.h projdefs.h
	g++ -O3 -fomit-frame-pointer -c crc.cpp
