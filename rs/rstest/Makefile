CXX = g++
CC = gcc
CXXFLAGS = -O3 -finline-functions -Wall
CFLAGS = -O3 -Wall

# Main test program
CPP_SRCS = rstest.cpp rs.cpp rsref.cpp rs_common.cpp RsComparisonTest.cpp \
           RsDecodingTest.cpp RsEncodingTest.cpp RsTestConfig.cpp RsVerification.cpp \
           rs_factory.cpp rs_tables.cpp rs_encoder_t1.cpp rs_decoder_general.cpp
C_SRCS = rskarn.c
OBJS = $(CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)
TARGET = rstest

# Encoder comparison program
ENCODER_CPP_SRCS = rs_common.cpp rs.cpp rsref.cpp rskarn.cpp compare_encoders.cpp RsVerification.cpp
ENCODER_C_SRCS =
ENCODER_OBJS = $(ENCODER_CPP_SRCS:.cpp=.o) $(ENCODER_C_SRCS:.c=.o)
ENCODER_TEST = compare_encoders

.PHONY: all clean

all: $(TARGET) $(ENCODER_TEST)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(ENCODER_TEST): $(ENCODER_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(ENCODER_OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -x c -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(ENCODER_OBJS) $(TARGET) $(ENCODER_TEST)
