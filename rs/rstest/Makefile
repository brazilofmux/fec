CXX = g++
CC = gcc
CXXFLAGS = -O3 -finline-functions -Wall
CFLAGS = -O3 -Wall

# Main test program
CPP_SRCS = rstest.cpp rs.cpp rsref.cpp rs_common.cpp RsComparisonTest.cpp \
           RsDecodingTest.cpp RsEncodingTest.cpp RsTestConfig.cpp RsVerification.cpp
C_SRCS = rskarn.c
OBJS = $(CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)
TARGET = rstest

# Encoder comparison program
ENCODER_CPP_SRCS = rs_common.cpp rs.cpp rsref.cpp rskarn.cpp compare_encoders.cpp RsVerification.cpp
ENCODER_OBJS = $(ENCODER_CPP_SRCS:.cpp=.o)
ENCODER_TEST = compare_encoders
FACTORY_LIB = ../librs/librs.a

.PHONY: all clean

all: $(TARGET) $(ENCODER_TEST)

$(TARGET): $(OBJS) ../librs/librs.a
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) ../librs/librs.a

$(ENCODER_TEST): $(ENCODER_OBJS) $(FACTORY_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(ENCODER_OBJS) $(FACTORY_LIB)

%.o: %.c
	$(CC) $(CFLAGS) -x c -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(ENCODER_OBJS) $(TARGET) $(ENCODER_TEST)
