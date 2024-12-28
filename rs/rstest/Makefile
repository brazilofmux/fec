CXX = g++
CXXFLAGS = -O3 -finline-functions -Wall

SRCS = rstest.cpp rs.cpp rsref.cpp rs_common.cpp RsComparisonTest.cpp RsDecodingTest.cpp RsEncodingTest.cpp RsTestConfig.cpp RsVerification.cpp rskarn.c
OBJS = $(SRCS:.cpp=.o)
TARGET = rstest

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
