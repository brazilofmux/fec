CXX = g++
CXXFLAGS = -O3 -finline-functions -Wall

SRCS = rstest.cpp rs.cpp rsref.cpp rs_common.cpp RsComparisonTest.cpp RsTestConfig.cpp RsVerification.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = rstest

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
