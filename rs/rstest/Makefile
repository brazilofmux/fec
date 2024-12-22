CXX = g++
CXXFLAGS = -O3 -finline-functions -Wall

SRCS = rstest.cpp RsDebug.cpp RsVerify.cpp rs.cpp rsref.cpp rs_common.cpp RsVerification.cpp
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
