CXX = g++
CC = gcc
AS = as
CXXFLAGS = -Wall -O2 -std=c++11 -I.
CFLAGS = -Wall -O2 -I.
LDFLAGS = -lm

# Source files
CPP_SRCS = viterbi27.cpp viterbi27_port.cpp viterbi27_sse2.cpp cpu_mode_x86.cpp fec.cpp sim.cpp
# ASM_SRCS = cpu_features.s sse2bfly27.s
ASM_SRCS = cpu_features.s
TEST_SRC = vtest27.cpp

# Object files
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
ASM_OBJS = $(ASM_SRCS:.s=.o)
TEST_OBJ = $(TEST_SRC:.cpp=.o)

# Library and test executable
LIB_NAME = libfec.a
TEST_NAME = vtest27

all: $(LIB_NAME) $(TEST_NAME)

$(LIB_NAME): $(CPP_OBJS) $(ASM_OBJS)
	ar rcs $@ $^

$(TEST_NAME): $(TEST_OBJ) $(LIB_NAME)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Rule for building .cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule for building .s files
%.o: %.s
	$(AS) $< -o $@

clean:
	rm -f $(CPP_OBJS) $(ASM_OBJS) $(TEST_OBJ) $(LIB_NAME) $(TEST_NAME)

test: $(TEST_NAME)
	./$(TEST_NAME)

.PHONY: all clean test
