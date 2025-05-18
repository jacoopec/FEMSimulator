# Makefile for FEM 1D Project

# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

# Files
SRCS = main.cpp fem_utils.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = fem_program

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
