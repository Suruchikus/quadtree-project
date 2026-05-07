CXX = g++
CXXFLAGS = -std=c++17 -O3 -DNDEBUG -march=native -Wall -Iinclude
LIBS = -lsdsl -ldivsufsort -ldivsufsort64

SRC = src/quadtree.cpp src/io.cpp experiments/experiments.cpp src/rank_support.cpp src/membership.cpp
OBJ = $(SRC:.cpp=.o)

all: run_experiments

run_experiments: $(OBJ)
	$(CXX) $(CXXFLAGS) -o run_experiments $(OBJ) $(LIBS)

clean:
	rm -f $(OBJ) run_experiments
