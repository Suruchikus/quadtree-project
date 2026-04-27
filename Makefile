CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Iinclude

SRC = src/quadtree_v2.cpp src/io.cpp experiments/experiments_v2.cpp src/rank_support.cpp
OBJ = $(SRC:.cpp=.o)

all: run_experiments

run_experiments: $(OBJ)
	$(CXX) $(CXXFLAGS) -o run_experiments $(OBJ)

clean:
	rm -f $(OBJ) run_experiments
