CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Iinclude

SRC = src/quadtree.cpp src/io.cpp experiments/quadtree_exp.cpp
OBJ = $(SRC:.cpp=.o)

all: run_experiments

run_experiments: $(OBJ)
	$(CXX) $(CXXFLAGS) -o run_experiments $(OBJ)

clean:
	rm -f $(OBJ) run_experiments
