MPICXX ?= $(shell which mpicxx)
CXX = $(MPICXX)
CXXFLAGS = -Wall -O3 -std=c++20 -fopenmp

# Paths
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
RESULTS_DIR = results

# Libraries
PETSC_DIR ?= $(shell echo $$PETSC_DIR)
PETSC_ARCH ?= $(shell echo $$PETSC_ARCH)
LIBS = -L$(PETSC_DIR)/lib -lpetsc 

# Include directories
CXXFLAGS += -I$(INCLUDE_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/$(PETSC_ARCH) 

# Source and object files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Executable name
EXEC = $(BIN_DIR)/main

# Build executable
$(EXEC): $(OBJS) | $(BIN_DIR) $(RESULTS_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

# Build object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create directories if not exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(RESULTS_DIR):
	mkdir -p $(RESULTS_DIR)

# Clean rule
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXEC)

.PHONY: clean
