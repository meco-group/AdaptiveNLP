# Additional libraries
ADDLIBS = -lcasadi -lpython3.10

# Internal directories
OBJ_DIR = src/obj
SRC_DIR = src
INCL_DIR = include

# Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = #-I /mnt/c/users/u0143705/Documents/PhD/matplotlib-cpp -I /usr/include/python3.10 -I /home/u0143705/.local/lib/python3.10/site-packages/numpy/core/include

OBJS_CORRIDOR = adaptiveCorridorExample.o buildingBlocks.o makeBuildingBlocks.o bookkeeper.o NLPInterface.o adaptiveNLP.o interfaceTester.o plotter.o adaptiveCorridorExampleHelpers.o
OBJS_MOONLANDER = moonlanderExample.o buildingBlocks.o makeBuildingBlocks.o bookkeeper.o NLPInterface.o adaptiveNLP.o interfaceTester.o plotter.o errorEstimator.o moonlanderExampleHelper.o
##########################################################################

OBJ_CORRIDOR = $(patsubst %,$(OBJ_DIR)/%, $(OBJS_CORRIDOR))
OBJ_MOONLANDER = $(patsubst %,$(OBJ_DIR)/%, $(OBJS_MOONLANDER))


# C++ Compiler command
CXX = g++

# C++ Compiler options
# CXXFLAGS = -O2 -DNDEBUG
# CXXFLAGS = -O2 -pg
CXXFLAGS = -O2

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/usr/local/lib

prefix = /usr/local
exec_prefix = ${prefix}

# Include directories
INCL = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --libs ipopt`

# Compilation rule to create object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCL_DIR)/%.hpp
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<

MPC-Example: $(OBJ_CORRIDOR)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $(OBJ_DIR)/$@ $(OBJ_CORRIDOR) $(ADDLIBS) $(LIBS)

MoonlanderExample: $(OBJ_MOONLANDER)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $(OBJ_DIR)/$@ $(OBJ_MOONLANDER) $(ADDLIBS) $(LIBS)

all: MPC-Example MoonlanderExample

.PHONY: clean

clean:
	rm -rf $(OBJ_DIR)/*.o
