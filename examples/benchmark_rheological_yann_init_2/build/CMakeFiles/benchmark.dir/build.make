# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build

# Include any dependencies generated for this target.
include CMakeFiles/benchmark.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/benchmark.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/benchmark.dir/flags.make

CMakeFiles/benchmark.dir/benchmark.cpp.o: CMakeFiles/benchmark.dir/flags.make
CMakeFiles/benchmark.dir/benchmark.cpp.o: ../benchmark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/benchmark.dir/benchmark.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmark.dir/benchmark.cpp.o -c /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/benchmark.cpp

CMakeFiles/benchmark.dir/benchmark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmark.dir/benchmark.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/benchmark.cpp > CMakeFiles/benchmark.dir/benchmark.cpp.i

CMakeFiles/benchmark.dir/benchmark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmark.dir/benchmark.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/benchmark.cpp -o CMakeFiles/benchmark.dir/benchmark.cpp.s

CMakeFiles/benchmark.dir/benchmark.cpp.o.requires:

.PHONY : CMakeFiles/benchmark.dir/benchmark.cpp.o.requires

CMakeFiles/benchmark.dir/benchmark.cpp.o.provides: CMakeFiles/benchmark.dir/benchmark.cpp.o.requires
	$(MAKE) -f CMakeFiles/benchmark.dir/build.make CMakeFiles/benchmark.dir/benchmark.cpp.o.provides.build
.PHONY : CMakeFiles/benchmark.dir/benchmark.cpp.o.provides

CMakeFiles/benchmark.dir/benchmark.cpp.o.provides.build: CMakeFiles/benchmark.dir/benchmark.cpp.o


CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o: CMakeFiles/benchmark.dir/flags.make
CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o: /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o -c /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp > CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.i

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp -o CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.s

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.requires:

.PHONY : CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.requires

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.provides: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.requires
	$(MAKE) -f CMakeFiles/benchmark.dir/build.make CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.provides.build
.PHONY : CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.provides

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.provides.build: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o


CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o: CMakeFiles/benchmark.dir/flags.make
CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o: /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o -c /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp > CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.i

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp -o CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.s

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.requires:

.PHONY : CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.requires

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.provides: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.requires
	$(MAKE) -f CMakeFiles/benchmark.dir/build.make CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.provides.build
.PHONY : CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.provides

CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.provides.build: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o


# Object files for target benchmark
benchmark_OBJECTS = \
"CMakeFiles/benchmark.dir/benchmark.cpp.o" \
"CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o" \
"CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o"

# External object files for target benchmark
benchmark_EXTERNAL_OBJECTS =

../benchmark: CMakeFiles/benchmark.dir/benchmark.cpp.o
../benchmark: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o
../benchmark: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o
../benchmark: CMakeFiles/benchmark.dir/build.make
../benchmark: libpalabos.a
../benchmark: /usr/lib/openmpi/lib/libmpi_cxx.so
../benchmark: /usr/lib/openmpi/lib/libmpi.so
../benchmark: /home/pedro/singularity/singularity-ce-3.8.1/workspace/LIGGGHTS-PUBLIC/src/liblmp_auto.so
../benchmark: CMakeFiles/benchmark.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable ../benchmark"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/benchmark.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/benchmark.dir/build: ../benchmark

.PHONY : CMakeFiles/benchmark.dir/build

CMakeFiles/benchmark.dir/requires: CMakeFiles/benchmark.dir/benchmark.cpp.o.requires
CMakeFiles/benchmark.dir/requires: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/liggghtsCouplingWrapper.cpp.o.requires
CMakeFiles/benchmark.dir/requires: CMakeFiles/benchmark.dir/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/src/latticeDecomposition.cpp.o.requires

.PHONY : CMakeFiles/benchmark.dir/requires

CMakeFiles/benchmark.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/benchmark.dir/cmake_clean.cmake
.PHONY : CMakeFiles/benchmark.dir/clean

CMakeFiles/benchmark.dir/depend:
	cd /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2 /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2 /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init_2/build/CMakeFiles/benchmark.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/benchmark.dir/depend

