# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/eddy/sssp/sssp-v2/testing/RandGraph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/eddy/sssp/sssp-v2/testing/RandGraph/build

# Include any dependencies generated for this target.
include CMakeFiles/randGraph.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/randGraph.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/randGraph.dir/flags.make

CMakeFiles/randGraph.dir/main.cpp.o: CMakeFiles/randGraph.dir/flags.make
CMakeFiles/randGraph.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/eddy/sssp/sssp-v2/testing/RandGraph/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/randGraph.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/randGraph.dir/main.cpp.o -c /home/eddy/sssp/sssp-v2/testing/RandGraph/main.cpp

CMakeFiles/randGraph.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/randGraph.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/eddy/sssp/sssp-v2/testing/RandGraph/main.cpp > CMakeFiles/randGraph.dir/main.cpp.i

CMakeFiles/randGraph.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/randGraph.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/eddy/sssp/sssp-v2/testing/RandGraph/main.cpp -o CMakeFiles/randGraph.dir/main.cpp.s

# Object files for target randGraph
randGraph_OBJECTS = \
"CMakeFiles/randGraph.dir/main.cpp.o"

# External object files for target randGraph
randGraph_EXTERNAL_OBJECTS =

randGraph: CMakeFiles/randGraph.dir/main.cpp.o
randGraph: CMakeFiles/randGraph.dir/build.make
randGraph: CMakeFiles/randGraph.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/eddy/sssp/sssp-v2/testing/RandGraph/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable randGraph"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/randGraph.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/randGraph.dir/build: randGraph

.PHONY : CMakeFiles/randGraph.dir/build

CMakeFiles/randGraph.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/randGraph.dir/cmake_clean.cmake
.PHONY : CMakeFiles/randGraph.dir/clean

CMakeFiles/randGraph.dir/depend:
	cd /home/eddy/sssp/sssp-v2/testing/RandGraph/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/eddy/sssp/sssp-v2/testing/RandGraph /home/eddy/sssp/sssp-v2/testing/RandGraph /home/eddy/sssp/sssp-v2/testing/RandGraph/build /home/eddy/sssp/sssp-v2/testing/RandGraph/build /home/eddy/sssp/sssp-v2/testing/RandGraph/build/CMakeFiles/randGraph.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/randGraph.dir/depend
