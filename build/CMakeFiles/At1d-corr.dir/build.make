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
CMAKE_SOURCE_DIR = /home/dhem/deom_git/HASSIM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dhem/deom_git/HASSIM/b

# Include any dependencies generated for this target.
include CMakeFiles/At1d-corr.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/At1d-corr.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/At1d-corr.dir/flags.make

CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.o: CMakeFiles/At1d-corr.dir/flags.make
CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.o: ../src/apps/At1d-corr.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.o -c /home/dhem/deom_git/HASSIM/src/apps/At1d-corr.cpp

CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dhem/deom_git/HASSIM/src/apps/At1d-corr.cpp > CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.i

CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dhem/deom_git/HASSIM/src/apps/At1d-corr.cpp -o CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.s

# Object files for target At1d-corr
At1d__corr_OBJECTS = \
"CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.o"

# External object files for target At1d-corr
At1d__corr_EXTERNAL_OBJECTS =

../bin/At1d-corr: CMakeFiles/At1d-corr.dir/src/apps/At1d-corr.cpp.o
../bin/At1d-corr: CMakeFiles/At1d-corr.dir/build.make
../bin/At1d-corr: ../lib/libdeom.a
../bin/At1d-corr: ../lib/libblockdeom.a
../bin/At1d-corr: ../lib/libdeom2.a
../bin/At1d-corr: ../lib/libideom.a
../bin/At1d-corr: ../lib/libblockdeom2.a
../bin/At1d-corr: CMakeFiles/At1d-corr.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/At1d-corr"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/At1d-corr.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/At1d-corr.dir/build: ../bin/At1d-corr

.PHONY : CMakeFiles/At1d-corr.dir/build

CMakeFiles/At1d-corr.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/At1d-corr.dir/cmake_clean.cmake
.PHONY : CMakeFiles/At1d-corr.dir/clean

CMakeFiles/At1d-corr.dir/depend:
	cd /home/dhem/deom_git/HASSIM/b && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b/CMakeFiles/At1d-corr.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/At1d-corr.dir/depend

