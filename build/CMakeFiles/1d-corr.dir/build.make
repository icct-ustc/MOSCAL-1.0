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
include CMakeFiles/1d-corr.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/1d-corr.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/1d-corr.dir/flags.make

CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.o: CMakeFiles/1d-corr.dir/flags.make
CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.o: ../src/apps/1d-corr.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.o -c /home/dhem/deom_git/HASSIM/src/apps/1d-corr.cpp

CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dhem/deom_git/HASSIM/src/apps/1d-corr.cpp > CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.i

CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dhem/deom_git/HASSIM/src/apps/1d-corr.cpp -o CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.s

# Object files for target 1d-corr
1d__corr_OBJECTS = \
"CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.o"

# External object files for target 1d-corr
1d__corr_EXTERNAL_OBJECTS =

../bin/1d-corr: CMakeFiles/1d-corr.dir/src/apps/1d-corr.cpp.o
../bin/1d-corr: CMakeFiles/1d-corr.dir/build.make
../bin/1d-corr: ../lib/libdeom.a
../bin/1d-corr: ../lib/libblockdeom.a
../bin/1d-corr: ../lib/libdeom2.a
../bin/1d-corr: ../lib/libideom.a
../bin/1d-corr: ../lib/libblockdeom2.a
../bin/1d-corr: CMakeFiles/1d-corr.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/1d-corr"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/1d-corr.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/1d-corr.dir/build: ../bin/1d-corr

.PHONY : CMakeFiles/1d-corr.dir/build

CMakeFiles/1d-corr.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/1d-corr.dir/cmake_clean.cmake
.PHONY : CMakeFiles/1d-corr.dir/clean

CMakeFiles/1d-corr.dir/depend:
	cd /home/dhem/deom_git/HASSIM/b && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b/CMakeFiles/1d-corr.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/1d-corr.dir/depend

