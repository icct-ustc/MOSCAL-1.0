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
include CMakeFiles/weak_measure_corr_1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/weak_measure_corr_1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/weak_measure_corr_1.dir/flags.make

CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.o: CMakeFiles/weak_measure_corr_1.dir/flags.make
CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.o: ../src/apps/weak_measure_corr_1.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.o -c /home/dhem/deom_git/HASSIM/src/apps/weak_measure_corr_1.cpp

CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dhem/deom_git/HASSIM/src/apps/weak_measure_corr_1.cpp > CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.i

CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dhem/deom_git/HASSIM/src/apps/weak_measure_corr_1.cpp -o CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.s

# Object files for target weak_measure_corr_1
weak_measure_corr_1_OBJECTS = \
"CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.o"

# External object files for target weak_measure_corr_1
weak_measure_corr_1_EXTERNAL_OBJECTS =

../bin/weak_measure_corr_1: CMakeFiles/weak_measure_corr_1.dir/src/apps/weak_measure_corr_1.cpp.o
../bin/weak_measure_corr_1: CMakeFiles/weak_measure_corr_1.dir/build.make
../bin/weak_measure_corr_1: ../lib/libdeom.a
../bin/weak_measure_corr_1: ../lib/libblockdeom.a
../bin/weak_measure_corr_1: ../lib/libdeom2.a
../bin/weak_measure_corr_1: ../lib/libideom.a
../bin/weak_measure_corr_1: ../lib/libblockdeom2.a
../bin/weak_measure_corr_1: CMakeFiles/weak_measure_corr_1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/weak_measure_corr_1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/weak_measure_corr_1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/weak_measure_corr_1.dir/build: ../bin/weak_measure_corr_1

.PHONY : CMakeFiles/weak_measure_corr_1.dir/build

CMakeFiles/weak_measure_corr_1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/weak_measure_corr_1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/weak_measure_corr_1.dir/clean

CMakeFiles/weak_measure_corr_1.dir/depend:
	cd /home/dhem/deom_git/HASSIM/b && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b/CMakeFiles/weak_measure_corr_1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/weak_measure_corr_1.dir/depend

