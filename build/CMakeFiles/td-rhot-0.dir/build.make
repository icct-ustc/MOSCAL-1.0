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
include CMakeFiles/td-rhot-0.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/td-rhot-0.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/td-rhot-0.dir/flags.make

CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.o: CMakeFiles/td-rhot-0.dir/flags.make
CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.o: ../src/apps/td-rhot-0.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.o -c /home/dhem/deom_git/HASSIM/src/apps/td-rhot-0.cpp

CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dhem/deom_git/HASSIM/src/apps/td-rhot-0.cpp > CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.i

CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dhem/deom_git/HASSIM/src/apps/td-rhot-0.cpp -o CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.s

# Object files for target td-rhot-0
td__rhot__0_OBJECTS = \
"CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.o"

# External object files for target td-rhot-0
td__rhot__0_EXTERNAL_OBJECTS =

../bin/td-rhot-0: CMakeFiles/td-rhot-0.dir/src/apps/td-rhot-0.cpp.o
../bin/td-rhot-0: CMakeFiles/td-rhot-0.dir/build.make
../bin/td-rhot-0: ../lib/libdeom.a
../bin/td-rhot-0: ../lib/libblockdeom.a
../bin/td-rhot-0: ../lib/libdeom2.a
../bin/td-rhot-0: ../lib/libideom.a
../bin/td-rhot-0: ../lib/libblockdeom2.a
../bin/td-rhot-0: CMakeFiles/td-rhot-0.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/td-rhot-0"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/td-rhot-0.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/td-rhot-0.dir/build: ../bin/td-rhot-0

.PHONY : CMakeFiles/td-rhot-0.dir/build

CMakeFiles/td-rhot-0.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/td-rhot-0.dir/cmake_clean.cmake
.PHONY : CMakeFiles/td-rhot-0.dir/clean

CMakeFiles/td-rhot-0.dir/depend:
	cd /home/dhem/deom_git/HASSIM/b && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b/CMakeFiles/td-rhot-0.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/td-rhot-0.dir/depend

