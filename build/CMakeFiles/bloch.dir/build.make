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
include CMakeFiles/bloch.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bloch.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bloch.dir/flags.make

CMakeFiles/bloch.dir/src/apps/bloch.cpp.o: CMakeFiles/bloch.dir/flags.make
CMakeFiles/bloch.dir/src/apps/bloch.cpp.o: ../src/apps/bloch.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bloch.dir/src/apps/bloch.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bloch.dir/src/apps/bloch.cpp.o -c /home/dhem/deom_git/HASSIM/src/apps/bloch.cpp

CMakeFiles/bloch.dir/src/apps/bloch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bloch.dir/src/apps/bloch.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dhem/deom_git/HASSIM/src/apps/bloch.cpp > CMakeFiles/bloch.dir/src/apps/bloch.cpp.i

CMakeFiles/bloch.dir/src/apps/bloch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bloch.dir/src/apps/bloch.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dhem/deom_git/HASSIM/src/apps/bloch.cpp -o CMakeFiles/bloch.dir/src/apps/bloch.cpp.s

# Object files for target bloch
bloch_OBJECTS = \
"CMakeFiles/bloch.dir/src/apps/bloch.cpp.o"

# External object files for target bloch
bloch_EXTERNAL_OBJECTS =

../bin/bloch: CMakeFiles/bloch.dir/src/apps/bloch.cpp.o
../bin/bloch: CMakeFiles/bloch.dir/build.make
../bin/bloch: ../lib/libdeom.a
../bin/bloch: ../lib/libblockdeom.a
../bin/bloch: ../lib/libdeom2.a
../bin/bloch: ../lib/libideom.a
../bin/bloch: ../lib/libblockdeom2.a
../bin/bloch: CMakeFiles/bloch.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/bloch"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bloch.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bloch.dir/build: ../bin/bloch

.PHONY : CMakeFiles/bloch.dir/build

CMakeFiles/bloch.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bloch.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bloch.dir/clean

CMakeFiles/bloch.dir/depend:
	cd /home/dhem/deom_git/HASSIM/b && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b/CMakeFiles/bloch.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bloch.dir/depend

