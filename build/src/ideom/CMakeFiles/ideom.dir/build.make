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
include src/ideom/CMakeFiles/ideom.dir/depend.make

# Include the progress variables for this target.
include src/ideom/CMakeFiles/ideom.dir/progress.make

# Include the compile flags for this target's objects.
include src/ideom/CMakeFiles/ideom.dir/flags.make

src/ideom/CMakeFiles/ideom.dir/ideomEquilibrium.cpp.o: src/ideom/CMakeFiles/ideom.dir/flags.make
src/ideom/CMakeFiles/ideom.dir/ideomEquilibrium.cpp.o: ../src/ideom/ideomEquilibrium.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/ideom/CMakeFiles/ideom.dir/ideomEquilibrium.cpp.o"
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ideom.dir/ideomEquilibrium.cpp.o -c /home/dhem/deom_git/HASSIM/src/ideom/ideomEquilibrium.cpp

src/ideom/CMakeFiles/ideom.dir/ideomEquilibrium.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ideom.dir/ideomEquilibrium.cpp.i"
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dhem/deom_git/HASSIM/src/ideom/ideomEquilibrium.cpp > CMakeFiles/ideom.dir/ideomEquilibrium.cpp.i

src/ideom/CMakeFiles/ideom.dir/ideomEquilibrium.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ideom.dir/ideomEquilibrium.cpp.s"
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dhem/deom_git/HASSIM/src/ideom/ideomEquilibrium.cpp -o CMakeFiles/ideom.dir/ideomEquilibrium.cpp.s

src/ideom/CMakeFiles/ideom.dir/ideomRem.cpp.o: src/ideom/CMakeFiles/ideom.dir/flags.make
src/ideom/CMakeFiles/ideom.dir/ideomRem.cpp.o: ../src/ideom/ideomRem.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/ideom/CMakeFiles/ideom.dir/ideomRem.cpp.o"
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ideom.dir/ideomRem.cpp.o -c /home/dhem/deom_git/HASSIM/src/ideom/ideomRem.cpp

src/ideom/CMakeFiles/ideom.dir/ideomRem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ideom.dir/ideomRem.cpp.i"
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dhem/deom_git/HASSIM/src/ideom/ideomRem.cpp > CMakeFiles/ideom.dir/ideomRem.cpp.i

src/ideom/CMakeFiles/ideom.dir/ideomRem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ideom.dir/ideomRem.cpp.s"
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dhem/deom_git/HASSIM/src/ideom/ideomRem.cpp -o CMakeFiles/ideom.dir/ideomRem.cpp.s

# Object files for target ideom
ideom_OBJECTS = \
"CMakeFiles/ideom.dir/ideomEquilibrium.cpp.o" \
"CMakeFiles/ideom.dir/ideomRem.cpp.o"

# External object files for target ideom
ideom_EXTERNAL_OBJECTS =

../lib/libideom.a: src/ideom/CMakeFiles/ideom.dir/ideomEquilibrium.cpp.o
../lib/libideom.a: src/ideom/CMakeFiles/ideom.dir/ideomRem.cpp.o
../lib/libideom.a: src/ideom/CMakeFiles/ideom.dir/build.make
../lib/libideom.a: src/ideom/CMakeFiles/ideom.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dhem/deom_git/HASSIM/b/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library ../../../lib/libideom.a"
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && $(CMAKE_COMMAND) -P CMakeFiles/ideom.dir/cmake_clean_target.cmake
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ideom.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/ideom/CMakeFiles/ideom.dir/build: ../lib/libideom.a

.PHONY : src/ideom/CMakeFiles/ideom.dir/build

src/ideom/CMakeFiles/ideom.dir/clean:
	cd /home/dhem/deom_git/HASSIM/b/src/ideom && $(CMAKE_COMMAND) -P CMakeFiles/ideom.dir/cmake_clean.cmake
.PHONY : src/ideom/CMakeFiles/ideom.dir/clean

src/ideom/CMakeFiles/ideom.dir/depend:
	cd /home/dhem/deom_git/HASSIM/b && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dhem/deom_git/HASSIM /home/dhem/deom_git/HASSIM/src/ideom /home/dhem/deom_git/HASSIM/b /home/dhem/deom_git/HASSIM/b/src/ideom /home/dhem/deom_git/HASSIM/b/src/ideom/CMakeFiles/ideom.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/ideom/CMakeFiles/ideom.dir/depend

