# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /gpfs/home/pgersberg/libeigen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /gpfs/home/pgersberg/libeigen/buikd_dir

# Include any dependencies generated for this target.
include doc/examples/CMakeFiles/class_VectorBlock.dir/depend.make

# Include the progress variables for this target.
include doc/examples/CMakeFiles/class_VectorBlock.dir/progress.make

# Include the compile flags for this target's objects.
include doc/examples/CMakeFiles/class_VectorBlock.dir/flags.make

doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o: doc/examples/CMakeFiles/class_VectorBlock.dir/flags.make
doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o: ../doc/examples/class_VectorBlock.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /gpfs/home/pgersberg/libeigen/buikd_dir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o -c /gpfs/home/pgersberg/libeigen/doc/examples/class_VectorBlock.cpp

doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.i"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /gpfs/home/pgersberg/libeigen/doc/examples/class_VectorBlock.cpp > CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.i

doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.s"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /gpfs/home/pgersberg/libeigen/doc/examples/class_VectorBlock.cpp -o CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.s

doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.requires:
.PHONY : doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.requires

doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.provides: doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.requires
	$(MAKE) -f doc/examples/CMakeFiles/class_VectorBlock.dir/build.make doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.provides.build
.PHONY : doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.provides

doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.provides.build: doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o

# Object files for target class_VectorBlock
class_VectorBlock_OBJECTS = \
"CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o"

# External object files for target class_VectorBlock
class_VectorBlock_EXTERNAL_OBJECTS =

doc/examples/class_VectorBlock: doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o
doc/examples/class_VectorBlock: doc/examples/CMakeFiles/class_VectorBlock.dir/build.make
doc/examples/class_VectorBlock: doc/examples/CMakeFiles/class_VectorBlock.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable class_VectorBlock"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/class_VectorBlock.dir/link.txt --verbose=$(VERBOSE)
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && ./class_VectorBlock >/gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples/class_VectorBlock.out

# Rule to build all files generated by this target.
doc/examples/CMakeFiles/class_VectorBlock.dir/build: doc/examples/class_VectorBlock
.PHONY : doc/examples/CMakeFiles/class_VectorBlock.dir/build

doc/examples/CMakeFiles/class_VectorBlock.dir/requires: doc/examples/CMakeFiles/class_VectorBlock.dir/class_VectorBlock.cpp.o.requires
.PHONY : doc/examples/CMakeFiles/class_VectorBlock.dir/requires

doc/examples/CMakeFiles/class_VectorBlock.dir/clean:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && $(CMAKE_COMMAND) -P CMakeFiles/class_VectorBlock.dir/cmake_clean.cmake
.PHONY : doc/examples/CMakeFiles/class_VectorBlock.dir/clean

doc/examples/CMakeFiles/class_VectorBlock.dir/depend:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/home/pgersberg/libeigen /gpfs/home/pgersberg/libeigen/doc/examples /gpfs/home/pgersberg/libeigen/buikd_dir /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples/CMakeFiles/class_VectorBlock.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/examples/CMakeFiles/class_VectorBlock.dir/depend

