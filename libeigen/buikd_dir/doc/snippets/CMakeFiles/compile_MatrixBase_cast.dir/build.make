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
include doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/depend.make

# Include the progress variables for this target.
include doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/progress.make

# Include the compile flags for this target's objects.
include doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/flags.make

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o: doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/flags.make
doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o: doc/snippets/compile_MatrixBase_cast.cpp
doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o: ../doc/snippets/MatrixBase_cast.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /gpfs/home/pgersberg/libeigen/buikd_dir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o -c /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets/compile_MatrixBase_cast.cpp

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.i"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets/compile_MatrixBase_cast.cpp > CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.i

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.s"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets/compile_MatrixBase_cast.cpp -o CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.s

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.requires:
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.requires

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.provides: doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.requires
	$(MAKE) -f doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/build.make doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.provides.build
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.provides

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.provides.build: doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o

# Object files for target compile_MatrixBase_cast
compile_MatrixBase_cast_OBJECTS = \
"CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o"

# External object files for target compile_MatrixBase_cast
compile_MatrixBase_cast_EXTERNAL_OBJECTS =

doc/snippets/compile_MatrixBase_cast: doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o
doc/snippets/compile_MatrixBase_cast: doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/build.make
doc/snippets/compile_MatrixBase_cast: doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable compile_MatrixBase_cast"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_MatrixBase_cast.dir/link.txt --verbose=$(VERBOSE)
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets && ./compile_MatrixBase_cast >/gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets/MatrixBase_cast.out

# Rule to build all files generated by this target.
doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/build: doc/snippets/compile_MatrixBase_cast
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/build

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/requires: doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/compile_MatrixBase_cast.cpp.o.requires
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/requires

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/clean:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_MatrixBase_cast.dir/cmake_clean.cmake
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/clean

doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/depend:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/home/pgersberg/libeigen /gpfs/home/pgersberg/libeigen/doc/snippets /gpfs/home/pgersberg/libeigen/buikd_dir /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets /gpfs/home/pgersberg/libeigen/buikd_dir/doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_cast.dir/depend

