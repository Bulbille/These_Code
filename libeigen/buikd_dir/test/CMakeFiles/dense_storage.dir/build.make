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
include test/CMakeFiles/dense_storage.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/dense_storage.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/dense_storage.dir/flags.make

test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o: test/CMakeFiles/dense_storage.dir/flags.make
test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o: ../test/dense_storage.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /gpfs/home/pgersberg/libeigen/buikd_dir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dense_storage.dir/dense_storage.cpp.o -c /gpfs/home/pgersberg/libeigen/test/dense_storage.cpp

test/CMakeFiles/dense_storage.dir/dense_storage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dense_storage.dir/dense_storage.cpp.i"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /gpfs/home/pgersberg/libeigen/test/dense_storage.cpp > CMakeFiles/dense_storage.dir/dense_storage.cpp.i

test/CMakeFiles/dense_storage.dir/dense_storage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dense_storage.dir/dense_storage.cpp.s"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /gpfs/home/pgersberg/libeigen/test/dense_storage.cpp -o CMakeFiles/dense_storage.dir/dense_storage.cpp.s

test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.requires:
.PHONY : test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.requires

test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.provides: test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/dense_storage.dir/build.make test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.provides.build
.PHONY : test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.provides

test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.provides.build: test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o

# Object files for target dense_storage
dense_storage_OBJECTS = \
"CMakeFiles/dense_storage.dir/dense_storage.cpp.o"

# External object files for target dense_storage
dense_storage_EXTERNAL_OBJECTS =

test/dense_storage: test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o
test/dense_storage: test/CMakeFiles/dense_storage.dir/build.make
test/dense_storage: test/CMakeFiles/dense_storage.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable dense_storage"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dense_storage.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/dense_storage.dir/build: test/dense_storage
.PHONY : test/CMakeFiles/dense_storage.dir/build

test/CMakeFiles/dense_storage.dir/requires: test/CMakeFiles/dense_storage.dir/dense_storage.cpp.o.requires
.PHONY : test/CMakeFiles/dense_storage.dir/requires

test/CMakeFiles/dense_storage.dir/clean:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/test && $(CMAKE_COMMAND) -P CMakeFiles/dense_storage.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/dense_storage.dir/clean

test/CMakeFiles/dense_storage.dir/depend:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/home/pgersberg/libeigen /gpfs/home/pgersberg/libeigen/test /gpfs/home/pgersberg/libeigen/buikd_dir /gpfs/home/pgersberg/libeigen/buikd_dir/test /gpfs/home/pgersberg/libeigen/buikd_dir/test/CMakeFiles/dense_storage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/dense_storage.dir/depend

