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
include doc/examples/CMakeFiles/Cwise_erfc.dir/depend.make

# Include the progress variables for this target.
include doc/examples/CMakeFiles/Cwise_erfc.dir/progress.make

# Include the compile flags for this target's objects.
include doc/examples/CMakeFiles/Cwise_erfc.dir/flags.make

doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o: doc/examples/CMakeFiles/Cwise_erfc.dir/flags.make
doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o: ../doc/examples/Cwise_erfc.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /gpfs/home/pgersberg/libeigen/buikd_dir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o -c /gpfs/home/pgersberg/libeigen/doc/examples/Cwise_erfc.cpp

doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.i"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /gpfs/home/pgersberg/libeigen/doc/examples/Cwise_erfc.cpp > CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.i

doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.s"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /gpfs/home/pgersberg/libeigen/doc/examples/Cwise_erfc.cpp -o CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.s

doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.requires:
.PHONY : doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.requires

doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.provides: doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.requires
	$(MAKE) -f doc/examples/CMakeFiles/Cwise_erfc.dir/build.make doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.provides.build
.PHONY : doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.provides

doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.provides.build: doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o

# Object files for target Cwise_erfc
Cwise_erfc_OBJECTS = \
"CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o"

# External object files for target Cwise_erfc
Cwise_erfc_EXTERNAL_OBJECTS =

doc/examples/Cwise_erfc: doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o
doc/examples/Cwise_erfc: doc/examples/CMakeFiles/Cwise_erfc.dir/build.make
doc/examples/Cwise_erfc: doc/examples/CMakeFiles/Cwise_erfc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable Cwise_erfc"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Cwise_erfc.dir/link.txt --verbose=$(VERBOSE)
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && ./Cwise_erfc >/gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples/Cwise_erfc.out

# Rule to build all files generated by this target.
doc/examples/CMakeFiles/Cwise_erfc.dir/build: doc/examples/Cwise_erfc
.PHONY : doc/examples/CMakeFiles/Cwise_erfc.dir/build

doc/examples/CMakeFiles/Cwise_erfc.dir/requires: doc/examples/CMakeFiles/Cwise_erfc.dir/Cwise_erfc.cpp.o.requires
.PHONY : doc/examples/CMakeFiles/Cwise_erfc.dir/requires

doc/examples/CMakeFiles/Cwise_erfc.dir/clean:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && $(CMAKE_COMMAND) -P CMakeFiles/Cwise_erfc.dir/cmake_clean.cmake
.PHONY : doc/examples/CMakeFiles/Cwise_erfc.dir/clean

doc/examples/CMakeFiles/Cwise_erfc.dir/depend:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/home/pgersberg/libeigen /gpfs/home/pgersberg/libeigen/doc/examples /gpfs/home/pgersberg/libeigen/buikd_dir /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples/CMakeFiles/Cwise_erfc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/examples/CMakeFiles/Cwise_erfc.dir/depend

