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
include doc/examples/CMakeFiles/TemplateKeyword_simple.dir/depend.make

# Include the progress variables for this target.
include doc/examples/CMakeFiles/TemplateKeyword_simple.dir/progress.make

# Include the compile flags for this target's objects.
include doc/examples/CMakeFiles/TemplateKeyword_simple.dir/flags.make

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o: doc/examples/CMakeFiles/TemplateKeyword_simple.dir/flags.make
doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o: ../doc/examples/TemplateKeyword_simple.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /gpfs/home/pgersberg/libeigen/buikd_dir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o -c /gpfs/home/pgersberg/libeigen/doc/examples/TemplateKeyword_simple.cpp

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.i"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /gpfs/home/pgersberg/libeigen/doc/examples/TemplateKeyword_simple.cpp > CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.i

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.s"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /gpfs/home/pgersberg/libeigen/doc/examples/TemplateKeyword_simple.cpp -o CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.s

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.requires:
.PHONY : doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.requires

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.provides: doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.requires
	$(MAKE) -f doc/examples/CMakeFiles/TemplateKeyword_simple.dir/build.make doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.provides.build
.PHONY : doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.provides

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.provides.build: doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o

# Object files for target TemplateKeyword_simple
TemplateKeyword_simple_OBJECTS = \
"CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o"

# External object files for target TemplateKeyword_simple
TemplateKeyword_simple_EXTERNAL_OBJECTS =

doc/examples/TemplateKeyword_simple: doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o
doc/examples/TemplateKeyword_simple: doc/examples/CMakeFiles/TemplateKeyword_simple.dir/build.make
doc/examples/TemplateKeyword_simple: doc/examples/CMakeFiles/TemplateKeyword_simple.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable TemplateKeyword_simple"
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TemplateKeyword_simple.dir/link.txt --verbose=$(VERBOSE)
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && ./TemplateKeyword_simple >/gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples/TemplateKeyword_simple.out

# Rule to build all files generated by this target.
doc/examples/CMakeFiles/TemplateKeyword_simple.dir/build: doc/examples/TemplateKeyword_simple
.PHONY : doc/examples/CMakeFiles/TemplateKeyword_simple.dir/build

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/requires: doc/examples/CMakeFiles/TemplateKeyword_simple.dir/TemplateKeyword_simple.cpp.o.requires
.PHONY : doc/examples/CMakeFiles/TemplateKeyword_simple.dir/requires

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/clean:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples && $(CMAKE_COMMAND) -P CMakeFiles/TemplateKeyword_simple.dir/cmake_clean.cmake
.PHONY : doc/examples/CMakeFiles/TemplateKeyword_simple.dir/clean

doc/examples/CMakeFiles/TemplateKeyword_simple.dir/depend:
	cd /gpfs/home/pgersberg/libeigen/buikd_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/home/pgersberg/libeigen /gpfs/home/pgersberg/libeigen/doc/examples /gpfs/home/pgersberg/libeigen/buikd_dir /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples /gpfs/home/pgersberg/libeigen/buikd_dir/doc/examples/CMakeFiles/TemplateKeyword_simple.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/examples/CMakeFiles/TemplateKeyword_simple.dir/depend

