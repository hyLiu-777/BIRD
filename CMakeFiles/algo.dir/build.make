# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/haoyu/ppr/BDPush

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/haoyu/ppr/BDPush

# Include any dependencies generated for this target.
include CMakeFiles/algo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/algo.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/algo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/algo.dir/flags.make

CMakeFiles/algo.dir/src/algo.cc.o: CMakeFiles/algo.dir/flags.make
CMakeFiles/algo.dir/src/algo.cc.o: src/algo.cc
CMakeFiles/algo.dir/src/algo.cc.o: CMakeFiles/algo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haoyu/ppr/BDPush/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/algo.dir/src/algo.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/algo.dir/src/algo.cc.o -MF CMakeFiles/algo.dir/src/algo.cc.o.d -o CMakeFiles/algo.dir/src/algo.cc.o -c /home/haoyu/ppr/BDPush/src/algo.cc

CMakeFiles/algo.dir/src/algo.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/algo.dir/src/algo.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haoyu/ppr/BDPush/src/algo.cc > CMakeFiles/algo.dir/src/algo.cc.i

CMakeFiles/algo.dir/src/algo.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/algo.dir/src/algo.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haoyu/ppr/BDPush/src/algo.cc -o CMakeFiles/algo.dir/src/algo.cc.s

# Object files for target algo
algo_OBJECTS = \
"CMakeFiles/algo.dir/src/algo.cc.o"

# External object files for target algo
algo_EXTERNAL_OBJECTS =

build/libalgo.a: CMakeFiles/algo.dir/src/algo.cc.o
build/libalgo.a: CMakeFiles/algo.dir/build.make
build/libalgo.a: CMakeFiles/algo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/haoyu/ppr/BDPush/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library build/libalgo.a"
	$(CMAKE_COMMAND) -P CMakeFiles/algo.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/algo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/algo.dir/build: build/libalgo.a
.PHONY : CMakeFiles/algo.dir/build

CMakeFiles/algo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/algo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/algo.dir/clean

CMakeFiles/algo.dir/depend:
	cd /home/haoyu/ppr/BDPush && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/haoyu/ppr/BDPush /home/haoyu/ppr/BDPush /home/haoyu/ppr/BDPush /home/haoyu/ppr/BDPush /home/haoyu/ppr/BDPush/CMakeFiles/algo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/algo.dir/depend

