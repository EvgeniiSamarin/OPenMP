# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/evgeniisamarin/CLionProjects/openMP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/evgeniisamarin/CLionProjects/openMP/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/MyOpenCVTest.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/MyOpenCVTest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MyOpenCVTest.dir/flags.make

CMakeFiles/MyOpenCVTest.dir/main.cpp.o: CMakeFiles/MyOpenCVTest.dir/flags.make
CMakeFiles/MyOpenCVTest.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/evgeniisamarin/CLionProjects/openMP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MyOpenCVTest.dir/main.cpp.o"
	/usr/local/Cellar/llvm/6.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MyOpenCVTest.dir/main.cpp.o -c /Users/evgeniisamarin/CLionProjects/openMP/main.cpp

CMakeFiles/MyOpenCVTest.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MyOpenCVTest.dir/main.cpp.i"
	/usr/local/Cellar/llvm/6.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/evgeniisamarin/CLionProjects/openMP/main.cpp > CMakeFiles/MyOpenCVTest.dir/main.cpp.i

CMakeFiles/MyOpenCVTest.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MyOpenCVTest.dir/main.cpp.s"
	/usr/local/Cellar/llvm/6.0.0/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/evgeniisamarin/CLionProjects/openMP/main.cpp -o CMakeFiles/MyOpenCVTest.dir/main.cpp.s

# Object files for target MyOpenCVTest
MyOpenCVTest_OBJECTS = \
"CMakeFiles/MyOpenCVTest.dir/main.cpp.o"

# External object files for target MyOpenCVTest
MyOpenCVTest_EXTERNAL_OBJECTS =

MyOpenCVTest: CMakeFiles/MyOpenCVTest.dir/main.cpp.o
MyOpenCVTest: CMakeFiles/MyOpenCVTest.dir/build.make
MyOpenCVTest: CMakeFiles/MyOpenCVTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/evgeniisamarin/CLionProjects/openMP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MyOpenCVTest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MyOpenCVTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MyOpenCVTest.dir/build: MyOpenCVTest
.PHONY : CMakeFiles/MyOpenCVTest.dir/build

CMakeFiles/MyOpenCVTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MyOpenCVTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MyOpenCVTest.dir/clean

CMakeFiles/MyOpenCVTest.dir/depend:
	cd /Users/evgeniisamarin/CLionProjects/openMP/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/evgeniisamarin/CLionProjects/openMP /Users/evgeniisamarin/CLionProjects/openMP /Users/evgeniisamarin/CLionProjects/openMP/cmake-build-debug /Users/evgeniisamarin/CLionProjects/openMP/cmake-build-debug /Users/evgeniisamarin/CLionProjects/openMP/cmake-build-debug/CMakeFiles/MyOpenCVTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MyOpenCVTest.dir/depend

