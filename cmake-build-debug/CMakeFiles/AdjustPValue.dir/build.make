# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/EL/CLionProjects/AdjustPValue

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/EL/CLionProjects/AdjustPValue/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/AdjustPValue.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/AdjustPValue.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/AdjustPValue.dir/flags.make

CMakeFiles/AdjustPValue.dir/main.cpp.o: CMakeFiles/AdjustPValue.dir/flags.make
CMakeFiles/AdjustPValue.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/EL/CLionProjects/AdjustPValue/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/AdjustPValue.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AdjustPValue.dir/main.cpp.o -c /Users/EL/CLionProjects/AdjustPValue/main.cpp

CMakeFiles/AdjustPValue.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AdjustPValue.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/EL/CLionProjects/AdjustPValue/main.cpp > CMakeFiles/AdjustPValue.dir/main.cpp.i

CMakeFiles/AdjustPValue.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AdjustPValue.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/EL/CLionProjects/AdjustPValue/main.cpp -o CMakeFiles/AdjustPValue.dir/main.cpp.s

CMakeFiles/AdjustPValue.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/AdjustPValue.dir/main.cpp.o.requires

CMakeFiles/AdjustPValue.dir/main.cpp.o.provides: CMakeFiles/AdjustPValue.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/AdjustPValue.dir/build.make CMakeFiles/AdjustPValue.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/AdjustPValue.dir/main.cpp.o.provides

CMakeFiles/AdjustPValue.dir/main.cpp.o.provides.build: CMakeFiles/AdjustPValue.dir/main.cpp.o


CMakeFiles/AdjustPValue.dir/InputFile.cpp.o: CMakeFiles/AdjustPValue.dir/flags.make
CMakeFiles/AdjustPValue.dir/InputFile.cpp.o: ../InputFile.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/EL/CLionProjects/AdjustPValue/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/AdjustPValue.dir/InputFile.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AdjustPValue.dir/InputFile.cpp.o -c /Users/EL/CLionProjects/AdjustPValue/InputFile.cpp

CMakeFiles/AdjustPValue.dir/InputFile.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AdjustPValue.dir/InputFile.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/EL/CLionProjects/AdjustPValue/InputFile.cpp > CMakeFiles/AdjustPValue.dir/InputFile.cpp.i

CMakeFiles/AdjustPValue.dir/InputFile.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AdjustPValue.dir/InputFile.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/EL/CLionProjects/AdjustPValue/InputFile.cpp -o CMakeFiles/AdjustPValue.dir/InputFile.cpp.s

CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.requires:

.PHONY : CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.requires

CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.provides: CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.requires
	$(MAKE) -f CMakeFiles/AdjustPValue.dir/build.make CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.provides.build
.PHONY : CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.provides

CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.provides.build: CMakeFiles/AdjustPValue.dir/InputFile.cpp.o


CMakeFiles/AdjustPValue.dir/PValue.cpp.o: CMakeFiles/AdjustPValue.dir/flags.make
CMakeFiles/AdjustPValue.dir/PValue.cpp.o: ../PValue.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/EL/CLionProjects/AdjustPValue/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/AdjustPValue.dir/PValue.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AdjustPValue.dir/PValue.cpp.o -c /Users/EL/CLionProjects/AdjustPValue/PValue.cpp

CMakeFiles/AdjustPValue.dir/PValue.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AdjustPValue.dir/PValue.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/EL/CLionProjects/AdjustPValue/PValue.cpp > CMakeFiles/AdjustPValue.dir/PValue.cpp.i

CMakeFiles/AdjustPValue.dir/PValue.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AdjustPValue.dir/PValue.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/EL/CLionProjects/AdjustPValue/PValue.cpp -o CMakeFiles/AdjustPValue.dir/PValue.cpp.s

CMakeFiles/AdjustPValue.dir/PValue.cpp.o.requires:

.PHONY : CMakeFiles/AdjustPValue.dir/PValue.cpp.o.requires

CMakeFiles/AdjustPValue.dir/PValue.cpp.o.provides: CMakeFiles/AdjustPValue.dir/PValue.cpp.o.requires
	$(MAKE) -f CMakeFiles/AdjustPValue.dir/build.make CMakeFiles/AdjustPValue.dir/PValue.cpp.o.provides.build
.PHONY : CMakeFiles/AdjustPValue.dir/PValue.cpp.o.provides

CMakeFiles/AdjustPValue.dir/PValue.cpp.o.provides.build: CMakeFiles/AdjustPValue.dir/PValue.cpp.o


# Object files for target AdjustPValue
AdjustPValue_OBJECTS = \
"CMakeFiles/AdjustPValue.dir/main.cpp.o" \
"CMakeFiles/AdjustPValue.dir/InputFile.cpp.o" \
"CMakeFiles/AdjustPValue.dir/PValue.cpp.o"

# External object files for target AdjustPValue
AdjustPValue_EXTERNAL_OBJECTS =

AdjustPValue: CMakeFiles/AdjustPValue.dir/main.cpp.o
AdjustPValue: CMakeFiles/AdjustPValue.dir/InputFile.cpp.o
AdjustPValue: CMakeFiles/AdjustPValue.dir/PValue.cpp.o
AdjustPValue: CMakeFiles/AdjustPValue.dir/build.make
AdjustPValue: CMakeFiles/AdjustPValue.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/EL/CLionProjects/AdjustPValue/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable AdjustPValue"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/AdjustPValue.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/AdjustPValue.dir/build: AdjustPValue

.PHONY : CMakeFiles/AdjustPValue.dir/build

CMakeFiles/AdjustPValue.dir/requires: CMakeFiles/AdjustPValue.dir/main.cpp.o.requires
CMakeFiles/AdjustPValue.dir/requires: CMakeFiles/AdjustPValue.dir/InputFile.cpp.o.requires
CMakeFiles/AdjustPValue.dir/requires: CMakeFiles/AdjustPValue.dir/PValue.cpp.o.requires

.PHONY : CMakeFiles/AdjustPValue.dir/requires

CMakeFiles/AdjustPValue.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/AdjustPValue.dir/cmake_clean.cmake
.PHONY : CMakeFiles/AdjustPValue.dir/clean

CMakeFiles/AdjustPValue.dir/depend:
	cd /Users/EL/CLionProjects/AdjustPValue/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/EL/CLionProjects/AdjustPValue /Users/EL/CLionProjects/AdjustPValue /Users/EL/CLionProjects/AdjustPValue/cmake-build-debug /Users/EL/CLionProjects/AdjustPValue/cmake-build-debug /Users/EL/CLionProjects/AdjustPValue/cmake-build-debug/CMakeFiles/AdjustPValue.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/AdjustPValue.dir/depend

