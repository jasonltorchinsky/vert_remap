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
CMAKE_SOURCE_DIR = "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build"

# Include any dependencies generated for this target.
include src/vertremap_mod/CMakeFiles/vertremap_mod.dir/depend.make

# Include the progress variables for this target.
include src/vertremap_mod/CMakeFiles/vertremap_mod.dir/progress.make

# Include the compile flags for this target's objects.
include src/vertremap_mod/CMakeFiles/vertremap_mod.dir/flags.make

src/vertremap_mod/CMakeFiles/vertremap_mod.dir/vertremap_base.F90.o: src/vertremap_mod/CMakeFiles/vertremap_mod.dir/flags.make
src/vertremap_mod/CMakeFiles/vertremap_mod.dir/vertremap_base.F90.o: ../src/vertremap_mod/vertremap_base.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/vertremap_mod/CMakeFiles/vertremap_mod.dir/vertremap_base.F90.o"
	cd "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod" && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/src/vertremap_mod/vertremap_base.F90" -o CMakeFiles/vertremap_mod.dir/vertremap_base.F90.o

src/vertremap_mod/CMakeFiles/vertremap_mod.dir/vertremap_base.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/vertremap_mod.dir/vertremap_base.F90.i"
	cd "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod" && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/src/vertremap_mod/vertremap_base.F90" > CMakeFiles/vertremap_mod.dir/vertremap_base.F90.i

src/vertremap_mod/CMakeFiles/vertremap_mod.dir/vertremap_base.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/vertremap_mod.dir/vertremap_base.F90.s"
	cd "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod" && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/src/vertremap_mod/vertremap_base.F90" -o CMakeFiles/vertremap_mod.dir/vertremap_base.F90.s

# Object files for target vertremap_mod
vertremap_mod_OBJECTS = \
"CMakeFiles/vertremap_mod.dir/vertremap_base.F90.o"

# External object files for target vertremap_mod
vertremap_mod_EXTERNAL_OBJECTS =

src/vertremap_mod/libvertremap_mod.a: src/vertremap_mod/CMakeFiles/vertremap_mod.dir/vertremap_base.F90.o
src/vertremap_mod/libvertremap_mod.a: src/vertremap_mod/CMakeFiles/vertremap_mod.dir/build.make
src/vertremap_mod/libvertremap_mod.a: src/vertremap_mod/CMakeFiles/vertremap_mod.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran static library libvertremap_mod.a"
	cd "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod" && $(CMAKE_COMMAND) -P CMakeFiles/vertremap_mod.dir/cmake_clean_target.cmake
	cd "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vertremap_mod.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/vertremap_mod/CMakeFiles/vertremap_mod.dir/build: src/vertremap_mod/libvertremap_mod.a

.PHONY : src/vertremap_mod/CMakeFiles/vertremap_mod.dir/build

src/vertremap_mod/CMakeFiles/vertremap_mod.dir/clean:
	cd "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod" && $(CMAKE_COMMAND) -P CMakeFiles/vertremap_mod.dir/cmake_clean.cmake
.PHONY : src/vertremap_mod/CMakeFiles/vertremap_mod.dir/clean

src/vertremap_mod/CMakeFiles/vertremap_mod.dir/depend:
	cd "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap" "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/src/vertremap_mod" "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build" "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod" "/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod/CMakeFiles/vertremap_mod.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : src/vertremap_mod/CMakeFiles/vertremap_mod.dir/depend

