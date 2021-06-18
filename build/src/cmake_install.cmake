# Install script for directory: /mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/mass_borrow_mod/cmake_install.cmake")
  include("/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/output_mod/cmake_install.cmake")
  include("/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/utils_mod/cmake_install.cmake")
  include("/mnt/c/Users/Owner/Education and Research/Graduate School/2020 - 2021 - Third Year/2020 - 2021 - Other/2021 - Summer Practicum/vert_remap/build/src/vertremap_mod/cmake_install.cmake")

endif()

