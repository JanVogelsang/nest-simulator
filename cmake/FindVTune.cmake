# FindVTune.cmake
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

# - Find VTune header and library
#
# This module defines
#  VTUNE_FOUND, if false, do not try to use VTune.
#  VTUNE_INCLUDE, where to find header files.
#  VTUNE_LIBRARIES, the libraries to link against to use VTune.
#
# As a hint allows VTUNE_ROOT_DIR.

# set( VTUNE_ROOT_DIR "/opt/intel/oneapi/vtune/2025.0/" )

# find include dir
find_path(VTUNE_INCLUDE
        NAMES ittnotify.h
        HINTS ${VTUNE_ROOT_DIR}/sdk/include
)

# find lib dir
find_library(VTUNE_LIBRARIES
        NAMES ittnotify
        HINTS ${VTUNE_ROOT_DIR}/sdk/lib64
)

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( VTune
FOUND_VAR
VTUNE_FOUND
REQUIRED_VARS
VTUNE_INCLUDE
VTUNE_LIBRARIES
)

mark_as_advanced( VTUNE_FOUND VTUNE_INCLUDE VTUNE_LIBRARIES )
