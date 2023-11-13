# _FindMetavisionSDK.cmake
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

# - Find Metavision SDK header and library
#
# This module defines
#  METAVISIONSDK_FOUND, if false, do not try to use Metavision SDK.
#  METAVISIONSDK_INCLUDE, where to find includes.
#  METAVISIONSDK_LIBRARIES, the libraries to link against to use Metavision SDK.
#  METAVISIONSDK_BINARIES, where to find binaries.
#  METAVISIONSDK_SOFTWARE_INFO, Metavision executable to retrieve version info.
#  METAVISIONSDK_VERSION, the library version.
#
# Requires METAVISIONSDK_ROOT_DIR to be set.

set( METAVISIONSDK_INCLUDE ${METAVISIONSDK_ROOT_DIR}/include )

find_library( METAVISIONSDK_LIBRARY_CORE
        NAMES metavision_sdk_core
        HINTS ${METAVISIONSDK_ROOT_DIR}/lib
)
find_library( METAVISIONSDK_LIBRARY_DRIVER
        NAMES metavision_sdk_driver
        HINTS ${METAVISIONSDK_ROOT_DIR}/lib
)
set( METAVISIONSDK_LIBRARIES ${METAVISIONSDK_LIBRARY_CORE} ${METAVISIONSDK_LIBRARY_DRIVER} )

set( METAVISIONSDK_BINARIES ${METAVISIONSDK_ROOT_DIR}/bin )

find_program( METAVISIONSDK_SOFTWARE_INFO
    NAMES metavision_software_info
    HINTS ${METAVISIONSDK_ROOT_DIR}/bin
)

execute_process(
    COMMAND ${METAVISIONSDK_SOFTWARE_INFO} --version
    RESULT_VARIABLE RESULT
    OUTPUT_VARIABLE METAVISIONSDK_VAR_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if ( RESULT EQUAL 0 )
    set( METAVISIONSDK_VERSION ${METAVISIONSDK_VAR_OUTPUT} )
endif ()

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( MetavisionSDK
  FOUND_VAR
    METAVISIONSDK_FOUND
  REQUIRED_VARS
    METAVISIONSDK_LIBRARY_CORE
    METAVISIONSDK_LIBRARY_DRIVER
    METAVISIONSDK_INCLUDE
  VERSION_VAR
    METAVISIONSDK_VERSION
  )

mark_as_advanced( METAVISIONSDK_ROOT_DIR METAVISIONSDK_LIBRARY_CORE METAVISIONSDK_LIBRARY_DRIVER METAVISIONSDK_BINARIES )
