# ---------------------------------------------------------
#
# Root folder CMakeLists.txt
#
# --------------------------------------------------------

cmake_minimum_required(VERSION 3.30)

# ---------------------------------------------------------
# Top-level meta-project
# ---------------------------------------------------------
project(MDNAproject LANGUAGES Fortran)

# ---------------------------------------------------------
# Global settings (safe defaults)
# ---------------------------------------------------------

# --- Centralize all build outputs into the root binary directory ---
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Track and link the Fortran runtime dynamically across all subprojects by default.
# Override explicitly for specific stand-alone projects using release compile flags within the project CMakeLists.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")

# Ensure all Fortran modules go into a predictable place
# Each target can override this if needed
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Position-independent code is generally safe for DLLs
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---------------------------------------------------------
# Add subprojects
# ---------------------------------------------------------
# Each subproject keeps its own flags and logic
add_subdirectory(Main_PGLdll)
add_subdirectory(Main_PGLdllTest)
add_subdirectory(Main_PGLEOS)

set(FORTLS_SOURCES "")

#-----------------------------------------------------------
# Collect file list for fortls to scan
# Recursively collect all targets from this directory and subdirectories
# ---------------------------------------------------------
function(collect_targets outVar dir)
    # Get targets in this directory
    get_property(localTargets DIRECTORY "${dir}" PROPERTY BUILDSYSTEM_TARGETS)
    set(allTargets ${localTargets})

    # Recurse into subdirectories
    get_property(subdirs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES)
    foreach(sub ${subdirs})
        collect_targets(subTargets "${sub}")
        list(APPEND allTargets ${subTargets})
    endforeach()

    set(${outVar} "${allTargets}" PARENT_SCOPE)
endfunction()

list(REMOVE_DUPLICATES FORTLS_SOURCES)


collect_targets(ALL_TARGETS "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Collected targets: ${ALL_TARGETS}")

foreach(var
    PGLDLL_SOURCES
    PGLDLLTEST_SOURCES
    PGLEOS_SOURCES
    GAMMAPA_SOURCES
)
    if(DEFINED ${var})
        list(APPEND FORTLS_SOURCES ${var})
        list(APPEND FORTLS_SOURCES ${${var}})
    endif()
endforeach()

message(STATUS "FORTLS_SOURCES: ${FORTLS_SOURCES}")

file(WRITE ${CMAKE_BINARY_DIR}/fortls_sources.txt "")

foreach(s ${FORTLS_SOURCES})
    # Normalize path by removing leading ../
    string(REGEX REPLACE "^\\.\\./" "" clean "${s}")
    # Modify the string ",\n" for the desired concatention
    file(APPEND ${CMAKE_BINARY_DIR}/fortls_sources.txt "${clean},\n")
endforeach()

file(COPY
    "${CMAKE_BINARY_DIR}/fortls_sources.txt"
    DESTINATION "${CMAKE_SOURCE_DIR}"
)

# Create portable path
file(COPY
    "${CMAKE_BINARY_DIR}/fortls_sources.txt"
    DESTINATION "${CMAKE_SOURCE_DIR}"
)
