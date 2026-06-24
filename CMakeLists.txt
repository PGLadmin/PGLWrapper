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

# When you add Main_GammaPA:
# add_subdirectory(Main_GammaPA)
