# ---------------------------------------------------------
# PGLdllTest — Executable that links against PGLdll
# ---------------------------------------------------------

# --- Sources ---------------------------------------------------------------
# The executable will inherit the compiled code and modules via target_link_libraries.
set(SRC_TEST
    dllTestMain.f90
)

# --- Target ----------------------------------------------------------------

add_executable(PGLdllTest ${SRC_TEST})

# --- Link against the PGLdll target ----------------------------------------

target_link_libraries(PGLdllTest PRIVATE PGLdll)

# --- Include PGLdll module directory ---------------------------------------

# Modern CMake Note: If PGLdll defines Fortran_MODULE_DIRECTORY,
# target_link_libraries often propagates this automatically,
# but keeping this explicit ensures Ninja knows exactly where to look.
target_include_directories(PGLdllTest PRIVATE
    $<TARGET_PROPERTY:PGLdll,Fortran_MODULE_DIRECTORY>
)

# --- Intel Fortran flags (Windows) -----------------------------------------

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    target_compile_options(PGLdllTest PRIVATE)
endif()

# --- Debugger working directory --------------------------------------------

set_property(TARGET PGLdllTest PROPERTY
    DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
