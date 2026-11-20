# ---------------------------------------------------------
# PGLdllTest — Executable that links against PGLdll
# ---------------------------------------------------------

# --- Sources ---------------------------------------------------------------
# The executable will inherit the compiled code and modules via target_link_libraries.
set(PGLDLLTEST_SOURCES_LCL
    ../Main_PGLdllTest/dllTestMain.f90
)

# ---------------------------------------------------------
# Export source list to parent scope for fortls
# ---------------------------------------------------------
set(PGLDLLTEST_SOURCES
    ${PGLDLLTEST_SOURCES_LCL}
    PARENT_SCOPE
)

# --- Target ----------------------------------------------------------------

add_executable(PGLdllTest ${PGLDLLTEST_SOURCES_LCL})

# --- Link against the PGLdll target ----------------------------------------

target_link_libraries(PGLdllTest PRIVATE PGLdll)


# --- Intel Fortran flags (Windows) -----------------------------------------

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    target_compile_options(PGLdllTest PRIVATE)
endif()

# --- Debugger working directory --------------------------------------------

set_property(TARGET PGLdllTest PROPERTY
    DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
