# ---------------------------------------------------------
# PGLEOS — Executable
# ---------------------------------------------------------

# --- Explicit ordered source lists (order matters) ----------------------------

set(SRC_UA
    ../SourceCodeUA/0ModuleUAData.f90
    ../SourceCodeUA/Wertheim23a.f90
    ../SourceCodeUA/FuEsd23a.f90
    ../SourceCodeUA/FuSpeadmd23a.f90
    ../SourceCodeUA/FuFloryWert2a.f90
    ../SourceCodeUA/FugiLSG-MEM2.f90
    ../SourceCodeUA/FugiNRTL.f90
    ../SourceCodeUA/FugiPr2a.f90
    ../SourceCodeUA/PGLFugiDefs.f90
    ../SourceCodeUA/FugiPrTc.f90
    ../SourceCodeUA/fugiPRtcWS.f90
    ../SourceCodeUA/fugiPRWS1a.for
    ../SourceCodeUA/MEM2.f90
    ../SourceCodeUA/UaDataTools.f90
    ../SourceCodeUA/UaEosTools.f90
)

set(SRC_2A
    ../SourceCode2a/PcSaftEosGross.f90
    ../SourceCode2a/PrLorraine.f90
)

set(SRC_FS
    ../FlashSubs/CR2e.f90
    ../FlashSubs/Flashsubs1b.f90
)

set(SRC_MIN
    ../MinPack/LmdifEzCov2.for
    ../MinPack/MinTools.for
)

set(SRC_LOCAL
    MainMenu2d.f90
    MenuItems2c.f90
)

# --- Target -------------------------------------------------------------------

add_executable(PGLEOS
    ${SRC_MIN}
    ${SRC_UA}
    ${SRC_2A}
    ${SRC_FS}
	${SRC_LOCAL}
)

# --- Intel Fortran flags (Windows) -------------------------------------------

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    # For Release/RelWithDebInfo, force the standalone executable to bake in the runtime
    target_compile_options(PGLEOS PRIVATE
        $<$<NOT:$<CONFIG:Debug>>:/libs:static>
    )
endif()

# --- Debugger working directory ----------------------------------------------

# This forces Visual Studio to run from the root repository directory
set_property(TARGET PGLEOS PROPERTY
    DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)

# This forces VS Code (CMake Tools) to also run from the root repository directory
set_property(TARGET PGLEOS PROPERTY
    VSCODE_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)

# --- Fixed-form handling for .for and .f files --------------------------------

# --- Fixed-form handling for .for files --------------------------------------
foreach(src ${SRC_MIN} ${SRC_UA} ${SRC_2A} ${SRC_FS})
    string(TOLOWER "${src}" src_lower)
    if(src_lower MATCHES "\\.(for|f)$")
        set_source_files_properties(${src} PROPERTIES
            COMPILE_OPTIONS "/fixed;/warn:general;/Qdiag-disable:8290;/Qdiag-error-limit:5"
        )
    endif()
endforeach()
