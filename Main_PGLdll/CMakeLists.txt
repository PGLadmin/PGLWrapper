# ---------------------------------------------------------
# PGLdll — Shared library project
# ---------------------------------------------------------

# --- Source lists -------------------------------------------------------------

set(SRC_UA
    ../SourceCodeUA/0ModuleUAData.f90
    ../SourceCodeUA/PGLFugiDefs.f90
    ../SourceCodeUA/PGLDLL.f90
    ../SourceCodeUA/UaDataTools.f90
    ../SourceCodeUA/UaEosTools.f90
    ../SourceCodeUA/MEM2.f90
    ../SourceCodeUA/FuEsd23a.f90
    ../SourceCodeUA/FuSpeadmd23a.f90
    ../SourceCodeUA/FugiLSG-MEM2.f90
    ../SourceCodeUA/FugiNRTL.f90
    ../SourceCodeUA/FugiPr2a.f90
    ../SourceCodeUA/FugiPrTc.f90
    ../SourceCodeUA/fugiPRtcWS.f90
    ../SourceCodeUA/fugiPRWS1a.for
    ../SourceCodeUA/Wertheim23a.f90
)

set(SRC_2A
    ../SourceCode2a/PcSaftEosGross.f90
    ../SourceCode2a/PrLorraine.f90
)

set(SRC_LOCAL
    PGLWrapper.f90
)

# --- Target -------------------------------------------------------------------

add_library(PGLdll SHARED
    ${SRC_UA}
    ${SRC_2A}
    ${SRC_LOCAL}
)

# --- Intel Fortran flags (Windows) -------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    target_compile_options(PGLdll PRIVATE)
endif()

# --- Fixed-form handling for .for files --------------------------------------
foreach(src ${SRC_UA})
    string(TOLOWER "${src}" src_lower)
    if(src_lower MATCHES "\\.(for|f)$")
        set_source_files_properties(${src} PROPERTIES
            COMPILE_OPTIONS "/fixed;/warn:general;/Qdiag-disable:8290;/Qdiag-error-limit:5"
        )
    endif()
endforeach()
