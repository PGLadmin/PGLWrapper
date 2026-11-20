@echo off
REM --- Load Intel oneAPI environment ---
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"

if not exist "%VSWHERE%" (
    echo Visual Studio Installer's vswhere.exe was not found. Confirm if Visual Studio is installed.
    exit /b 1
)

set "VS_INSTALL_DIR="
set "VS_MAJOR="
set "VS_VERSION="

for /f "tokens=1 delims=." %%v in (
    '"%VSWHERE%" -latest -products * -property installationVersion'
) do set "VS_MAJOR=%%v"

for /f "tokens=*" %%i in (
    '"%VSWHERE%" -latest -products * -property installationPath'
) do set "VS_INSTALL_DIR=%%i"

REM Require Visual Studio 2019 (16.x) or newer
if not defined VS_MAJOR (
    echo No supported Visual Studio installation found.
    exit /b 1
)

if %VS_MAJOR% LSS 16 (
    echo Visual Studio 2019 or newer is required.
    exit /b 1
)

if %VS_MAJOR% EQU 16 (
    set "VS2019INSTALLDIR=%VS_INSTALL_DIR%"
    set "VS_VERSION=vs2019"
) else if %VS_MAJOR% EQU 17 (
    set "VS2022INSTALLDIR=%VS_INSTALL_DIR%"
    set "VS_VERSION=vs2022"
) else if %VS_MAJOR% EQU 18 (
    set "VS2026INSTALLDIR=%VS_INSTALL_DIR%"
    set "VS_VERSION=vs2026"
) else (
    echo Unsupported Visual Studio version %VS_MAJOR%.
    exit /b 1
)

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 %VS_VERSION%
if errorlevel 1 exit /b %errorlevel%

REM --- Change directory to the folder containing this .bat file ---
cd /d "%~dp0"

REM --- Launch VS Code in this folder ---
code .
