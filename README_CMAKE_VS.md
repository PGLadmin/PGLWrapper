# Visual Studio 2022 Collaborator CMake Setup Guide (Internal Documentation)

This internal guide ensures collaborators can reliably configure, build, and debug the MDNA project using Visual Studio 2022 with CMake and Intel ifx.

To encourage code development and refinement, compile warnings are issued applying the standard for f18. Only in cases where the standard prevents compiling is this omitted.

## Table of Contents
[1. Required Software and Installer Configuration](#1-required-software-and-installer-configuration)<br>
    1.1 Visual Studio 2022 Installation and Configuration<br>
    1.2 Required Workload<br>
    1.3 Required Individual Components<br>

[2. Intel oneAPI Requirements](#2-intel-oneapi-requirements)

[3. Required Environment Settings in CMakePresets.json](#3-required-environment-settings-in-cmakepresetsjson)

[4. Opening the Project Correctly](#4-opening-the-project-correctly)

[5. Verifying VS CMake Preset Integration](#5-verifying-vs-cmake-preset-integration)

[6. Fixing Preset Detection Issues](#6-fixing-preset-detection-issues)

[7. Building and Debugging Targets](#7-building-and-debugging-targets)

[8. Recommended Workflow](#8-recommended-workflow)

[Appendix: Common Errors and Fixes](#appendix-common-errors-and-fixes)


## 1. Required Software and Installer Configuration

1.1 Visual Studio 2022 Installation and Configuration

Install Visual Studio 2022 (Community, Professional, or Enterprise).

### 1.2 Required Workload

Under Workloads:

Desktop development with C++ ✔

> This workload provides:

* MSVC toolchain

* CMake integration

* Ninja build support

* Native debugging tools

### 1.3 Required Individual Components

Under Individual components, ensure the following are selected:

CMake & Build Tools

* C++ CMake tools for Windows ✔

* CMake ✔

* Ninja ✔

* MSBuild support for CMake ✔

* Windows 10/11 SDK ✔

Debugging Tools

* Just-In-Time Debugger ✔

* Native Debugging Tools ✔

These components ensure Visual Studio can configure, build, and debug CMake-based Fortran projects.

#### 1.4 Configuration

After installation start VS and select the option to 'Continue without code'.
* Tools>Options...
* CMake>General
* Select to 'Always use CMake Presets'
* Click OK and you may close VS

This selection prevents VS from falling back to a 32-bit configuration which results in build failures. The CMakePresets.json provided with this project enforce 64-bit builds when loaded properly.

[Return to TOC](#table-of-contents)

## 2. Intel oneAPI Requirements

Install the following Intel toolkits:

Version 2024

>Intel oneAPI Base Toolkit

>Intel HPC Toolkit

Version 2025

>Intel oneAPI HPC Toolkit

These provide:

* ifx (Intel Fortran compiler)

* Intel runtime libraries

* setvars.bat environment setup script

[Return to TOC](#table-of-contents)

## 3. Required Environment Settings in CMakePresets.json

These settings are already configured in the provided CMakePresets.json. The settings ensure Visual Studio and VS Code can:

* Detect the Intel Fortran compiler and compile 64-bit, overriding the default 32-bit configuration of Visual Studio

* Run compiler tests during CMake auto-configuration when the project loads

* Load the correct runtime libraries for debugging

* Link to the PGLdll library for that specific project

[Return to TOC](#table-of-contents)

## 4. Opening the Project Correctly

Collaborators must open the top-level repo folder, not a subfolder, nor a .sln nor .vfproj file. The CMakePresets.json is only in the top repo folder and is inherited by the subfolders.

### Visual Studio Recommended Method for First-Time Project Loading.

On the welcome screen, bypass all options for opening Solutions or Projects, or creating new items.

Use the link to 'Continue without code'.

* File → Open → Folder…
* Select: MDNAproject/

Alternative (forces CMake mode)

* File → Open → CMake…
* Select: MDNAproject/CMakeLists.txt

This ensures Visual Studio enters CMake Project Mode, not Workspace Mode.

Be patient when the project loads because CMake will build the Cache if necessary and write the results to the 'Output' window. CMake will also create the 'out' folder if it does not exist in the repo root. The 'out' folder contains the CMake workspace. The VS output window will report CMake success or errors if present.

### Visual Studio Method for Subsequent Loading

You may use the option above, but you also may use the welcome screen options to open the project folder for the root of the repo. Be patient if CMake recreates the cache.

[Return to TOC](#table-of-contents)

## 5. Verifying VS CMake Preset Integration

When Visual Studio is correctly using the CMake presets, the following UI elements will appear:

Ribbon

> Compile options will list: debug | relwithdebinfo | release

> If instead the compile options lists: x64 Debug | x64 Release | x86 Debug | x86 Release, then Visual Studio is not using the CMake presets and corrective action is needed.

Project Menu will list

> Configure MDNA Project

> Configure MDNA Project with CMake Debugger

> Delete Cache and Reconfigure

> View CMakeCache.txt

> Edit Cmake Presets for MDNAProject

> CMake Workspace Settings

Solution Explorer Window

>The solution explorer defaults to 'Folder View'. Right-click the top folder MDNA Project and select the option for CMake Targets View. Use the accordion links to find the source code used for each project. Notice that PGLdllTest does not duplicate source code in PGLdll. As with the Folder View, files are opened in the editor by double-clicking.

The toolbar should show a startup item dropdown listing: PGLdllTest | PGLEOS. Note that PGLdll is a library and thus not a startup item.

[Return to TOC](#table-of-contents)

## 6. Fixing Preset Detection Issues

If the CMake Cache fails the automatic tests when it builds the cache, you will find an option in the Targets View at the top of the Solution Explorer to Debug. Often that leads to the point of the tests where CMake crashed but not the cause of the crash.

If Visual Studio Menus/Ribbon show:

* Missing CMake-related menu items

* Ribbon shows x64/x86 configurations

Then the errors must be fixed:

> Step 1 — Fix JSON syntax errors. Open CMakePresets.json in the VS editor. Look at the bottom for errors or at the top for "Your settings are not being used due to a parsing error", or "CMakePresets.json could not be parsed". VS will indicate the offending line in the Errors window.

> Even a single missing comma invalidates the entire preset file. Fix the errors. Often the CMake cache can be Deleted and Reconfigured from the project window, but VS copies the CMake cache for its own use, and sometimes the VS cache does not refresh so if the cache does not rebuild, proceeding with steps 2 and 3 are recommended.

Step 2 — **Close VS completely.** Locate the normally hidden .vs folder in the repo root. Delete the .vs folder. This clears stale cached configuration. Also, the 'out' folder can be deleted and CMake will regenerate it.

Step 3 — Restart VS and reopen the folder. Visual Studio will re-detect the CMake project and load the presets.

If necessary, use Detailed Debugging -- By default, CMake 3.30 distributed with VS 2022 has detailed logging disabled. To create the log,

> 1. Use the Project Menu to open CMakeCache.txt
> 2. Search for (part of) CMAKE_VERBOSE_MAKEFILE:BOOL=FALSE and set to TRUE
> 3. Save the file and from the 'Project' menu delete the cache and reconfigure.
> 4. Browse in File Explorer to out\build\debug\CMakeFiles and open CMakeConfigureLog.yaml in a text editor and search the bottom portion for the fatal error.

[Return to TOC](#table-of-contents)

## 7. Building and Debugging Targets

Many of the configuration tools within VS are NOT used. Using additional VS configuration tools beyond those discussed here may confuse CMake.

### Building a Target

**The !DEC$ export for masterDir,PGLinputDir must be toggled ON when compiling PGLDll and toggled OFF when compiling PGLEOS.** Also, some linker warnings are disabled in the .vfproj files. This inconvenience will be eliminated in the future.

* Select the desired compile configuration in the Ribbond dropdown
* In the Solution Explorer, open CMake Targets View
* Right-click a target (e.g., PGLdllTest)
* Select Build

### Debugging a Target

* In the CMake Targets view, open the source code and insert the desired breakpoint(s).
* In the toolbar, use the startup button dropdown to select the desired project and then click the green arrow.
* Right-clicking the project target in the Solution Targets View also provides an option for debugging. (Do not use the right-click option to create a debug configuration because that creates a VS launch.vs.json that will differ from the configurations in CMakePresets.json and cause confusion and may not load the 64 bit environment. The launch.vs.json can be deleted from the .vs folder if it is created.)

#### Debugging the DLL

* Set breakpoints inside PGLdll source files
* Set PGLdllTest as the startup item
* Debug normally

[Return to TOC](#table-of-contents)

## 8. Recommended Workflow

### New Project Creation

* The unique main program should be placed in a folder using the naming convention Main_\<proj name\>. Shared files can be placed in other shared folders.
* Copy an existing CMakeLists.txt as a template and paste into the project folder.
* After files are created, open a powershell and cd to the folder with the source code. Run scan-fortran.ps1 on the folder. The script will recurse through subfolders and report files where modules are defined and common blocks are used. An output file is created in addition to the display.
* Repeat as necessary on all folders that contain source code used in the project.
* Manually order the source code files into the CMakeLists.txt for proper compile order so that modules are compiled before code that uses the modules.
* Files that are static such as numerical recipes can be placed near the top of the list. CMake will typically compile only the files needed for
* Files that define subroutines are best placed before files that call subroutines. The main program is often compiled last.
* Link errors of missing symbols will clarify where the manual ordering needs refinement. Clear the project and rebuild when you rearrange files.

### Development - debug builds

Build individual targets and debug symbols will be included.  Debug using breakpoints. Select to build the project instead of selecting individual files to compile.

* Set the 'debug' build preset from the dropdown in the ribbon. See the VS Code README for details on the build configurations, and finding the files. The 'out' folder structure is the same as using CMake in VS Code.

### Performance Testing - relwithdebinfo

* Set the 'relwithdebinfo' preset. This build will implement optimization for speed but includes debug symbols. Benchmarks developed with a debug configuration should be repeated to assure that the optimization does not break loops or dependencies.

### Final Builds - release

* Use the 'release' preset

This build implements optimization and hides debug symbols.

[Return to TOC](#table-of-contents)

## Special Note - Switching between CMake and .vfproj builds

If you switch between a CMake and traditional .sln/.vfproj build, completely close Visual Studio and delete the hidden .vs folder from the repo folder. Then reopen Visual Studio and the desired cache will be rebuilt. If Visual Studio is acting weird, perform this 'reset'.

## Appendix: Common Errors and Fixes

### A. Visual Studio shows x64/x86 Debug/Release in the CMake environment.

Cause: Visual Studio ribbon is not using the CMake presets.

Fix:

> Open Project → CMake Settings and select a preset

> If missing, **close VS**, delete .vs (and optionally the 'out' folder) and reopen the folder

### B. "Your settings are not being used due to a parsing error"

Cause: Invalid JSON in CMakePresets.json.

Fix:

> Correct the JSON syntax

> As soon as the file is saved, the cache should rebuild.

> If necessary, **Close VS**, delete .vs (and optionally the 'out' folder) and reopen the folder

### C. Intel Fortran compiler not detected

Cause: Missing environment variables.

Fix:

> Ensure your preset includes: "FC": "ifx" or use the full path to ifx.exe

### D. CMake configuration fails with missing DLLs

Cause: Intel runtime not on PATH.

Fix:

> Ensure your preset includes:

> "PATH": "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin;$env{PATH}"

### E. CMake Targets View is empty

Cause: Visual Studio opened the folder in Workspace (.sln, .vfproj) Mode.

Fix:

>Close VS, delete the .vs folder.

> Use File → Open → CMake… Select the top-level CMakeLists.txt<br>
> or Open the Folder

### F. Breakpoints show "unresolved"

Cause: Using Release preset.

Fix: Use debug or relwithdebinfo presets.

### G. Visual Studio fails to configure after Intel install

Cause: Stale cached configuration.

Fix: Close VS and delete the .vs folder.

### H. VS continues to fail to implement the environment variables for FORTRAN

Cause: Various

Fix:

This is the 'nuclear' option and should assure that all Intel oneAPI settings are inherited by VS.

* Click Windows Start Menu
* Search for: Intel oneAPI command prompt for Intel 64 for Visual Studio 2022
* cd to the repo root folder
* enter 'devenv .' This loads VS from within the Intel environment, so all FORTRAN environment variables should exist within VS.

[Return to TOC](#table-of-contents)
