# Onboarding Guide: VS Code Fortran Development Setup

This repository offers an option to use **VS Code**, the **Modern FORTRAN** and **CMake** VS Code extensions, working in a **Conda** environment to manage the `fortls` language server, the Intel compiler as an interactive FORTRAN linter, and visual debugging. In this environment, compiling, linking, debugging is provided using the CMake extension. Follow these steps to set up your VS Code environment. There is a separate .md file for working with Visual Studio with CMake.

# Table of Contents

Click the corresponding link to visit a topic. Links are provided at the end of each section to jump back to the TOC.

[Step 1: Install Visual Studio, Intel Fortran and VS Code](#step-1-install-visual-studio-intel-fortran-and-vs-code)

[Step 2: Install a Python Distribution](#step-2-install-a-python-distribution)

[Step 3: Create the Shared Fortran-Tools Python Environment](#step-3-create-the-shared-fortran-tools-python-environment)

[Step 4: Open the Project in VS Code](#step-4-open-the-project-in-vs-code)

[Step 5.  Introduction to the VS Code Menus and Icons](#step-5-introduction-to-the-vs-code-menus-and-icons)

[Step 6. Understanding Linting and Language Services](#step-6-understanding-linting-and-language-services)

[Step 7. Browsing Projects, Building Projects](#step-7-browsing-projects-building-projects)

[Step 8. Cleaning the projects](#step-8-cleaning-the-projects)

[Step 9. Building a project](#step-9-building-a-project)

[Step 10. What is the working directory and where is the compiled code?](#step-10-what-is-the-working-directory-and-where-is-the-compiled-code)

[Step 11. Starting the Debugger](#step-11-starting-the-debugger)

[Step 12. Opening Source Code and Setting Breakpoints](#step-12-opening-source-code-and-setting-breakpoints)

[13. A Few Hints on Getting Started with VS Code](#13-a-few-hints-on-getting-started-with-vs-code)

[14. Adding Files to a Project](#14-adding-files-to-a-project)

[15. Adjusting Compile Flags](#15-adjusting-compile-flags)

[16. Troubleshooting](#16-troubleshooting)

---
## Step 1. Install Visual Studio, Intel Fortran and VS Code

### Install Visual Studio
The configuration of this repo is intended to be used with either Visual Studio (using the `.sln` and `.vfproj` files or with the CMake files) or Visual Studio Code (using the CMake files). The VS Code workflow relies on the environment variables configured by Intel oneAPI and Visual Studio, so Visual Studio must be installed first.

Visual Studio 2019 or greater is required. Visual Studio 2019 will require an upgrade of CMake, and details are provided below. Visual Studio 2022 is recommended. Version 2026 is fine, but may have more bugs and stability issues until it matures. To download the 2022 version, [download and install Visual Studio 2022](https://visualstudio.microsoft.com/downloads).

> Follow the guide in README_CMAKE_VS.md Section 1-1.3 to install the needed Visual Studio components. The debugger, C++ environment, CMake and the Ninja compiler components are included in the components. The other parts of that file regard using cmake within Visual Studio and are not neede for VS Code.

### Install Intel FORTRAN onAPI

The FORTRAN 2025 version is recommended as being very stable compared to the most recent release. Download and install Intel Fortran oneAPI 2025.3 from the [Intel Archive](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-archive.html). Note that releases are available in the archive only two years after initial release.

### Install VS Code
[Download VS Code](https://code.visualstudio.com/download) and install as a **'User'**, not as an administrator. Run the installation on the user account that will be used for working with this repository.

[Return to TOC](#table-of-contents)

## Step 2. Install a Python Distribution

Python is available via several distributions including Anaconda, Miniconda, CondaForge, MiniForge, and any will work. To check if you have a distribution,
* Look in the start menu for a shortcut to a prompt for one of the distributions. Depending on the installation settings, the configuration of the shortcut may have been skipped, so you still may have an installation even if you don't have an Conda item in the Start Menu.
* If you don't find a Start Menu shortcut, navigate to your user profile directory C:\Users\\<username\> which can also be found by entering %USERPROFILE% in the Windows Explorer file path bar. Look for any of the folders with the distribution names. Looking for a folder is important instead of opening a command window or powershell, because depending on settings during installation the conda commands may not work by default in a command prompt or powershell.
* If you have installed the stand-alone Jupyter Notebook Desktop, it is installed in the hidden AppData folder and relying on that installation is fragile. A separate installation of Conda is recommended. Jupyter notebook can be added to the separate conda installation later and the stand-alone Notebook installation can be removed to save diskspace if desired.

If you do not have a conda installation, you can search the Internet so select the version to install. The 'Mini' versions are smaller and have all the features needed for this project. Instructions are provided here for MiniConda3. All of the files distributed for the project are configured for Miniconda3, but instructions are provided below for customizing the files to work with a different conda distribution and for maintaining the customizations with Git version tracking. The default installation settings are acceptable, though pay attention to installation options if you want to load Python by default in every command or powershell window, or if you want a Start Menu shortcut to open a Python prompt in the (base) environment:

* **Windows Installer:** [Miniconda Windows 64-bit](https://docs.conda.io/en/latest/miniconda.html)
* *Note:* Install it as a **user** on the same account as VS Code. Ensure it installs to the default directory (`%USERPROFILE%\Miniconda3`), for example `C:\Users\<your username>\Miniconda3`. If you install one of the other distributions, the distribution name will appear in place of Miniconda3.

[Return to TOC](#table-of-contents)

## Step 3. Create the Shared Fortran-Tools Python Environment

Open a _Command_ (cmd) terminal and run the following commands to create the fortran-tools needed for this repository (substitute the distribution name in the path if you have installed a different distribution):

1. Open a command window to activate the Python base environment, or you can skip the command if you open a Start Menu Anaconda/Miniconda/CondaForge/Miniforge Prompt which will load the (base) environment.<br>
    ```cmd
    %USERPROFILE%\miniconda3\Scripts\activate.bat
    ```
2. Once you see (base) indicated on the command prompt, create the isolated fortran-tools environment with the FORTRAN language server fortls. Note that you can optionally install in the (base), but keeping the installation isolated is considered best practice to avoid dependency conflicts and assure reliable package updates.<br>
    ```cmd
    conda create --name fortran-tools python fortls -c conda-forge --yes
    ```
3. Activate the new environment and install the fprettify FORTRAN formatter into the fortran-tools environment<br>
    ```cmd
    conda activate fortran-tools
    ```
    The prompt should change to (fortran-tools) and fprettify will be installed in that environment using the preferred pip method.
    ```cmd
    pip install fprettify
    ```
4. If the system is running Visual Studio 2019, the CMake version distributed with VS2019 is too old and ninja can also be updated. One option is to install CMake from the official KitWare distribution, which makes it available outside of fortran-tools, but for this project, it is suitable to run the following command to add a minimum of version 3.30 to the fortran-tools environment. If not already in the fortran-tools environment, run
    ```cmd
    conda activate fortran-tools
    ```
    Then, at the (fortran-tools) prompt, install the newer CMake and ninja
    ```cmd
    conda install -n fortran-tools -c conda-forge "cmake>=3.30" ninja
    ```
    You may need to add the cmake.exe path to VS Code, so determine that path for your installation, often something like %USERPROFILE%\Miniconda3\envs\fortran-tools\Library\bin\cmake.exe. You can determine installation path using 'cmake --version' which should report 3.30 or newer

You may close the Command window. If you want to confirm installation of fortls and fprettify, navigate in Windows File Explorer to %USERPROFILE%\Miniconda3\envs\fortran-tools\Scripts to find fortls.exe, fprettify.exe.

[Return to TOC](#table-of-contents)

## Step 4. Open the Project in VS Code

### A. Customize settings.json if necessary

Navigate to the project root folder. If you have installed a different distribution of conda instead of Miniconda3, open `.vscode/settings.json` and globally replace `Miniconda3` with the path to your distribution. Keep this machine-specific customization local; do not commit it.

Continue with step B, but return here if you discover errors in cmake. If you have VS2019 and updated cmake, and cmake produces and error, you will need to add the path to settings.json or the VS Code user interface.

1. If you don't see the terminal window, use "Cntl + `" to display the terminal frame.
2. Click the tab for 'Output'.
3. Adjust the filter to 'CMake/Build' and scroll to the top.
4. Look for the command that loads cmake. If it is not your newly installed cmake, one option is to add a custom path to settings.json adjusting the path below to your path (include the comma at the end):

   ```cmd
   "cmake.cmakePath": "${env:USERPROFILE}/Miniconda3/envs/fortran-tools/Library/bin/cmake.exe",
   ```
   The other option is

   1. Use "Cntl + ," to open settings.
   2. On the user tab, search for 'cmake.cmakePath'.
   3. Enter the path to your newly installed cmake.exe.

See [13. A Few Hints on Getting Started with VS Code](#13-a-few-hints-on-getting-started-with-vs-code) to merge updates from the repo without overwriting the customizations.
### B. Start VS Code from the Intel oneAPI environment

* **Preferred - Navigate to the repo directory in Command Prompt (cmd) and run `launch_vscode.bat`**. The script uses Visual Studio Installer's `vswhere.exe` to locate the latest installed Visual Studio. Keep the Command Prompt open if you may need to relaunch VS Code from the initialized oneAPI environment.

  **or**

1. Start the Intel oneAPI command prompt for your installed Visual Studio version.
2. Change directory to the root folder of this repo.
3. Launch VS Code from that prompt by entering `code .`.

Use the script or these steps for subsequent sessions so VS Code receives the oneAPI environment variables.

4. If you use other profiles in VSCode, it is recommended to create a profile for this folder so that settings.json is not modified by other extensions (such as Python). See [VS Help on Profiles](https://code.visualstudio.com/docs/configure/profiles).

On the first launch, when prompted, install the recommended extensions:

* **Modern Fortran**

* **CMake (by Microsoft)**. (This is the VS Code interface to the cmake installation, not the cmake application.)

* The **C/C++** extension will be automatically installed as a dependency.

If you miss the installation notification prompt, click the bell icon in the bottom right of VS Code. You can also install extensions by clicking the building blocks button on the far left panel.

* Optionally, install the **Bookmarks** extension by Alessandro Fragnani. This extension permits you to bookmark lines of code to return to work in progress. One caveat with this is a limitation of VS Code. The bookmarks interfere with setting breakpoints by clicking in the gutter, so if you want a add breakpoint on a line with a bookmark, right-click on the line to use the pop-up menu. Bookmarks can also be toggled by the pop-up menu created by a right-click.

Also configure the debugger for allocated and structured variables.

* VS Code needs to be configured to display FORTRAN structures or allocated arrays. The method is described in this [GitHub gist](https://gist.github.com/albertziegenhagel/6431811950864bd0009b6a1fa78e7f2b). The method requires launching VS Code from a command window where Intel's `setvars.bat` has initialized the Visual Studio environment. The `.vsix` file can be extracted with 7-Zip or another ZIP utility. Then follow the gist to create the `.vsdbg-config.json` file and the `~/.cppvsdbg/extensions/FEE.link` file. The `.link` file should contain the path to the folder where you placed `.vsdbg-config.json`.


---

### Why this workflow is seamless for the team:

* **Predictable Paths:** By ensuring everyone uses the default Miniconda installation path, `${env:USERPROFILE}/Miniconda3/...` (or Ananconda3/Condaforge3, etc.) will resolve flawlessly on every single team member's Windows machine.
* **No Git Pollution:** Since the configuration is entirely portable, nobody will accidentally modify or commit their own absolute paths (like `C:/Users/your-name/...`) to the shared repository.

The VS Code configuration uses `${env:USERPROFILE}` rather than a developer-specific absolute path. If your Conda distribution or installation location differs from the defaults above, adjust the relevant local settings without committing machine-specific paths.

#### Details on CMake Extensions

What Each Extension Handles
* fortran-lang.linter (Modern Fortran): Drives your syntax highlighting, hover documentation, and automatically communicates with the fortls.exe server tracking the .git exclusions.
* fortls language server analyzes code for relationships such as subroutines, variable types, definitions.
* ms-vscode.cmake-tools (CMake Tools): Provides the status bar UI to easily toggle between your debug, release, and relwithdebinfo configuration presets with one click.
*  ms-vscode.cpptools (C/C++): Supplies the native cppvsdbg debug engine required to read the Visual Studio-style embedded debugging symbols generated by the Intel ifx compiler.

[Return to TOC](#table-of-contents)

## Step 5. Introduction to the VS Code Menus and Icons

This section is for users new to VS Code. The icons in the left toolbar are toggles to open/close special panes.

### Folder View Icon

The folder icon in the top of the left toolbar toggles a explorer folder view. Single-clicking a file opens a preview window and the filename in the tab is displayed in italics. Double click a file name or the preview tab to open the file, and the filename in the tab will be displayed in regular font. The font color provides feedback on the state of a file. When editing a .json, the tab color indicates when all syntax is OK. Also look for squiggly lines which are warnings when editing a .json file.

### Search Icon

The search icon can be used to search for text within projects. See the search topic in the advanced features section to use this effectively.

### Integration with Git

The Git icon provides access to the repo to show modified, untracked and differences. It also provides a visual map of history. When viewing the file listings, the letter at the right of the upper panel shows if items are unmodified, etc. and hovering shows more details. To view a graph of history, in the line labeled Source Control, click the ... and select Views>Graph.

Double-clicking a modified item shows the 'diff' in the file window and provides a popup to 'undo' the changes in hunks.

Files can be staged using the GUI. The terminal accepts git line commands for advanced commands that are not provided in the GUI.

### Run and Debug

The run and debug icon permits running and debugging similar to the Cmake debug icons described below.

### Extensions

The extensions icon provides access to the settings for adding extensions or modifying extension settings.

### Cmake

These features are described below.

### Bookmarks

If you have installed the Bookmark extension, this icon toggles display of the _file_ bookmarks are inserted.

### Terminal Window

At the bottom of VS Code the terminal pane contains tabs for Problems, Output, Debug Console, Terminal subpanes, Ports. The terminal pane display is toggled using Cntl + \` or View>Terminal. (Don't confuse the 'Terminal' menu with 'View>Terminal'). The panes provide the following features:
* Problems - shows issues identified from the language server and linter. The initial display will be incorrect until the project is compiled, so ignore the initial problems until after a project is compiled (built). The language server runs in the background constantly, but the linter updates only when the code is compiled.
* Output - shows VS Code output including the extensions.
* Debug Console - this is the terminal window where the debugger will run. Output can be filtered for specific extensions, such as Modern Fortran to evaluate extension settings or errors. Expressions can be entered in the input line when debugging.
* Terminal - used for command line instructions, including git. Additional terminal subpanes can be opened using the + icon.
* Ports - not used by this project.

[Return to TOC](#table-of-contents)

## Step 6. Understanding Linting and Language Services

When a FORTRAN file is opened, VS Code will run a linter and language service. Here is the breakdown of how **Linters** and **Language Services** differ, how they overlap, and how they actually work together to help you code.

---

### The Short Version

* **A Linter** is like an editor or a proofreader. It looks for stylistic issues, anti-patterns, and bad habits (e.g., "You have an unused variable here," or "Your indentation is inconsistent"). It enforces specific code quality rules and requires that the compiler is run on the project.
* **A Language Service** is like a deeply knowledgeable assistant. It actually *understands* the architecture of your codebase (e.g., "What properties does this object have?", "Where is this function defined?", "What type of data does this variable hold?").

---

### Side-by-Side Comparison

---

| Feature | Linter (e.g., ifx, ESLint, Flake8) | Language Service (e.g., fortls, gopls, Pyright) |
| --- | --- | --- |
| **Primary Goal** | Catch code smells, style violations, and potential bugs. | Provide code intelligence, navigation, and autocomplete. |
| **Scope** | Usually analyzes **one file at a time** line-by-line using abstract syntax trees. | Analyzes the **entire project graph**, understanding how files link together. |
| **Features** | Warning markers, formatting fixes, enforcing naming conventions. | Go-to-definition, find all references, safe renaming (refactoring), IntelliSense/autocomplete. |
| **Customization** | Highly configurable via a local config file (e.g., `.eslintrc`, `.yaml`) where you turn rules on/off. | Usually configured via compiler options, build paths, or workspace settings (like `CMakePresets.json`). |

---

When you edit code file in VS Code, both tools  run in the background simultaneously:

1. **The Language Service** reads your project structure, indexes your modules or files, and figures out how they depend on each other. When you type a function name, it provides the autocomplete dropdown. If you pass an integer to a function expecting a string, it flags a type error.
2. **The Linter** scans the text you just wrote. It notices that you left a trailing whitespace at the end of line 12, or that you used a deprecated syntax that works but isn't best practice. It overlays its own warnings right alongside the language service's errors. The linter problems update only when the file is saved, which triggers the compile.

Many language services can also run linting rules inherently (for instance, a language server might warn you about an unreferenced variable because it already tracked every reference in the project), but dedicated linters exist so teams can customize and enforce strict style rules across an entire engineering team.

The idenfied issues will be displayed in the 'Problems' tab of the terminal pane as described below. **The initial tab may not be up to date when first opening the project.**

[Return to TOC](#table-of-contents)

## Step 7. Browsing Projects, Building Projects

This code distribution depends heavily on .json files to assure a similar development experience for all users. Only a subset of CMake integration settings are configured, and thus some of the available VS CMake controls are not used or addressed in these instructions.

**Open the Cmake panels** by clicking the CMake triangle (with the wrench) in the left toolbar.

This repository uses CMake presets with the Ninja generator for VS Code builds; no standalone CMake app is used in this workflow.

### The CMake Panels

**Note that only a subset of CMake features are used in this project. Thus not all buttons or options will work. Follow the directions below for the utilized features.**

#### Project Status

The Project Status frame is used to select the build configuration and build target.

Set the Build Configuration and Target Project:

* Project Status > Configure: This binds your Configuration/Build Type (e.g., debug, relwithdebinfo, release) directly from your CMakePresets.json.
    * The debug configuration provides unoptimized debugging.
    * The relwithdebinfo configureation is optimized to /O2 which can still be reliably debugged to assure that the optimization did not break anything.
    * The release build configuration is optimized agressively to /O3 which vectorizes DO loops, which makes debugging impractical. However, testing of the release code against benchmarks is still necessary to make sure the final level of optimization did not break anything.

* Project Status > Build: This is where you click the square brackets or current target to change your active target from [all] to a specific asset like PGLdll or PGLEOS. For this project, it is important to build subproject PGLdll separately from the PGLEOS or PGLdllTest because the PGLdll components are configured to link dynamically whereas the PGLEOS is configured to link statically. The CMake system can be easily confused by linking to intermediate files shared by the two methods if they don't have the correct build type within the intermediate files. Thus, don't use [all], and clean all projects when switching from a project that uses dynamic linking to a project that uses static linking. **Once you set the build field in the Project Status frame, the main Build button on the bottom VS Code status bar binds itself to that subproject choice.**

* Project Status > Debug: This is where you click to select the target that will be launched by the Debug button.

#### Project Outline

The top level 'Project Outline' line includes a flyout menu of icons when hovering in the title line.

* Project Outline > Configure Projects: Clears and rewrites the CMake cache. This can be used to 'refresh' the CMake settings.

* Project Outline > Configure all Projects with CMake Debugger: not used due to the potential confusion between the dynamic and static builds for some subprojects.

* Project Outline > Build Projects: Do not use at the overall project level since the dynamic/static build configurations of the subprojects conflict some common build folders are used.

* Project Outline > ... (ellipsis)> Clean All Projects: Useful to purge all cached intermediate build files (.obj, .mod, etc.) and it important to use when switching back and forth between a project with dynamic linking (PGLDll and PGLDllTest) to a project with static linking (PGLEOS).

* Project Outline > ... (ellipsis)> Clean Reconfigure All Projects: Cleans the CMake cache and refreshes across all projects.

##### Repo Listing

Hovering over MDNAProject provides the following links:

* MDNAproject > Configure: rewrites the CMake cache
* MDNAproject > Build: Don't use to avoild clashes with the different link modes in the repo. This is similar to the option in the Project Outline menu described above.
* MDNAproject > Clean: Cleans the previous compile/link for all projects in the repo, similar to the options from the Project Outline menu above.

##### Subproject Listings

* Hovering over a subproject will provide pop-up icons to build or set a bookmark.
* Right-clicking a subproject provides shortcuts to build, open the CMakeLists, or set as the build target. The active build subproject will have a rocket icon.
* Clicking the subproject accordion provides a project browser and code listing specific for the project. Using this menu for the project is the preferred method of browsing code because only the relevant source code is displayed. The subproject CMakelists.txt files are located in these subproject folders.

You may click the build icon beside the subproject target.

At the bottom of the Project Outline frame are the links to parent CMakeLists.txt and CMakePresets.json. These parent files hold the settings for the overall repo. These settings are inherited by the individual projects. VS Code expects you to open the parent repo folder so that CMake reads these parent files and builds the repo CMake cache as soon as the folder is opened.

#### Pinned Commands

The pinned commands panel is not preconfigured for this project.

#### Project Bookmarks

This pane holds **CMake** bookmarks. If the subproject is bookmarked, the file listings will be duplicated in the bookmark pane. (These are not the same as _file_ bookmarks provided by the Bookmarks extension.)

### Bottom toolbar

CMake adds a Build, Debug, Launch icon to the bottom toolbar. If the Build target is not properly set in the Project Status pane, the build icon will build 'all' subprojects, which is undesired. When the subproject is properly set in the Project Status pane, then the bottom toolbar build icon can be used.

The debug icon is used to launch the debugger for the project specified in Project Status Pane.

Once release code is written, it should be placed in the project parent folder and run from that location in the terminal window.

[Return to TOC](#table-of-contents)

## Step 8. Cleaning the projects

This step is not needed on the initial configuration, but it is needed to assure that the caches are not conflicting with previous projects.

**CMake and VS Code cache the project agressively. It is very important to perform a repo level 'clean' when switching between builds of the PGLDll and PGLEos/GammaPA projects.** PGLDll is built for dynamic linking and PGLEos and GammaPA are built for static linking (stand-alone app) so one build type can pollute the other due to the aggressive caching.

Occasionally, you may need to Reload the Window, which completely flushes the CMake settings in memory. To execute a reload, call up the VS Code command pallete using Cntl + Shift + P. Start typing 'Developer: Reload Window' and the dropdown will automatically filter.

**Developer: Reload Window** only restarts the **front-end UI layer** of VS Code. It is essentially like refreshing a tab in your web browser.

Here is exactly what happens behind the scenes when you trigger a reload versus what happens to your build system:

### 8.1. What a Window Reload Actually Does

* Shuts down and restarts the VS Code extension host process.
* Forces extensions like **Modern Fortran** and **CMake Tools** to reboot, re-read your `.vscode/settings.json`, and re-initialize their background language servers (`fortls.exe`).
* Clears out any stuck UI states, extension panics, or "server crashed 5 times" locking codes.

### 8.2. What Happens to CMake Natively

* **The CMake cache is entirely untouched.** When CMake runs its initial configuration step, it writes all its discoveries—like your compiler paths, library configurations, and architecture settings—directly to a physical file on your hard drive (`out/build/debug/CMakeCache.txt`).
* When you reload the window, the newly booted CMake Tools extension simply looks at your project directory, sees that `CMakeCache.txt` already exists, and says, *"Great, the configuration is already done!"* It then loads that exact cached state right back into memory.

---

### The Cheat Sheet: When to use which command

Because they target completely different layers of your development environment, you have to use them strategically depending on what you just modified:

| Action You Just Took | What You Need to Run | Why? |
| --- | --- | --- |
| **Modified `.vscode/settings.json**` or updated a Python/Conda tool path. | `Developer: Reload Window` | Forces VS Code's UI to load your new editor settings. |
| **Modified `CMakePresets.json**` or changed underlying compiler flags/strategies. | `CMake: Delete Cache and Reconfigure` | Forces CMake to delete `CMakeCache.txt` on disk and re-run its compiler sanity checks. |
| **Added a brand new `.f90` source file** to a directory tracked by CMake. | `CMake: Configure` | Tells CMake to scan the directories again and update the build targets without wiping its whole memory. |

If CMake is totally confused, the best option is to close VS Code completely. With VS Code closed, delete the 'out' folder. Nothing important is deleted, this forces CMake to start a fresh configuration when VS Code reloads the folder. When you reopen VS Code, then CMake will reconstruct the 'out' folder.

[Return to TOC](#table-of-contents)

## Step 9. Building a project

**The !DEC$ export for masterDir,PGLinputDir must be toggled ON when compiling PGLDll and toggled OFF when compiling PGLEOS.** This inconvenience will be eliminated in the future.

The selected project build is displayed in the **Project Status>Build** field. That selection is bound to the 'Build' button in the toolbar at the bottom of the CMake Panel.

The build button at the bottom VS Code taskbar will trigger display of the 'OUTPUT' terminal where any compile/link issues are displayed.

[Return to TOC](#table-of-contents)

## Step 10. What is the working directory and where is the compiled code?

The VS Code working directory is set to the repo root so that the code can properly read/write to the 'Input' and 'Output' folders. When code is launched with the 'Debug' button, VS Code will properly use this working directory. When working with source code in this project which specifies paths, typically you will set the directory as '.' within the code. In some files, this is configured by setting RELEASE = .TRUE. . JRE uses additional build folders and has customized some directories for his development work and those directories are not used in the distributed code.

When compiling with the provided CMake files, the compiled code is within the 'out/build/\<build config>/<project folder name>'. Using the VS Code launch button bypasses some settings and tries to run the code from this location, but then the 'Input' and 'Output' cannot be found. To run the project _outside_ of the CMake and VS Code environment, move the executable to the root of the repo.

* When using PGLTest.exe in this way, due to dynamic linking, the PGLDll.dll needs to also be in the repo root. If on a machine without the oneAPI environment, additional runtime .dll files are needed from the oneAPI redist folder. A dependency app can be used to determine which .dll files are needed.

[Return to TOC](#table-of-contents)

## Step 11. Starting the Debugger

After the project is built successfully, use the 'Debug' bug icon in the bottom CMake toolbar to start the debug process.

**Note that if the focus is on the Output terminal, you will need to click the 'DEBUG CONSOLE' tab to see the console.**

When the debugger runs, a pop-up toolbar is displayed. If you close VS Code with the debugger running, or don't properly shutdown the debugger, the breakpoint and debugger session may be in the backgroud when a new debug session is started. Look in the bottom of the CMake panel to see debug sessions and close them if the debug toolbar or the displays appear to be 'messed up'. The bottom status bar will be orange when the debugger is running.

[Return to TOC](#table-of-contents)

## Step 12. Opening Source Code and Setting Breakpoints

Browse in the CMake Project Outline frame to select a source code file to open. Expand the subproject name accordion and the 'References' list to find the code.

**Note that single clicking a file name triggers the 'preview' mode and the file name is displayed in italics. Single clicking a different file switches the preview to the new file but does not 'load' the file.**

Double click the file, or double-click the filename tab to open a file. Then click in the left line numbers to set a breakpoint.

When debugging the 'Run and Debug' window provides the local variables and provides a 'Watch Window' and 'Call Stack'. Clicking the Call Stack will switch to the variables at that point of the stack. For example, to find variables from shared modules, the module variables will be visible lower in the call stack.

[Return to TOC](#table-of-contents)

## 13. A Few Hints on Getting Started with VS Code

### Autocomplete in terminal

By default, pressing TAB during a command typing will suggest autocomplete with the first match. Press escape to cancel the suggestion. If multiple suggestions appear, use the arrow keys and then tab. Or type a few letters and then press Ctrl + Spacebar. VS Code documentation provides instructions to configure automatically.

### Jumping to Subroutine Definitions

When viewing code, right-clicking on a subroutine call or variable will provide a link to the definition, even if in a different file. The back arrow in the top menu bar will then return you to the previous cursor location.

### VS Code Whitespace Settings for Coding

White spaces at the end of lines create confusion in revision control because if the amount of white space changes the line is marked as 'changed'. Also, files should end with a new line. These settings are preconfigured in the .vscode/settings.json.

### Filtering Problems

The Problems window (use Ctrl + \`) is very helpful when viewing code files when linting is enabled. However, warnings can be overwhelming especially when using implicit variable types. On the Problems tab, use the filter box to select which items are shown. Click an item to jump to the location.

Be patient when VS Code first loads. After compiling, initally the problem count may be large, but it may gradually decrease as the linter finds dependent files, etc. as it works in the background. After making a code change, save the file to retrigger the compile.

### Searching within a subproject

CMake generates a file fortls_sources.txt. This is a scratch file - do not commit this file to the repo. Open in a text editor and copy the list of source files from the subproject that you wish to search. Use the magnifying glass in the left toolbar or Ctrl + Shift + F and paste the list into the box for 'files to include'. Then enter the search term in the top box.

### Command Palette

The command palette is opened using Ctrl + Shift + P. Commands are best entered by starting to type the command instead of browsing. Useful commands for this project are:

* CMake: Clean Rebuild
* CMake: Delete Cache and Reconfigure
* Fortran: Restart the Fortran Language Server - if the tool tips are not appearing.
* Developer: Reload Window - a clean reload of the .json files (But does not clean CMake FORTRAN build files)

The recently-used commands appear at the top of the dropdown.

### Viewing Command History

Within the terminal, command history is available.

* Use the up arrow to scroll through the history one command at a time.
* In the command window, use the shortcut keys Ctrl + Alt + R to bring up the same history as Ctrl + Shift + P -> Terminal: Run Recent command.
* History for the current session is available by typing 'history' in the command line.

### Linting Cleanup

The project linting compiles to check syntax and creates '\*\_genmod.f90', '\*\_genmod.smod', and '\*.mod' files in the source code folder. These are hidden in the VS Code view but will be apparent when viewing in File Explorer. These are ignored by Git. A powershell script is provided to clean up the folders. To clean the folders, using windows file explorer, right-click on the powershell script 'clean_genmod.ps1' and select 'Run in Powershell'.

Linting is not expected to work properly on the first load of a project that has no CMake cache. After you set the project and compile, you may need to [reload the window](#step-8-cleaning-the-projects) to eliminate inaccurate  linter issues in the 'Problem' window. The CMake system may evaluate faster than the compiler on first launch of a project without a CMake cache and thus the reload may fix that.

### Working with Customized settings.json

If you have customized settings.json, it should not be committed back to the repository. However, the local changes will prevent merging any updates into your local copy. The following workflow is recommended.

1. View the repo using a Git-aware directory tool that scans for modified files. For example, the VS Code Git folder viewer or GitHub desktop will show if the file is modified, and opening the file will show the diff. Instead of manually copying code diff, the following is easiest.
2. Copy the customized settings.json to a local untracked file, such as setttings_mine.json.
3. Revert settings.json to the distributed version.
4. Execute the merge.
5. Compare settings.json with settings_mine.json to transfer any needed customizations to settings.json.

[Return to TOC](#table-of-contents)

## 14. Adding Files to a Project

The project CMakelists.txt provides a list of files in the order of compiling and linking. To build the file list, a powershell script scan-fortran.ps1 is provided to crawl the folders for module definitions and common blocks. Files with module definitions should be listed before they are used. CMake is 'smart' and recompiles only the set of files needed. If the fourth from the last file is edited, only the last four files will be compiled again before linking. If missing symbols are found when linking, check the file order. When the files are sorted by folder, the order that the files are added is also important.

If source code is placed in new folders, add the new folders to the source code folders in .fortls.

[Return to TOC](#table-of-contents)

## 15. Adjusting Compile Flags

Compile flags are predominantly in CMakePresets.json with some in .vscode/settings.json. To view more warnings, edit to /warn:all. To enforce a specific FORTRAN standard, use /stand:f90 or adjust the Fortran standard, e.g. f18. If you use implicit variable, turn off the declaration warning by editing /warn:declaration to /warn:nodeclaration.

Compile commands for fixed format .f/.for files are within the project CMakeLists.txt. Project-specific flags can be added immediately after the target_compile_options. Flags are appended, and subsequent flags will take precedence in the case of conflicts. The exact compile commands are in out/\<build config>/build.ninja and .ninja.log. Some flags in the logs may use Linux format dashes instead of the Windows format slash. CMake adjusts the flags to the Windows format before passing to the compiler/linker.

[Return to TOC](#table-of-contents)

## 16. Troubleshooting

### Compiler issues
If After running `launch_vscode.bat`, the compiler is not found check that the Intel compiler is on `PATH` in the command window used to run the launch script:
```cmd
where ifx.exe
```
The result should point to `...\Intel\oneAPI\compiler\...\bin\ifx.exe`. If Intel oneAPI is installed outside its default location C:\Program Files (x86)\Intel\oneAPI\, update the `setvars.bat` path in `launch_vscode.bat` locally. Keep this machine-specific change local; do not commit it. Use the procedures in [13. A Few Hints on Getting Started with VS Code](#13-a-few-hints-on-getting-started-with-vs-code) for maintaining customized files.


### Git not working
Git should be in the environment used to launch VS Code. If 'where git' does not find git in the command window used to launch VS Code, then you can either add Git to your path, or you can open the VS Code settings (Cntl + ,). Search in the 'user settings for 'git.path' and insert the path. Close and reopen VS Code.

[Return to TOC](#table-of-contents)
