## Using `gfortran` on Windows

`gfortran` is part of the GNU compiler collection and runs in Unix-based operating systems, i.e. Linux or Mac OS. Windows users have two options to build maespa with `gfortran`:

* install WSL (Windows Sub-system for Linux) and use the Linux tools to build a Linux-compatible executable:
    1. follow [this tutorial](https://docs.microsoft.com/en-us/windows/wsl/install-win10) to install WSL and a Linux distribution (such as Ubuntu);  
    1. install `gfortran` and `make` from your linux distribution (e.g. open Ubuntu, or open a terminal and type `bash`):
    ```bash
    sudo apt install gfortran make
    ```
    1. open a terminal at the root of the maespa project, and type:
    ```bash
    make all
    ```
    `make` follows the instructions from the `Makefile` and builds the project using the Linux version of `gfortran`. The maespa executable that is produced works only for Linux, so you can only use it from your Linux sub-system.

* or/and install MinGW (or Cygwin). There are two main versions of MinGW, but we will use MinGW for Win64 as advised from the [GCC wiki](https://gcc.gnu.org/wiki/GFortranBinaries). Follow these steps to install it:
    1. Download the installer from [here](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe/download), and install it.  
    1. Add it to your path: press the `Windows key`, type `environment`, press `enter`, click on `environment variables`, click on `path,` and `New`, add the complete path of installation until the bin folder, e.g. if you installed `mingw-w64` in `Program Files`:  `C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin`. The bin folder may be at a different path on your computer.  
    1. Rename `mingw32-make.exe` from the bin folder into `make.exe`, and also `mingw32-gfortran.exe` into `gfortran.exe`.  
    Then to build maespa:
    1. open a command or a Powershell at the root of the maespa project, and type:
    ```shell
    make all
    ```
    `make` follows the instructions from the `Makefile` and builds the project using the MinGW version of `gfortran`. The maespa executable that is produced works directly for your windows operating system.
