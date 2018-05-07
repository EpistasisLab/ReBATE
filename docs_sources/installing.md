#Installing ReBATE

## Prerequisites
All of the necessary Python packages can be installed via the [Anaconda Python distribution](https://www.continuum.io/downloads), which we strongly recommend that you use. We also strongly recommend that you use Python 3 over Python 2 if you're given the choice.

ReBATE requires that the following external Python packages be installed (all included in Anaconda):

argparse, time, sys, os, IO, Common, numpy, math, pandas, scipy.spatial.distance, operator, csv, distutils.core, distutils.extension, Cython.Distutils, datetime

NumPy and SciPy can be installed in Anaconda via the command:

```
conda install numpy scipy
```
### Additional Steps For Windows Users
Here we discuss the additional prerequisites for running ReBATE on a Windows operating system. First you will need a command line terminal for Windows (e.g. [Cygwin](https://www.cygwin.com/install.html), or [GitBash](https://gitforwindows.org/)) with Anaconda properly installed. From this terminal you can compile the Cython code, and run ReBATE. Second you will need a C compiler for compiling the Cython code once before being able to run ReBATE on your Windows machine. 

There are a number of possible ways to get the C compiler working on your Windows terminal, all of which will depend on the version of windows, the version of python/Anaconda, and whether it is 32-bit or 64-bit. Be aware that there may be some troubleshooting in getting the C compiler operational in your command line terminal. Below we outline steps that worked for us in the spring of 2018 using Windows 10, and Python 3.5.2 with Anaconda 4.0 (64-bit), using GitBash as our terminal.

1.) Ensure that setuptools is updated. Run the following in your terminal: 
```
pip install -upgrade setuptools
```
2.) [Download/install Visual Studio Community 2017](https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2017). This includes the necessary C compiler. It is not necessary to install all of the visual studio components, as this can be slow and take up a good deal of space. However make sure that you download and install the 'Desktop Development with C++' and the 'Python Development' workloads as detailed at the following [link](https://docs.microsoft.com/en-us/visualstudio/python/working-with-c-cpp-python-in-visual-studio).  Within the Python Development workload, also select the box on the right for Python native development tools. 

3.) Once Visual Studio has been successfully installed, there was a remaining bug that required some troubleshooting for the C compiler to be found by the terminal. Try compiling Cython as described below and if it doesn't work try something like the following fix that worked for us: 

(A) [Add this to your PATH environment variables](https://www.java.com/en/download/help/path.xml): C:\Program Files (x86)\Windows Kits\10\bin\x64   

(B) Copy these two files (rc.exe and rcdll.dll) from (C:\Program Files (x86)\Windows Kits\8.1\bin\x86) to (C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin)

4.) At this point you should be able to compile Cython as described below.

## Compile Cython
Once these prerequisites are installed, it will be necessary to compile the Cython code on the respective operating system within which ReBATE will be run. It is only necessary to do this once, not every time ReBATE is run. This happens in two stages (1) a .pyx file is compiled by cython to a .c file, then (2) the .c file is compiled by a C compiler to a .so file (or a .pyd file for Windows). 

Simply run the following file included with ReBATE to produce the .c and (.so or .pyd) files:  
```
./make.sh
```
If there is need to recompile the Cython files, first remove the previous .c and (.so or .pyd) files by running: 
./clean.sh
```
```