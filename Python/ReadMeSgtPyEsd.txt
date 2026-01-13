Notes on SgtPy
0. Be careful to find the dir of the python installation that you will be using for SgtPy and do all your pip installs from a cmd window in that dir. There may be multiple pythons installed on your computer. It is best to reduce that number to 1-2 so you don't get confused.
1. You must pip install Cython before attempting to install SgtPy. This will require install of C++ if you have not already. There are other installs like numpy, scipy, listed on GitHub, but these are usually installed already when you use anaconda.
2. SgtPyEsd also requires pip install ChemPy, keyboard.
3. Modify PGLInputDir in GlobConst.py to reflect the dir where you have installed PGLWrapper.
4. You may encounter "deprecation" warnings during install of SGTPy. These are inherent to SGTPy. If ignoring them leads to 
a serious failure, you will need to contact Gustavo Chaparro. 