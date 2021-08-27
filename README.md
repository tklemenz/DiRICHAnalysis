# Introduction
Analysis framework for the DiRICH based read-out to be used with the RadMap modules as a trigger for the Photon-THGEM detector.
Repository is initially forked from the RadMapAnalysis by tklemenz, which is designed to work with a Padiwa based read-out system and will be modified.


# Install guide

## Prerequisites

* root
I recommend installing root from source as described in https://root.cern/install/#build-from-source
* cmake
* fmt
* boost

### Recommended

* doxygen
* graphviz

## Installation

1. Go to your favorite destination
2. `git clone https://github.com/bulukutlu/DiRICHAnalysis`
3. install fmt (https://fmt.dev/latest/usage.html up to installing the library)
4. install boost (https://www.boost.org/doc/libs/1_55_0/doc/html/bbv2/installation.html)
5. source your `thisroot.sh`
6. `source RadMapAnalysis/compile.sh`
7. It's that easy

# How to use it

## Framework

When you want to open a root file containing objects from the framework you need to load the RadMapAnalysis.so (.dylib if you're a mac disciple) library.

The easiest way I found for now is to create an alias in your `~/.bashrc` where the library is loaded when you call `root`.
It could look something like this:

`alias radROOT='root -l -e "gSystem->Load(\"/path/to/DiRICHAnalysis/build/libRadMapAna.so\");"'`

Then you can simply call `radROOT` and the objects from the framework are known to root.

## Doxygen

The doxygen documentation has just started and is not complete by any means!

Go to the root folder of the project (should be `/path/to/DiRICHAnalysis`) and run `doxygen doxygen.conf`. This needs to be done whenever you made changes to the documentation and want to update the doxygen files.

Then open your browser (firefox works well in this case) and put `/path/to/DiRICHAnalysis/doxygen/html/index.html` in the URL line. One of the most useful things is the TODO list. Things can be added by putting `/// @todo This needs to be done.` in the corresponding file. But appart from that it is also pretty useful for general documentation.



# General stuff

* compile.sh should not be changed!
* Executables will be located in the bin folder.
* As the only style guide I would suggest to NOT use tab but only spaces. Otherwise the code looks horrible in different editors...

Please only commit the following files or files in the following folders if changes were made:

	- CMakeLists.txt
	- include
	- src
	- macros
	- README.md
	- doxygen.conf

Do NOT commit:

	- bin
	- build
	- doxygen


Feel free to not only add macros but also add missing functionality to the classes or create new classes that you think are necessary.

If you don't know how to add classes and/or macros just go along the lines of the Dummy (add to LinkDef and CMakeLists.txt etc...).

Thanks for contributing!!
