 [![Licence](https://img.shields.io/hexpm/l/plug.svg?style=flat-square)](http://www.apache.org/licenses/LICENSE-2.0)
[![Build Status](https://img.shields.io/travis/imminfo/ymir.svg?style=flat-square)](https://travis-ci.org/imminfo/ymir)


# Ymir
Probabilistic models of immune receptors assembling

---

### Overview

#### Features

#### Paper

#### Details
Undetailed dependencies. For more information on dependencies see `Installing and using Ymir` subsection.

---

## Installing and using Ymir

You can use Ymir in two ways:

- using pre-made scripts for the most common tasks - computation of generation probabilities (sometimes I call them "assembling probabilities",
estimation of parameters of an immune receptor assembling model from the experimental data or generation of pre-selection immune receptor sequences. 
For this tasks you need to compile some source code, but don't worry! It's very easy, just look at the next subsection.
After compiling, take a look at the `Examples / ready-to-use scripts` section where pre-made scripts are explained in details.

- as an included library in your project. 

### Compiling Ymir

> Note: you can see detailed installation of dependencies in next sections.

Dependencies:

- C++ compiler - either [Clang](http://clang.llvm.org/) or [GCC](https://gcc.gnu.org/). On benchmarks Clang is performing better then GCC.

- [CMake](http://www.cmake.org/download/) - the build system which used in Ymir.

- [Python 2](https://www.python.org/downloads/) - the JsonCPP library which Ymir use need this version of Python.

- [Python 3](https://www.python.org/downloads/) - we need it for some useful scripts like converting input files to the Ymir's format and wrapping calls to Ymir C++ scripts.

- [JsonCPP](https://github.com/open-source-parsers/jsoncpp) - JSON files used for storing metadata about models. 


#### Installation on Ubuntu

Open your Terminal.

Install CMake and Pythons with

    sudo apt-get install cmake
    sudo apt-get install python
    sudo apt-get install python3

Other installation steps (JsonCPP and compiling) are equal to those for Mac users.

#### Installation on Mac OS

Open your Terminal.

Install [Homebrew](http://brew.sh/) by typing this into the Terminal:

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

After that you need to install CMake and Pythons with the following commands:

    brew install cmake
    brew install python
    brew install python3

At this point you can either install the JsonCPP library and download benchmark files by yourself (see next paragraphs) or 
simply go to [Ymir's GitHub](https://github.com/imminfo/ymir/releases), download the latest release source code called `Source_code_FULL.zip` and unzip it.
You will see the new folder called `ymir`. In both ways to compile Ymir you need to make your Terminal point to Ymir's folder.

To install JsonCPP you should clone the repository and install JsonCPP library by the following commands:

    git clone https://github.com/imminfo/ymir.git
    cd ymir
    curl -sL https://github.com/open-source-parsers/jsoncpp/archive/master.zip > jsoncpp.zip
    unzip jsoncpp.zip
    rm jsoncpp.zip
    cd jsoncpp-master
    python amalgamate.py
    mv ./dist/jsoncpp.cpp ../Ymir/src/jsoncpp.cpp
    mv ./dist/json ../include/json
    cd ..

And now, finally when you are in Ymir's folder, you can compile Ymir:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    cmake --build . && cd ..

It will build a number of Ymir programs to run: tests, benchmarks and ready-to-use scripts.

You can test Ymir with running

    python3 ytest.py

For benchmarks you need to download benchmark files (if you haven't download `Source_code_FULL.zip`):

    curl -sL https://github.com/imminfo/ymir/releases/download/v1.0-pre3/benchmark-files.zip > benchmark-files.zip 
    unzip benchmark-files.zip
    rm benchmark-files.zip 
    mv benchmark-files/ ./benchmark/data/
    rm -r benchmark-files

To run benchmarks execute

    python3 ybenchmark.py

To see how you can use pre-made scripts go to the `Examples / ready-to-use scripts` section.

### Using Ymir as a library

Ymir is a header only library, however, you need to put the `jsoncpp.cpp` file to the `Ymir/src` folder.

In your CMakeLists.txt files you should put `include_directories(include)` if Ymir's folder is in the `include` directory.

Next text is assumed that your include path is setted to `$YMIR_FOLDER/`.

For using only the core data structures MAAG and MAAGBuilder add:

#include <Ymir/Graph>

For computing and generation assembling probabilities from files:

#include <Ymir/Model>

For the statistical inference of marginal parameters of generation models:

#include <Ymir/Inference>

For super secret unimplemented things:

#include <Ymir/Pattern>

Wow-wow, easy here. They are not implemented yet.

---

## Main helper script

python3 pyymir.py

* `python3 pyymir.py -v` - print the current version of Ymir.
* `python3 pyymir.py -a` - list all available algorithms for the statistical inference of marginal probabilities and their aliases.
* `python3 pyymir.py -m` - list all available models and their aliases.
* `python3 pyymir.py -f` - list all available converters for input formats and their aliases.
* `python3 pyymir.py -s` - list all available pre-made scripts.

---

## Examples / ready-to-use scripts

After compiling an executable files will appear in the `build` folder.

You can download sample files with outputs from [tcR](https://imminfo.github.io/tcr/) package and [MiTCR](http://mitcr.milaboratory.com/) software with

curl -sL https://github.com/imminfo/ymir/releases/download/v1.0-pre3/example-files.zip > examples.zip 
unzip examples.zip 

### Compute generation probabilities of human TCR-alpha data from tcR output files

The most common task with generation probabilities is to compute them. 

python3 compute.py -i <INPUT> -f <FORMAT> -m <MODEL> [-o <OUTPUT>] [-p] [-l]

> Note: run `python3 compute.py -h` to see a help message.

* `<INPUT>` - input file (text or gzipped) or a folder with input files (of the same format) or a list of space-separated files and/or folders in any combinations.
* `<FORMAT>` - format of input files (tcR, MiTCR, MiGEC, etc.) as an alias or as a Python 3 class from your module linked to the '$YMIR_HOME/converters' package. For a list of possible input formats and their aliases in this Ymir distribution run `$python3 pyymir.py -f`.
* `<MODEL>` - either an alias of the one from available models in  Ymir or a path to a folder with a model's .json file. For a list of available models in this Ymir distribution run `$python3 pyymir.py -m`.
* `<OUTPUT>` - optional path to the output folder for output files (default is `./ymir_genprob/`).
* `-p` - recompute or use predefined model's gene usage (default is to recompute, to change it add `-p` to your script call).
* `-l` - add this to leave converted files (default is to remove converted files).

Example command to compute generation probabilities from tcR file:

python3 compute.py -i ./examples/tcr.htra.txt -f tcr -m htra


### Estimate human TCR-beta generation model parameters using EM-algorithm from MiTCR output files

python3 inference.py -i ??? -f mitcr -m htrb

> Note: run `python3 inference.py -h` to see a help message.

Example command to perform statistical inference from MiTCR file:

python3 inference.py -i ./examples/mitcr.htrb.txt -f mitcr -m htrb

Run `python3 pyymir.py -m` to view all available models and their aliases.

### Generate artificial human TCR-alpha repertoire before selection

python3 generate.py -c 500000 -m htra

> Note: run `python3 generate.py -h` to see a help message.

asdasd

---

## Types of input files and their formats

### Model files

#### Main JSON file

#### Gene segments files

#### Event probability files

### Cloneset files

#### Available parsers

> Reminder: run `python3 pyymir.py -f` to view all available formats and their aliases.

---

## Benchmarks

---

For full documentation visit [mkdocs.org](http://mkdocs.org).

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs help` - Print this help message.

mkdocs.yml    # The configuration file.
docs/
index.md  # The documentation homepage.
...       # Other markdown pages, images and other files.
