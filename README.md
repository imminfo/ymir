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

- [Python 3](https://www.python.org/downloads/) - we need it for some useful scripts like converting input files to
the Ymir's format and wrapping calls to Ymir C++ scripts.

- [JsonCPP](https://github.com/open-source-parsers/jsoncpp) - JSON files used for storing metadata about models. 

You can load Ymir to [CLion](https://www.jetbrains.com/clion/) and compile with it.

#### Installation on Ubuntu

#### Installation on Mac OS

Open your Terminal.

Install [Homebrew](http://brew.sh/) by typing this into the Terminal:

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

After that you need to install CMake and Pythons with the following commands:

    brew install cmake
    brew install python
    brew install python3

Go to [Ymir's GitHub](https://github.com/imminfo/ymir/releases) and download the latest release source code called `Source code + JSON` and unzip it.
You will see the new folder called `ymir`. Now you need to go via Terminal to Ymir's folder with `cd <path to ymir's folder here>/ymir`, e.g.
`cd ./ymir`.

Alternatively you could clone the repository and install JsonCPP library by the following commands:

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
    cmake --build .

Now it will build a number of programs to run Ymir: tests, benchmarks and pre-made scripts.

You can test Ymir with:
    
    ./test/Test

And run benchmarks with:

    curl -sL https://github.com/imminfo/ymir/releases/download/v1.0-pre2/benchmark-files.zip > benchmark-files.zip 
    unzip benchmark-files.zip
    rm benchmark-files.zip 
    mv benchmark-files/ ./benchmark/
    rm -r benchmark-files
    ./benchmark/Benchmark    

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

* `python3 pyymir.py -v` - print the version of Ymir.
* `python3 pyymir.py -a` - list all available algorithms for the statistical inference of marginal probabilities and their aliases.
* `python3 pyymir.py -m` - list all available models and their aliases.
* `python3 pyymir.py -f` - list all available converters for input formats and their aliases.
* `python3 pyymir.py -s` - list all available pre-made scripts.

---

## Examples / ready-to-use scripts

After compiling an executable files will appear in the `build` folder.

You can download sample files with outputs from [tcR](https://imminfo.github.io/tcr/) package and [MiTCR](http://mitcr.milaboratory.com/) software from here [here]().

### Compute generation probabilities of human TCR-alpha data from tcR output files

    python3 compute.py -i ??? -f tcr -m htra
    
sdfsdf

### Estimate human TCR-beta generation model parameters using EM-algorithm from MiTCR output files

    python3 inference.py -i ??? -f mitcr -m htrb
    
Run `python3 pyymir.py -f` to view all available models and their aliases.

### Generate artificial human TCR-alpha repertoire before selection

    python3 generate.py -c 500000 -m htra

asdasd

---

## Types of input files and their formats

### Model files

### Gene segments files

### Cloneset files

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
