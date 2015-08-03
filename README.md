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

- as a library in your project or 

- using pre-made scripts for the most common tasks - computation of generation probabilities,
estimation of parameters of an assembling model from the experimental data or
generation of pre-selected immune receptor sequences. For this tasks you need to compile
some source code, but don't worry! It's very easy, just look at the next subsection. After compiling,
take a look at the `Examples / ready-to-use scripts` section where pre-made scripts are explained in details.

### Compiling Ymir

Dependencies:

- [Python 3](https://www.python.org/downloads/)

- [JsonCPP](https://github.com/open-source-parsers/jsoncpp)

- [MPFR (non-necessary due to license)]()

assembling / generation

Targets: scripts, tests, benchmarking, lib

### Using Ymir as a library

Ymir is a header only library

For using only the core data structure MAAG and MAAGBuilder use:

    #include <Ymir/Graph>

For computing and generation from files:

    #include <Ymir/Model>

For statistical inference of parameters:

    #include <Ymir/Inference>

For super secret unimplemented things:
    
    #include <Ymir/Pattern>

wow-wow, easy here. There are not implemented yet.

---

## Examples / ready-to-use scripts

After compiling an executable files will appear in the `build` folder.

### Compute generation probabilities of human TCR-alpha data

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs help` - Print this help message.

### Estimate human TCR-alpha generation model parameters using EM-algorithm

### Generate artificial human TCR-alpha repertoire before selection

---

## Input file formats

### Model files

### Gene segments files

### Cloneset files

---

## Testing

---

## Benchmarking

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
