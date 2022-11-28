# ParMETIS 

ParMETIS is an MPI-based library for partitioning graphs, partitioning finite element meshes, 
and producing fill reducing orderings for sparse matrices. The algorithms implemented in 
ParMETIS are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint 
partitioning schemes developed in our lab.

##  Downloading ParMETIS

You can download ParMETIS by simply cloning it using the command:
```
git clone https://github.com/KarypisLab/ParMETIS.git
```

## Building the ParMETIS library

To build ParMETIS you can follow the instructions below:

### Dependencies

General dependencies for building ParMETIS are: gcc, cmake, build-essential, and an MPI library. 
In Ubuntu systems these can be obtained from the apt package manager (e.g., apt-get install cmake, mpich, etc) 

```
sudo apt-get install build-essential
sudo apt-get install cmake
```

In addition, you need to download and install
[GKlib](https://github.com/KarypisLab/GKlib) and 
[METIS](https://github.com/KarypisLab/METIS) by following the instructions there. 


### Building and installing ParMETIS  

ParMETIS is primarily configured by passing options to make config. For example:

```
make config cc=mpicc prefix=~/local
make install
```

will configure ParMETIS to be built using mpicc and then install the binaries, header files, and libraries at 

```
~/local/bin
~/local/include
~/local/lib
```

directories, respectively.

### Common configuration options are:

    cc=[compiler]     - The C compiler to use [default is determined by CMake]
    shared=1          - Build a shared library instead of a static one [off by default]
    prefix=[PATH]     - Set the installation prefix [~/local by default]
    gklib_path=[PATH] - Set the prefix path where GKlib has been installed. You can skip
                        this if GKlib's installation prefix is the same as that of ParMETIS.
    metis_path=[PATH] - Set the prefix path where METIS has been installed. You can skip
                        this if METIS' installation prefix is the same as that of ParMETIS.

### Advanced debugging related options:

    gdb=1           - Build with support for GDB [off by default]
    debug=1         - Enable debugging support [off by default]
    assert=1        - Enable asserts [off by default]
    assert2=1       - Enable very expensive asserts [off by default]

### Other make commands

    make uninstall
         Removes all files installed by 'make install'.

    make clean
         Removes all object files but retains the configuration options.

    make distclean
         Performs clean and completely removes the build directory.


### Definitions of supported data types

ParMETIS uses the same data types for integers and floating point numbers (32/64 bit
integers and single/double precision floating point numbers) as used when configuring
and building METIS.


## Copyright & License Notice
Copyright 1998-2020, Regents of the University of Minnesota


