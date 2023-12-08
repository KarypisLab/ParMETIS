# ParMETIS

[![ci](https://github.com/scivision/ParMETIS/actions/workflows/ci.yml/badge.svg)](https://github.com/scivision/ParMETIS/actions/workflows/ci.yml)

ParMETIS is an MPI-based library for partitioning graphs, partitioning finite element meshes,
and producing fill reducing orderings for sparse matrices. The algorithms implemented in
ParMETIS are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint
partitioning schemes developed in our lab.

## Downloading ParMETIS

You can download ParMETIS by simply cloning it using the command:
```
git clone https://github.com/KarypisLab/ParMETIS.git
```


To build ParMETIS you can follow the instructions below:

```sh
cmake --workflow --preset default
```



### Definitions of supported data types

ParMETIS uses the same data types for integers and floating point numbers (32/64 bit
integers and single/double precision floating point numbers) as used when configuring
and building METIS.

## Copyright & License Notice

Copyright 1998-2020, Regents of the University of Minnesota
