# MolcasTo47

This repository provides a simple Fortran code to extract data from 
calculations with Molcas/OpenMolcas (https://gitlab.com/Molcas/OpenMolcas)
and create a generic `FILE47` input for the NBO program 
(https://nbo.chem.wisc.edu/)

Copyright 2024 Aleksandr Zaichenko and Jochen Autschbach

## Installation

Adjust settings in `Makefile` for your compiler/operating system
and type `make` in the top-level directory.

Requires BLAS and HDF5 libraries. On an Ubuntu-based Linux system, you can
install them with 
```
 apt install libhdf5-dev libblas-dev
```

## Usage

Molcas needs to be run such that the basis set information and the
overlap matrix are stored in an `h5` file, which will be processed by
`molcasto47`. An optional argument can be given for an orbital file in
Molcas format. If the latter is provided, the density matrix to be used
in the NBO run will be crated from the orbital file data.

```
 molcasto47 H5FILE [orbitalfile] [-fock] 
```

Output is `File.47`

Optional third flag `-fock` for writting of the Fock matrix.  
