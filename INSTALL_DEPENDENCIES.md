# Installing External Dependencies

This guide covers the installation of external libraries required to build **Aero-S** (nonlinear dynamics) and **M2C** (compressible fluid dynamics) as part of SOFICS. The two solvers communicate through MPI, so a parallel build is assumed throughout.

Installation of Dakota, Gmsh, PETSc, and Boost is not handled by this guide; refer to their official instructions: [Dakota](https://dakota.sandia.gov/), [Gmsh](https://gmsh.info/), [PETSc](https://petsc.org/release/install/), and [Boost](https://www.boost.org/doc/user-guide/getting-started.html).

**PETSc and Boost are required by M2C.** Both are widely used in the computational physics community and are available as modules on most HPC clusters, so check what your system already provides (e.g. `module avail petsc`, `module avail boost`) before building either from source. On your local machine, you will likely need to build them yourself by following their respective installation guides.

M2C aborts configuration unless `PETSC_DIR` is exported:

```sh
export PETSC_DIR=/path/to/petsc
```

**Note:** These instructions target **Ubuntu/Debian** systems. Package names and commands may differ on other distributions.

---

## Prerequisites

|  Requirement |               Minimum Version                |
|--------------|----------------------------------------------|
| GCC / G++    | 7+ (C++17 support required by M2C)           |
| gfortran-10  | same version as GCC (required by Aero-S)     |
| CMake        | 3.31+ (required by Aero-S)                   |
| Make         | GNU Make                                     |
| flex / bison | any recent version (required by M2C parser)  |
| PETSc        | no minimum enforced; `PETSC_DIR` must be set |
| Boost        | 1.72+ (required by M2C)                      |

```sh
sudo apt-get update
sudo apt-get install build-essential gfortran-10 cmake flex bison
```

---

## Required Libraries

### 1. MPI (REQUIRED)

```sh
sudo apt-get install libopenmpi-dev openmpi-bin
```

Verify:

```sh
which mpicc mpicxx
mpiexec --version
```

### 2. Eigen3 (>= 3.3) (REQUIRED)

Header-only linear algebra library used by both Aero-S and M2C.

```sh
sudo apt-get install libeigen3-dev
```

If installed to a non-standard location:

```sh
cmake -B build -DEIGEN3_INCLUDE_DIR=/path/to/eigen3 .
```

### 3. BLAS and LAPACK (REQUIRED)

Numerical linear algebra routines used by Aero-S.

```sh
sudo apt-get install libopenblas-dev liblapack-dev
```

---

## Recommended Libraries

At least one sparse direct solver is recommended for Aero-S. Without one, the code falls back to a built-in solver.

### 4. MUMPS (Sparse Direct Solver)

The parallel variant is used since MPI is enabled. Parallel MUMPS also requires ScaLAPACK and BLACS (see below).

```sh
sudo apt-get install libmumps-dev
```

If installed to a non-standard location:

```sh
cmake -B build -DMUMPS_INCLUDE_PATH=/path/to/MUMPS/include .
```

### 5. ScaLAPACK and BLACS

Required for parallel MUMPS. Without these, MUMPS will not be enabled even if its headers and libraries are found.

```sh
sudo apt-get install libscalapack-openmpi-dev libblacs-openmpi-dev
```

### 6. SPOOLES (Alternative Sparse Direct Solver)

An alternative to MUMPS. Aero-S can use either or both.

```sh
sudo apt-get install libspooles-dev
```

If installed to a non-standard location:

```sh
cmake -B build \
  -DSPOOLES_INCLUDE_PATH=/path/to/spooles \
  -DSPOOLES_spooles_LIBRARY=/path/to/spooles.a \
  -DSPOOLES_spoolesMT_LIBRARY=/path/to/spoolesMT.a \
  .
```

### 7. METIS (Mesh Partitioning)

Useful for domain decomposition in parallel runs.

```sh
sudo apt-get install libmetis-dev
```

---

## Build

### Quick setup (all required + recommended packages):

```sh
sudo apt-get install build-essential gfortran cmake flex bison \
  libopenmpi-dev openmpi-bin \
  libeigen3-dev libopenblas-dev liblapack-dev \
  libmumps-dev libscalapack-openmpi-dev libblacs-openmpi-dev
```

### Configure and compile:

```sh
cd sofics
cmake -B build -DCMAKE_INSTALL_PREFIX="$PWD/build/sofics" .
cd build
make -j$(nproc)
```

This builds both Aero-S (`build/bin/aeros`) and M2C (`build/bin/m2c`).

Set `CMAKE_INSTALL_PREFIX` at configure time, omitting it defaults to `/usr/local`.

### Install (required before running the examples):

```sh
cd ..
cmake --install build
```

This creates `build/sofics/bin`, including `driver.sh`, which the `dakota.in`
files in `examples/` reference via `analysis_driver`. `./build.sh` (see
`README.md`) does this automatically.

### Verify:

```sh
./sofics/bin/aeros
./sofics/bin/m2c
```

### Check the Aero-S build summary:

At the end of the CMake configuration step, a summary table is printed:

```
=================================================
           Summary of build options
-------------------------------------------------
MPI                       YES
Distributed FETI:         YES
Aeroelastic:              YES
Mumps:                    YES
Scalapack:                YES
Arpack:                   NO
Spooles:                  NO     (or YES if installed)
Cholmod:                  NO
Pardiso:                  NO
Acme:                     YES
Metis:                    NO     (or YES if installed)
Eigen template library:   YES
OpenMP:                   YES
Build type:               Release
=================================================
```
