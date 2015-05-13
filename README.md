Planet\_INC
===========

Planet\_INC is a planetary code for atmospheric calculations using [Antioch](https://github.com/libantioch/antioch)
for the kinetics chemical computations and [GRINS](https://github.com/grinsfem/grins) for the
solving part.

Dependencies
============

Planet\_INC requires [Antioch](https://github.com/libantioch/antioch), [GRINS](https://github.com/grinsfem/grins),
[Eigen](http://eigen.tuxfamily.org/), GSL, and GRINS dependencies.

The required dependencies are, in the order of building:
  - [Boost](http://www.boost.org/)
  - [PETSc](http://www.mcs.anl.gov/petsc/). Advised configation options are
```
        ./configure --with-clanguage=C++ --with-shared-libraries \
        --with-mpi-dir=$MPI_DIR \
        --with-mumps=true --download-mumps=1 \
        --with-metis=true --download-metis=1 \
        --with-parmetis=true --download-parmetis=1 \
        --with-blacs=true --download-blacs=1 \
        --with-scalapack=true --download-scalapack=1 \
        --with-superlu-dist=true --download-superlu-dist=1 \
        --with-x11=0 \
        --with-x=0 \
        --with-blas-lapack=true --download-blas-lapack \
        --with-errorchecking=1 \
        --with-debugging=1 \
        --with-clib-autodetect=0 \
        --prefix=/the/right/path
```
  - [SLEPc](http://slepc.upv.es/). 
  - [libmesh](https://github.com/libMesh/libmesh). Configure with the following options:
```
        ../libmesh/configure --prefix=/where/to/install/libmesh --enable-everything --with-metis=PETSc
```

GSL, Antioch and Eigen can be built independently. Note that Antioch is
also a dependency of GRINS.

Installation and running
========================

Download the code in a convenient place:
```
$ cd /convenient/place/
$ mkdir src
$ cd src
$ git clone https://Planet-INC/Planet-INC .
```

Then compile and install the code in another folder:
```
$ mkdir ../built-dir
$ cd ../built-dir
$ ../src/configure --prefix=/where/the/executable/should/be
$ make
$ make test (optional)
$ make install (might require root provilege)
```

To run the executable, the advised command is:
```$ mpiexec -p NUMBER_OF_CORES planet INPUT_FILE -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package superlu-dist```

The options -ksp\_type, -pc\_type, and -pc\_factor\_mat\_solver\_package are PETSc options. See PETSc documentation for
details and other desired options.
