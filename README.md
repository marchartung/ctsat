# ctsat
Template based SAT solver. Basic heuristics regarding branching, reduce and restart are runtime dynamic but compiled into different solver versions. Supports parallel execution using PThread and MPI.

The sequential version can be build in the directory core using the given make file, the Pthread parallel version in the parallel and the MPI version in the mpi directory. Please notice that in order to build the mpi version, a mpicxx compiler must be present.

Build targets: d ... debug, s ... standard with assertions, r ... release without assertions and highest optimization

Currently builds only with GCC C++ compiler supporting C++11 on linux.
