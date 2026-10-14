# Contributing to SOFICS

SOFICS is a research code developed and maintained by a small academic group. We welcome bug reports, questions, and contributions, and we are glad to hear from anyone using the software. Please keep in mind that development happens alongside our research, so responses may take a few days.

## Reporting a bug

A report is most useful when it includes:

* the SOFICS version or git tag, and the commit hashes of the two submodules (`git submodule status`),
* the compiler, MPI implementation, and CMake versions (`gcc --version`, `mpiexec --version`, `cmake --version`),
* which example you were running, or the input files if it is your own case,
* the command you ran and the error message, with the relevant part of `log.out`.

## Contributing code

Fork the repository, make your changes on a branch, and open a pull request against `main`. Please keep a pull request focused on one issue, and describe how you tested the change. If your change is in the fluid or structural solver rather than in SOFICS itself, it belongs in that solver's repository, since SOFICS only records which commit of each solver it builds against.

## Contact

See the [Contact and Support](README.md#contact-and-support) section of the README for whom to approach about each component.
