- [Report new Bug Issue](https://github.com/RSE-Cambridge/qe-gpu/issues/new)
- [Guidelines for Contributing](CONTRIBUTING.md)
- [Pull Request Template](.github/PULL_REQUEST_TEMPLATE.md)
- [Project License](License)


## GPU-accelerated Quantum ESPRESSO (QE-GPU)

This is an open-source custom version of Quantum ESPRESSO with embedded GPU
support based on CUDA FORTRAN. This product has been made possible thanks to
the effort of the [NVIDIA](http://www.nvidia.com/page/home.html) HPC Software
and Benchmarks Group. This version is maintained by
[Filippo Spiga](https://github.com/fspiga), contributions are encouraged. Partial 
support was provided by [E4 Computer Engineering SpA](https://www.e4company.com/en/)
via the European PRACE Pre-Commercial Procurement project (Phase 3). To contribute
please refer to the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md)


### Requirements

The [PGI](http://www.pgroup.com/products) compiler version 17.4 or above is required to use QE-GPU. 
We suggest the latest community edition 17.10 (freely available from [PGI](http://www.pgroup.com/products/community.htm).
It containes CUDA SDK and pre-built Open MPI for parallel execution (check the
[PGI Installation Guide](http://www.pgroup.com/doc/pgiinstall174.pdf) how to 
install it). **No other compilers are supported**

You need NVIDIA TESLA Kepler (K20, K40, K80) or Pascal (P100) or Volta (V100). 
No other cards are supported. NVIDIA TESLA P100 and V100 are strongly recommend 
for their on-board memory capacity and douple precision performance.

This version of QE-GPU it is based on Quantum ESPRESSO v6.1. It runs **exclusively** 
in parallel. For x86 architecture, you also also recent version Intel Math Kernel 
Library (MKL). For POWER architecture, make sure you can access IBM Engineering and 
Scientific Subroutine Library (ESSL) for maximum performance.


### Installation

To compile QE-GPU there is no automatic procedure. You must copy a `make.inc` template from "install/" directory into the main directory, edit it based and then run make. Make sure `GPU_ARCH` and `CUDA_RUNTIME` are specified correctly and various paths to libraries and header files point into the correct locations.

By invoking _make_ alone a list of acceptable targets will be displayed. Binaries go in "bin/". Read comments in the `make.inc` templates to customize it further based on your ebvironment and where math libraries are located. The architectures/environments supported are x86-64, POWER and CRAY.

These templates are available:
* `make.inc_x86-64` to compile on a generic x86-64 machine with NVIDIA GPU
* `make.inc_CRAY_PizDaint` to compile on Piz Daint at CSCS, CRAY XC30 with P100 GPU (`GPU_ARCH=60`)
* `make.inc_POWER_DAVIDE*` to compile on PRACE "DAVIDE" machine at CINECA, based on POWER8 with GPU (`GPU_ARCH=60`)
* `make.inc_POWER_SUMMITDEV` to compile on early access system SUMMITDEV at ORNL, based on POWER8 with GPU (`GPU_ARCH=60`)

The QE-GPU package has been reduced in size to the minimum essential. For more
information, please refer to the general documentation provided with the full
Quantum ESPRESSO suite or visit the official web site
[http://www.quantum-espresso.org/](http://www.quantum-espresso.org/)


### Citation

If you use the code for science or any form of scientific and technical dissemination activity, we kindly ask to cite the code using the two following references:
* Romero, J., Phillips, E. Fatica, M., Spiga, F.: GPU-accelerated Quantum ESPRESSO, Version 1.0 (November 2017), http://doi.org/10.5281/zenodo.823200 
* Romero, J., Phillips, E. Fatica, M., Spiga, F., Giannozzi, P.: _A performance study of Quantum ESPRESSO's PWscf code on multi-core and GPU systems_, 8th IEEE International Workshop on Performance Modeling, Benchmarking and Simulation of High Performance Computer Systems (PMBS17), Lecture Notes in Computer Science, Springer, Denver (2017)


### Benchmarks

Benchmarks will be collected in a separate [repository](https://github.com/romerojosh/qe-gpu-benchmarks) 

### License

All the material included in this distribution is free software; you can
redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 675 Mass
Ave, Cambridge, MA 02139, USA.
