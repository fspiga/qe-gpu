- [Report new Bug Issue](.github/ISSUE_TEMPLATE.md)
- [Guidelines for Contributing](CONTRIBUTING.md)
- [Pull Request Template](.github/PULL_REQUEST_TEMPLATE.md)
- [Project License](License)


## GPU-accelerated Quantum ESPRESSO (QE-GPU)

This is an open-source custom version of Quantum ESPRESSO v6.x with embedded GPU support based on CUDA FORTRAN. This product has been made possible thanks to the effort of the [NVIDIA](http://www.nvidia.com/page/home.html) HPC Software and Benchmarks Group. This version is maintained by [Filippo Spiga](https://github.com/fspiga) and hosted on University of Cambridge Research Software Engineering GitHub page. Partial support was provided by E4 Computer Engineering SpA via the European PRACE Pre-Commercial Procurement project (Phase 3).

To contribute please refer to the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md)


### Requirements

The feely available compiler suite [PGI Community Edition 17.4](http://www.pgroup.com/products/community.htm) is required to use QE-GPU. It containes CUDA SDK 8.0 and pre-built Open MPI for parallel execution (follow the [PGI Instalation Guide](http://www.pgroup.com/doc/pgiinstall174.pdf) how to install it). **No other compilers are supported**.

You need data-centre grade NVIDIA TESLA Kepler (K20, K20, K40, K80) and NVIDIA TESLA Pascal (P100) compute GPUs. No other cards are supported. NVIDIA TESLA P100 is strongly recommend for its memory capacity and performance.

This version of QE-GPU run **exclusively** in parallel, Open MPI is required and also Intel MKL.

### Installation

Installation process follow the same procedure as Quantum ESPRESSO suite:

```
./configure --enable-gpu
make pw
```

"make" alone prints a list of acceptable targets with GPU acceleration. Binaries go in "bin/". Additional configure options are made available to customize the building process:

* `--enable-gpu=<kepler|pascal>` to enable GPU support (default: on, "pascal" architecture selected)

The QE-GPU package has been reduced in size to the minimum essential. For more information, please refer to the general documentation provided with the Quantum ESPRESSO suite or visit the web site [http://www.quantum-espresso.org/](http://www.quantum-espresso.org/)


### Cite as:

A Zenodo citation will be uploaded soon


### License

All the material included in this distribution is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
