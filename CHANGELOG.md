# Changelog
All notable changes to this project will be documented in this file. This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixes
- Workaround for compiler bug on POWER systems (see 92327a54)

## [1.0.0] - 2017-11-03

### Added
- GPU support using CUDA FORTRAN to PWscf v6.1, k-point calculation only

### Changed
- Simplified build system 

### Removed
- Auxiliary CPU-only packages: XSpectra, TDDFPT, PWCOND, PP, NEB, GWW, EPW, COUPLE, atomic
- Autoconf scripts
- Example files
- Documentation (please refer to official release)
