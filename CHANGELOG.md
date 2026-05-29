# Changelog

## How to Use This Changelog
- Pending changes for the next release go under **Unreleased**.
- Released changes are grouped under their version heading.
- Each section uses subsections **Added** / **Changed** / **Fixed** / **Deleted** (and **Known Issues** for Unreleased).
- Empty subsections may be omitted in released sections; in Unreleased they are kept as placeholders.


## Unreleased

### Added

### Changed
- Modernized the CMake build system: removed `cmake/IBTDefault.cmake` and its `ka*` macros in favor of standard `target_*()` commands; bumped to CMake 3.20 and C++17; switched PETSc discovery to pkg-config (only `PETSC_DIR` env var required, `kaRootDir` no longer needed); added `CMakePresets.json`; binaries now land in `_build/bin/`.
- Restructured `CMakePresets.json`: `default` is now a hidden base preset; `release` and `debug` are the user-facing configure/build presets, each with its own build directory (`_build/release/`, `_build/debug/`); `CMAKE_EXPORT_COMPILE_COMMANDS` enabled by default.

### Fixed
- Renamed `typedef DCCtrlPETSc Petsc` to `typedef DCCtrlPETSc DCPetsc` in `mechanics/src/DCTK/DCCtrlPETSc.h` and updated all call sites. The name `Petsc` collided with PETSc's own `::Petsc` C++ namespace, which is exposed in private headers included by debug PETSc builds, causing compilation failures.

### Known Issues


## Release 1.0

### Added
- OHaraRudyIso model: adds PKA phosphorylation to multiple different currents as specified in the paper by Heijman et al. 2011 (https://doi.org/10.1016/j.yjmcc.2011.02.007) to the OHaraRudy model. The fractions of PKA phosphorylation can be specified in the corresponding ev-File (...P_frac).
- Tomek model (https://doi.org/10.7554/eLife.48890) with the addition of a stretch activated current and troponin C coupling to the Land17 tension model (can be removed from the model via SAC and TRPN preprocessor macros in the parameter.h file).

### Fixed
- Changed init_d from e-132 to 0.0 in Tomek_endo.ev, Tomek_mid.ev and Tomek_epi.ev as simulation would stall otherwise.
- Changed PCa Multiplier in OHaraRudy mid ev-File from 1.8 to original value 2.5.

### Changed
- Initial values of the Land17 tension model can now be specified in the CardioMechanics config file.
- Local deformation energy can now be exported.
- Sequence of variables in the ev-Files of Tomek and OHaraIso now match the sequence of the output.
- Hardcoded Docker image name replaced with the GitHub variable GITHUB_REPOSITORY in lowercase so that the action can be run in forks.
- Docker actions: replaced specific commit pins with semantic versioning.

### Deleted
