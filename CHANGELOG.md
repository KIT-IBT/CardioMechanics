# Changelog

## How to Use This Changelog
- Pending changes for the next release go under **Unreleased**.
- Released changes are grouped under their version heading.
- Each section uses subsections **Added** / **Changed** / **Fixed** / **Deleted** (and **Known Issues** for Unreleased).
- Empty subsections may be omitted in released sections; in Unreleased they are kept as placeholders.


## Unreleased

### Added
- Inverse problem (active-stress estimator): given a target deformed surface over time, `Solver.Type=ActiveStressEstimator` recovers the scalar active stress in every element at each timestep via Gauss-Newton with Tikhonov regularization and a direct MUMPS solve. Comprises the solver `CBSolverActiveStressEstimator`, the tension model `TensionEstimator`, the plugin `Solver.Plugins.PointsCtrl`, the active-stress Jacobian for T4 elements, and the `ExtractSurfaceNodesFromVTU` tool. The estimator runs serially. Reinstates a feature of the pre-modernization code base; the forward problem is unchanged.
- `examples/inverseEllipsoid`: truncated-ellipsoid ventricle running the full forward -> target-extraction -> inverse round-trip, with a regression test (`tests/cardiomechanics/test_inverse.py`) covering the same chain at a shortened stop time.
- `tools/python/CreateTargetSurfaces.py`: reusable CLI that converts forward-run VTUs into the binary target-surface format consumed by the estimator.
- Unit tests for the pure-python tool logic under `tests/python/` (no binaries required).
- CardioMechanics copies the settings file it was started with next to the results, as `<Export.Prefix>_settings.xml`, so an output directory always records the parameters that produced it. Note that `-parameter` overrides given on the command line are not reflected in the copy.
- Manual sections covering the inverse problem: the regularized Gauss-Newton formulation and its regularization terms in "Mathematical Model", the `ActiveStressEstimator` settings and the `PointsCtrl` plugin in "Simulation Framework", the `ExtractSurfaceNodesFromVTU` tool and the `CreateTargetSurfaces` script in "Tools", and the target-surface file format in "File Formats".

### Changed
- Corrected the statement in the manual that two solver classes are available, which no longer held once `ActiveStressEstimator` was added.
- `docs/BUILD.md` and the CI PETSc image now recommend and use OpenBLAS (`--download-openblas`) instead of reference BLAS. Results are unaffected; the inverse problem is roughly 7x faster, pure mechanics roughly 2x, and EP-dominated runs largely unchanged.
- `tools/python/VTK2tetgen.py` rewritten to be Python 3 compatible and importable: the CLI and conversion moved into `main()` under a `__main__` guard and the `vtk` import is deferred, so the geometry helpers can be imported without VTK. Mesh output is unchanged; `.bases` differs only in whitespace between the row index and the values.

### Fixed
- Corrected the sign of the master-to-target gap vector in `CBContactHandling`. The gap was taken as `(ip - p).Norm()`, which discards the sign of the signed distance along the master normal, so the vector pointed the wrong way whenever the target lay on the negative-normal side. Only the estimator consumes this vector; the forward contact force computes its own distance and is unaffected.
- `CBDataFromFile`'s default constructor left `startTime_` and `period_` uninitialized, which produced NaN sampled values on the default-constructed path used by `CBPointsCtrl`.

### Known Issues
- The active-stress estimator is serial only and aborts if launched under `mpirun`.


## Release 1.1

Build-system release. Existing build workflows will not carry over: see
[docs/BUILD.md](docs/BUILD.md) for the new setup.

### Added
- Regression test suite under `tests/`: pytest golden-file characterization tests covering `CellModelTest` (all ionic models, including Land17-coupled runs), the Land 2015 mechanics benchmark (Problem 1), and the EM01 electromechanics pipeline (`BidomainMatrixGenerator`, `acCELLerate`, and the coupled `CardioMechanics` run). Comparisons use numerical tolerances rather than byte-exact matching.
- Install rules for all executables, so `cmake --install <build-dir> --prefix <dir>` copies the binaries into `<dir>/bin/`. Note that `CMAKE_INSTALL_PREFIX` defaults to `/usr/local`, so `--prefix` should always be passed explicitly.
- Warning baseline `-pedantic -Wall -Werror=vla -Wno-extra-semi`, applied at directory scope so that per-target `-Wno-*` opt-outs in legacy directories take precedence.

### Changed
- Modernized the CMake build system: removed `cmake/IBTDefault.cmake` and its `ka*` macros in favor of standard `target_*()` commands; bumped to CMake 3.20 and C++17; switched PETSc discovery to pkg-config (only `PETSC_DIR` env var required, `kaRootDir` no longer needed); added `CMakePresets.json`; binaries now land in `_build/<preset>/bin/`.
- Restructured `CMakePresets.json`: `default` is now a hidden base preset; `release` and `debug` are the user-facing configure/build presets, each with its own build directory (`_build/release/`, `_build/debug/`); `CMAKE_EXPORT_COMPILE_COMMANDS` enabled by default.
- The bundled model parameter files under `electrophysiology/data` are now resolved via `CM_SOURCE_DIR`, baked in at build time. The `kaRootDir` environment variable is no longer required at runtime, though it still overrides the compiled-in path when set.
- `find_package(VTK)` now requests only the components actually used instead of the whole package, which drops the Qt requirement from the build.
- Docker images: the CI image builds via the `release` preset and no longer installs Qt. `Dockerfile-thirdparty-petsc` source-builds PETSc v3.24.0, while Open MPI and VTK come from apt (`libopenmpi-dev`, `libvtk9-dev`).
- Resolved compiler warnings across the mechanics and electrophysiology sources exposed by the new warning baseline.

### Deleted
- `installRequirements.sh`, superseded by the dependency instructions in `docs/BUILD.md`.

### Fixed
- Suppress spurious PETSc "options left" warning for CardioMechanics' own CLI flags (`-settings`, `-verbose`, etc.). PETSc 3.21+ reports unused options at finalize by default; since the app parses its flags directly from `argv` rather than through the PETSc options API, they were never marked used. The fix removes them from PETSc's options database after parsing via a new `DCCtrl::ClearOption` abstraction backed by `PetscOptionsClearValue`.
- Renamed `typedef DCCtrlPETSc Petsc` to `typedef DCCtrlPETSc DCPetsc` in `mechanics/src/DCTK/DCCtrlPETSc.h` and updated all call sites. The name `Petsc` collided with PETSc's own `::Petsc` C++ namespace, which is exposed in private headers included by debug PETSc builds, causing compilation failures.
- Fixed segfaults on every parallel run (`np` ≥ 2). The target-based build linked only `-lpetsc` and omitted PETSc's private dependencies (MUMPS, ScaLAPACK, BLAS/LAPACK, MPI-Fortran); their symbols then bound to other providers pulled in transitively (e.g. VTK's Accelerate BLAS), corrupting MUMPS's distributed-RHS solve. The build now links PETSc's full pkg-config closure and adds PETSc's library directory to the rpath so `@rpath` dependencies (e.g. HYPRE) resolve at runtime.
- Link MPI into `PETScVec2VTK`, which failed to build without it.


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
