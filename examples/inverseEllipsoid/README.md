# Inverse problem: active-stress estimator (ellipsoid)

Given a target deformed surface over time, the estimator recovers the scalar
active stress in every element at each timestep (Gauss-Newton + Tikhonov
regularization, direct MUMPS solve). This example runs the full round-trip on a
truncated-ellipsoid ventricle (T4 mesh):

```
forwardT4.xml  ->  CreateTargetSurfaces.py  ->  inverseT4.xml
  (deform)          (extract target surface)     (recover active stress)
```

The forward run applies a known DoubleHill active tension and deforms the mesh.
`CreateTargetSurfaces.py` extracts the deformed surface nodes from each forward
VTU into the binary target format, and the inverse run recovers the per-element
active stress that reproduces those targets.

## Files
- `geometry/`: T4 mesh (`meshT4IP.node`, `meshT4IP.sur`, `meshT4.ele`, `meshT4.bases`).
- `forwardT4.xml`: forward run (`NewmarkBeta`, Guccione + DoubleHill tension).
- `inverseT4.xml`: inverse run (`Solver.Type=ActiveStressEstimator`, serial/direct).

Both use the modernized Guccione strain-energy convention `W = C/2 (e^Q - 1)`
with `C = 556`.

## Run

```sh
mkdir -p Results

# 1. forward: generate the target deformation (fast)
CardioMechanics -settings forwardT4.xml

# 2. extract the target surfaces (needs ExtractSurfaceNodesFromVTU on PATH)
python3 ../../tools/python/CreateTargetSurfaces.py \
    -results Results/forward_vtu \
    -node geometry/meshT4IP.node -sur geometry/meshT4IP.sur \
    -output TargetSurfaces -list TargetSurfaces.list -dt 0.01

# 3. inverse: recover the active stress (serial, direct MUMPS solve)
CardioMechanics -settings inverseT4.xml
```

The estimator runs serially (direct solve); do not launch it under `mpirun`.
`-dt` is the forward export timestep (`Export.TimeStep`, here `1e-2` s).

Recovered `ActiveStress` (cell data) and the deformed geometry are written to
`Results/Inverse_vtu/Inverse.<N>.vtu`. The regression test
`tests/cardiomechanics/test_inverse.py` runs this same chain shortened to
`StopTime=0.25`.
