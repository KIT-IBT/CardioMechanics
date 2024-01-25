# Single cell simulations

This folder contains the parameter files for the ionic models used in scenarios I to V of the paper.
To find stable initial values for the ODEs of the ionic models, single cell simulations using 0.8s basic cycle length were run for 1000 cycles.
You can replicate this by running the bash script `./runCellModels.sh` in this directory for the different scenarios by changing the parameter `MODE` in the script to the associated scenario, e.g. by using `MODE="I"`.

The $I_{SAC}$ sensitivity study in the appendix can be reproduced by using the bash scripts `./runSensitivityAnalysis.sh` and `./runSensitivityAnalysis_E.sh`, respectively.
For these, a clamp file is used to define the different stretch intervals.
The structure of the file is given by 

| VAR | TIME (s) | VALUE |
| --- | --- | ---|
| stretch | 0.45 | 1.0 |
| ... | ... | ... |

The variable `VAR` can be set to `Vm` and `Cai` if you want to clamp transmembrane potential or intracellular calcium, respectively.