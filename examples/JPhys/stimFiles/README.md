# Stimulation

The atria were stimulated in the sinus node, which corresponds to all points located in material 73 of the EP mesh.

The fascicular approach was used to stimulate the ventricles.
This means that nodes within a certain radius of five root points on the endocardial layer of the ventricles were stimulated.
These root points are given in terms of [Cobiveco](https://github.com/KIT-IBT/Cobiveco) coordinates in the following table.


| ab | tm | rt | tv | radius (mm) | $\Delta$ tm | LAT (s) | description |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.46 | 0.99 | 0.101 | 0 | 3 | 0.05 | 0.16  | lv posterior |
| 0.90 | 0.99 | 0.528 | 0 | 3 | 0.05 | 0.16  | lv anterior |
| 0.33 | 0.99 | 0.821 | 0 | 3 | 0.05 | 0.16  | lv septal |
| 0.50 | 0.99 | 0.875 | 1 | 3 | 0.05 | 0.165 | rv septal |
| 0.37 | 0.99 | 0.432 | 1 | 3 | 0.05 | 0.172 | rv moderator band |

At the point of doing the work for this publication, Cobiveco required the base to be cut off.
The file [clippedVentricles.vtp](./clippedVentricles.vtp) included here shows the chosen cutoff and contains the class labels required to run the Cobiveco algorithm.
