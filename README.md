# ColloidalCantilever

This repository contains Python code for simulating the drying of a colloidal suspension that has been deposited onto a cantilever.  As the film dries, mechanical stresses are generated.  The stress in the film is transmitted to the cantilever, causing it to deform.  

The film and the cantilever are assumed to be two dimensional. 
The generation of stress in the film is described using nonlinear poroelasticity in the thin-film limit.  The cantilever is described using a nonlinear beam theory.  The model is described in Hennessy et al (2026), *Stress development in attractive colloidal suspensions: the role of drying and gelation rates*, which is currently under review.

The model here uses an evaporation flux that is computed by solving
for the vapour concentration field around the film.  The problem for the vapour concentration is solved using the finite element method.

## Dependencies 

The code for simulating film drying and the deformation of the cantilver requires SciPy, NumPy, and Matplotlib to run.

The code for computing the evaporation flux requires Gmsh and legacy FEniCS.

## Getting started

The code is configured to simulate the drying of a colloidal suspension across four evaporation rates (relative humidity values).  The deflection of the substrate and the position of the solidification front are plotted.

The code can be run using the command

```
python3 simulate_residual.py
```

## Computing the evaporation flux

To run the code that computes the evaporation flux, first the mesh must be generated.  All of the files for mesh generation can be found in the [mesh](/mesh/) folder.  When the mesh is generated, the code in `compute_evap_flux.py` can be run.