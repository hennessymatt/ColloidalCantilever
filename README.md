# ColloidalCantilever

This repository contains Python code for simulating the drying of a colloidal suspension that has been deposited onto a cantilever.  As the film dries, mechanical stresses are generated.  The stress in the film is transmitted to the cantilever, causing it to deform.  

The film and the cantilever are assumed to be two dimensional. 
The generation of stress in the film is described using nonlinear poroelasticity in the thin-film limit.  The cantilever is described using a nonlinear beam theory.  The model is described in Hennessy et al (2026), *The fracture toughness of attractive colloidal suspensions is controlled by the relative rate of gelation to drying*, which is currently under review.



## Dependencies 

The code requires SciPy, NumPy, and Matplotlib to run.

## Getting started

The code is configured to simulate the drying of a colloidal suspension across four evaporation rates (relative humidity values).  The deflection of the substrate and the position of the solidification front are plotted.

The code can be run using the command

```
python3 simulate.py
```