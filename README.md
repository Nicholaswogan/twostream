# twostream

This repository contains the two-stream radiative transfer algorithm from [Toon et al. (1989)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JD094iD13p16287) written in Python, C++, and Fortran.

Each implementation has the following inputs:

- `nz` - The number of vertical layers in the atmosphere.
- `tau` - The optical depth of each layer in the atmosphere. The first element is the top layer and the last element is the surface layer.
- `w0` - The single scattering albedo of each layer in the atmosphere.
- `u0` - cosine of the solar zenith angle.
- `Rsfc` - Albedo of the surface.

and the following output:

- `amean` - The mean intensity at the edges of each layer in the atmosphere. The first element is the top edge op the top layer

Every term is unit-less, and the flux at the top of the atmosphere is taken to be 1.

I use this algorithm to compute photolysis rates in photochemical modeling. [Here is a link to an example](https://github.com/Nicholaswogan/Photochem/blob/11fcdda980e13983d816e9f8e42a4a4ea123be12/src/photochem.f90#L364).

