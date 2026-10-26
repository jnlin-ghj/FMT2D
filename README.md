# FMT2D

FMT2D (Fast Marching Tomography 2D) is a two-dimensional multi-phase traveltime
tomography package for ocean-bottom seismometer (OBS) wide-angle seismic data.
It solves the eikonal equation with the Fast Marching Method (FMM) and performs
simultaneous inversion of velocity structure and reflector geometry.

FMT2D is written in Fortran (core forward modeling and inversion) with Python
example scripts for generating models, observation geometries, running the
binaries, and plotting results.

## Features

- **Eikonal-based forward modeling.** A globally unique traveltime field is
  computed over the whole model space, avoiding shadow zones and multipathing
  of ray-by-ray tracing.
- **Five preset P-wave phases:**
  - `Pg`: refracted wave in the solid medium
  - `PsedP`: reflection from the base of the sedimentary layer
  - `PmP`: Moho reflection
  - `Pg_ss`: source-side water-layer multiple of Pg
  - `Pg_ss2`: receiver-side water-layer multiple of Pg
- **Simultaneous velocity–interface inversion** with four modes:
  - joint inversion
  - velocity-only inversion
  - interface-only inversion
  - automatic alternating inversion
- **Three smoothing regularization options:** Laplacian, Gaussian, and hybrid.
- **LSQR solver** with damping and smoothing regularization.
- **Python example scripts** reproducing the checkerboard tests in the paper.

## Repository structure

