# kicked-rotor-chaos
From KAM tori to quantum localization: Chirikov standard map + quantum kicked rotor (Floquet / split-operator / FFT).

This repository contains:
- **Fortran (gfortran)** code to simulate the classical standard map and compute diagnostics (Poincaré sections, Lyapunov/FTLE, MSD diffusion). 
- **Python (Jupyter)** code to simulate the quantum kicked rotor via a symmetric split-operator Floquet step using **unitary FFTs**, and to post-process/plot poster-ready figures. 

## Repository structure
- `src/`: source of simulations
  - `classical/`
    - `constants.f90`, `standard_map.f90`, `diagnostics.f90`, `main_standardmap.f90`, `Makefile`
    - `out_standardmap/` (generated)
  - `quantum/`
    - `QuantumFloquet.ipynb` (quantum simulation)
    - `out_quantum/` (generated)
- `notebooks/`: for generating the plots
  - `ClassicPlots.ipynb` (plots + poster figures)
  - `QuantumPlots.ipynb` (plots + poster figures)
  - `science.mplstyle`: used by the notebooks, style information
- `figs/`: poster figures and assets

## Quick start

### 1) Classical standard map (Fortran)
Build and run (bash):
`make`
`make run`

### 2) Quantum Kicked Rotor (Python)
Run all cells in a python environment.
