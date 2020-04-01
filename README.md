# Project 1: Molecular dynamics simulation of Argon atoms

This project simulates the dynamics of Argon atoms under the influence of a Lennard-Jones potential.
Currently, two particles moving in a box with periodic boundary conditions are simulated for a
user-defined number of timesteps. The particle positions are plotted in a separate window as
the simulation progresses.

Requisites: numpy,matplotlib

Authors:
- Brennan Undseth
- Ludwig Hendl
- Isacco Gobbi

The final version of the simulator is centered around the simulation_func_vector.py
script, which comprises the core of the code. Parameters and batch simulations are
run by manipulating parameters in the main.py script, and the diffusion_processor.py
file is a helper script for compiling diffusion data.