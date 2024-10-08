# Distributed Generator Placement in IEEE 33-Bus Test System Using Equilibrium Optimizer

This repository provides a MATLAB implementation for the optimal placement of Distributed Generators (DGs) in the IEEE 33-bus test system using the Equilibrium Optimizer (EO) algorithm. The primary objective is to minimize power losses and improve voltage profiles across the distribution network.

### Features:
- Y-Bus Formulation: The repository includes a comprehensive function for calculating the Y-Bus matrix.
- Newton-Raphson Power Flow: A robust implementation of the Newton-Raphson method is provided for solving power flow equations efficiently in the distribution network.
- Constant Generators as DGs: The DGs modeled in this work are constant power generators. This assumption simplifies the analysis by focusing on consistent generation without variability.
- Equilibrium Optimizer (EO) Algorithm: The placement of DGs is optimized using the EO algorithm.

### Repository Contents:
- MATLAB scripts for Y-Bus matrix formation and Newton-Raphson power flow calculation.
- The Equilibrium Optimizer implementation tailored for the DG placement problem.
- Example simulations for the IEEE 33-bus test system showcasing the optimal placement of DGs.

