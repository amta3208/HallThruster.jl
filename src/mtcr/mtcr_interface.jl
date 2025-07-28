"""
# MTCR Interface Module

This module provides the main interface for integrating MTCR (Multi-Temperature Chemical Reactor)
Fortran code with HallThruster.jl. MTCR is a thermochemical nonequilibrium solver that handles
detailed chemistry modeling with state-to-state kinetics.

## Main Components

- `MTCRSolver`: High-level solver interface
- `MTCRConfig`: Configuration management
- Data conversion utilities for Julia ↔ Fortran interoperability
- Low-level Fortran bindings

## Usage Example

```julia
using HallThruster.MTCR

# Create configuration for 0D nitrogen plasma
config = MTCRConfig(
    species = ["N", "N2", "N+", "N2+", "E-"],
    temperatures = (Tt=750.0, Tv=750.0, Te=115000.0),
    number_density = 1.0e13,
    time_params = (dt=0.5e-5, tlim=1e3)
)

# Set initial conditions (mole fractions)
initial_conditions = [1e-20, 0.9998, 1e-20, 0.0001, 0.0001]

# Solve the system
results = solve_mtcr_0d(config, initial_conditions)
```
"""
module MTCR

using DocStringExtensions

# Include submodules
include("fortran_wrapper.jl")
include("data_conversion.jl")
include("mtcr_config.jl")
include("mtcr_solver.jl")

# Re-export main functionality
export MTCRConfig, MTCRResults
export solve_mtcr_0d, nitrogen_10ev_example
export initialize_mtcr, finalize_mtcr

end # module MTCR
