# MTCR Interface Module

This module provides a Julia interface to MTCR (Multi-Temperature Collisional-Radiative Chemical Kinetics Solver), a Fortran code developed in the Nonequilibrium Gas and Plasma Dynamics Laboratory at the University of Colorado Boulder for detailed modeling of a reacting ionized gas in thermodynamic and chemical nonequilibrium.

## Overview

MTCR is a comprehensive thermochemical solver that handles:
- State-to-state vibrational and electronic kinetics
- Multi-temperature nonequilibrium effects
- Detailed chemical reaction mechanisms
- Collisional processes (electron-impact, heavy-particle)
- Radiative processes

## Quick Start

```julia
using HallThruster

# Initialize MTCR (requires compiled shared library)
HallThruster.MTCR.initialize_mtcr("/path/to/libmtcr.dylib")

# Run the reference nitrogen example
results = HallThruster.MTCR.nitrogen_10ev_example("/path/to/libmtcr.dylib")

# Access results
println("Final species densities: ", results.species_densities[:, end])
println("Final temperature: ", results.temperatures.tt[end])

# Clean up
HallThruster.MTCR.finalize_mtcr()
```

## Module Structure

- `mtcr_interface.jl` - Main module interface and exports
- `fortran_wrapper.jl` - Low-level Fortran bindings via ccall
- `data_conversion.jl` - Unit conversions and data format handling
- `mtcr_config.jl` - Configuration types and validation
- `mtcr_solver.jl` - High-level solver interface
- `IMPLEMENTATION_PLAN.md` - Detailed implementation documentation

## Configuration

Create configurations using the provided types:

```julia
config = HallThruster.MTCR.MTCRConfig(
    species = ["N", "N2", "N+", "N2+", "E-"],
    mole_fractions = [1e-20, 0.9998, 1e-20, 0.0001, 0.0001],
    total_number_density = 1.0e13,  # 1/cm³
    temperatures = HallThruster.MTCR.TemperatureConfig(750.0, 750.0, 750.0, 115000.0),
    time_params = HallThruster.MTCR.TimeIntegrationConfig(0.5e-5, 5.0, 1e3)
)

results = HallThruster.MTCR.solve_mtcr_0d(config)
```

## Requirements

### Fortran Library
- MTCR must be compiled as a shared library (libmtcr.dylib/so/dll)
- Requires Fortran compiler and MPI library
- See `IMPLEMENTATION_PLAN.md` for build instructions

### Julia Dependencies
- Base Julia ≥1.10
- DocStringExtensions.jl (already included in HallThruster.jl)

## Current Status

✅ **Complete**: Julia interface implementation
⚠️ **TODO**: Fortran shared library compilation
⚠️ **TODO**: Test case validation
⚠️ **TODO**: HallThruster.jl coupling

See `IMPLEMENTATION_PLAN.md` for detailed status and next steps.

## Example: 0D Nitrogen Plasma

The reference test case simulates a 0D nitrogen plasma reactor:

```julia
# This replicates /mtcr/examples/0D_Nitrogen_Te_10eV
results = HallThruster.MTCR.nitrogen_10ev_example("/path/to/libmtcr.dylib")

# Expected behavior:
# - Initial: 99.98% N2, 0.01% N2+, 0.01% E-
# - Evolution: Dissociation and ionization toward equilibrium
# - Temperatures: Tt=750K, Tv=750K, Te=115,000K (10 eV)
```

## API Reference

### Main Functions
- `initialize_mtcr(lib_path)` - Initialize MTCR system
- `solve_mtcr_0d(config)` - Run 0D simulation
- `nitrogen_10ev_example(lib_path)` - Reference test case
- `finalize_mtcr()` - Clean up resources

### Configuration Types
- `MTCRConfig` - Main configuration
- `TemperatureConfig` - Temperature settings
- `TimeIntegrationConfig` - Time integration parameters
- `MTCRResults` - Results container

### Utility Functions
- `validate_results(results)` - Validate simulation results
- `save_results(results, filename)` - Save results to CSV
- `get_molecular_weights(species)` - Get species molecular weights

## Error Handling

The interface includes comprehensive error handling:
- Configuration validation before simulation
- Fortran error catching and conversion to Julia exceptions
- Physical consistency checks on results
- Graceful degradation on library loading failures

## Units

- **Input**: SI units (kg/m³, K, Pa, s)
- **Internal**: CGS units (g/cm³, K, dyne/cm², s) - automatic conversion
- **Output**: SI units (kg/m³, K, Pa, s)

All unit conversions are handled automatically by the interface.

## Contributing

When modifying this module:
1. Update `IMPLEMENTATION_PLAN.md` with status changes
2. Add tests for new functionality
3. Maintain backward compatibility
4. Follow Julia style guidelines
5. Update documentation

## Support

For issues related to:
- **Julia interface**: Check this module's documentation
- **MTCR physics**: Refer to MTCR documentation in `/mtcr/docs/`
- **Integration**: See `IMPLEMENTATION_PLAN.md`
