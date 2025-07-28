# MTCR-HallThruster.jl Integration Implementation Plan

This document outlines the comprehensive plan for integrating MTCR (Multi-Temperature Chemical Reactor) Fortran code with HallThruster.jl to enable detailed thermochemical nonequilibrium modeling within Hall thruster simulations.

## Project Overview

**HallThruster.jl** is a Julia package for simulating Hall effect thrusters, focusing on plasma physics, collision processes, and thruster performance. It has a well-structured modular architecture with separate modules for physics, collisions, simulation, etc.

**MTCR** is a Fortran codebase for modeling thermochemical nonequilibrium in high-temperature gas flows. It provides detailed chemistry modeling with state-to-state kinetics, vibrational/electronic excitation, and various collision processes.

The target example case `/mtcr/examples/0D_Nitrogen_Te_10eV` is a zero-dimensional nitrogen plasma reactor with:
- 5 species: N, N2, N+, N2+, E-
- Multiple electronic states for each species
- Electron temperature of 10 eV (~115,000 K)
- Various thermochemical processes (dissociation, ionization, charge exchange)

## Architecture Overview

The integration follows a modular design within the existing HallThruster.jl structure:

```
src/
├── mtcr/                    # New MTCR interface module
│   ├── mtcr_interface.jl    # Main interface module
│   ├── fortran_wrapper.jl   # Low-level Fortran bindings
│   ├── data_conversion.jl   # Julia ↔ Fortran data conversion
│   ├── mtcr_config.jl       # Configuration management
│   ├── mtcr_solver.jl       # High-level solver interface
│   └── IMPLEMENTATION_PLAN.md # This document
```

## Implementation Phases

### Phase 1: Fortran Library Preparation ⚠️ TODO

**Objective**: Build MTCR as a shared library compatible with Julia's ccall interface.

**Tasks**:
1. **Modify MTCR Makefile**
   - Add target to build shared library (.dylib on macOS, .so on Linux, .dll on Windows)
   - Ensure proper symbol export for API functions
   - Handle MPI dependencies appropriately for library usage

2. **Test Library Compilation**
   ```bash
   cd mtcr/source
   make api  # Should create libmtcr.a
   # Need to modify to create shared library
   ```

3. **Verify API Function Exports**
   - Ensure functions from `interface.f90` are properly exported
   - Test basic library loading in Julia
   ```julia
   lib = dlopen("path/to/libmtcr.dylib")
   # Verify symbols are available
   ```

**Dependencies**:
- Fortran compiler (gfortran or Intel Fortran)
- MPI library (for MTCR's parallel capabilities)
- ODEPACK library (already included in MTCR)

**Deliverables**:
- Modified Makefile with shared library target
- Compiled shared library (libmtcr.dylib/so/dll)
- Basic Julia loading test

### Phase 2: Julia-Fortran Interface ✅ COMPLETE

**Objective**: Create robust low-level bindings between Julia and MTCR Fortran code.

**Completed Components**:

1. **Low-level bindings** (`fortran_wrapper.jl`) ✅
   - `ccall` interfaces to MTCR API functions
   - Memory management for array passing
   - Error handling for Fortran exceptions
   - Key functions implemented:
     - `initialize_api_wrapper()`
     - `finalize_api_wrapper()`
     - `calculate_sources_wrapper()`
     - `calculate_temperatures_wrapper()`
     - `get_species_names_wrapper()`

2. **Data conversion utilities** (`data_conversion.jl`) ✅
   - Unit conversions (SI ↔ CGS)
   - Array format transformations
   - Species data validation
   - Key functions:
     - `convert_state_si_to_cgs()` / `convert_state_cgs_to_si()`
     - `mole_fractions_to_mass_densities()`
     - `validate_species_data()`
     - `create_species_mapping()`

**Status**: Implementation complete, needs testing with actual Fortran library.

### Phase 3: High-Level Interface ✅ COMPLETE

**Objective**: Provide user-friendly Julia interface hiding Fortran complexity.

**Completed Components**:

1. **Configuration management** (`mtcr_config.jl`) ✅
   - Structured configuration types:
     - `MTCRConfig`: Main configuration
     - `TemperatureConfig`: Temperature settings
     - `TimeIntegrationConfig`: Time integration parameters
     - `PhysicsConfig`: Physics modeling options
     - `ProcessConfig`: Process flags
   - Input validation and default parameters
   - MTCR input file generation
   - `nitrogen_10ev_config()` for test case

2. **Solver interface** (`mtcr_solver.jl`) ✅
   - High-level simulation functions:
     - `solve_mtcr_0d()`: Main 0D solver
     - `initialize_mtcr()` / `finalize_mtcr()`: System management
     - `nitrogen_10ev_example()`: Reference test case
   - Result validation and file I/O
   - Error handling and logging

3. **Main interface module** (`mtcr_interface.jl`) ✅
   - Module organization and exports
   - Documentation and usage examples
   - Integration point with HallThruster.jl

**Status**: Implementation complete, ready for testing.

### Phase 4: Test Case Implementation ⚠️ TODO

**Objective**: Replicate the 0D Nitrogen Te=10eV example and validate results.

**Tasks**:

1. **Build and Test Fortran Library**
   - Complete Phase 1 tasks
   - Verify library can be loaded and initialized

2. **Implement Reference Test Case**
   ```julia
   using HallThruster.MTCR

   # Initialize MTCR
   initialize_mtcr("/path/to/libmtcr.dylib")

   # Run reference case
   results = nitrogen_10ev_example("/path/to/libmtcr.dylib")

   # Validate results
   @assert validate_results(results)
   ```

3. **Results Validation**
   - Compare with original MTCR output from `/mtcr/examples/0D_Nitrogen_Te_10eV/`
   - Verify species evolution matches expected behavior
   - Check temperature evolution and energy conservation
   - Validate final equilibrium state

4. **Performance Benchmarking**
   - Measure interface overhead
   - Compare performance with standalone MTCR
   - Optimize data conversion bottlenecks

**Expected Results**:
- Initial conditions: N2 dominant (99.98%), small N2+ and E- (0.01% each)
- Temperature evolution: Tt=750K, Tv=750K, Te=115,000K
- Chemical evolution toward equilibrium
- Final state with increased dissociation and ionization

**Deliverables**:
- Working test case function
- Validation report comparing with MTCR reference
- Performance benchmark results

### Phase 5: Integration with HallThruster.jl ⚠️ TODO

**Objective**: Enable MTCR usage within HallThruster.jl simulations.

**Tasks**:

1. **HallThruster.jl Integration**
   - Add MTCR module to main HallThruster.jl exports
   - Update Project.toml if needed
   - Create integration examples

2. **Coupling Interface Development**
   ```julia
   # Extract data from HallThruster simulation
   function extract_mtcr_conditions(ht_solution, cell_index)
       # Convert HallThruster state to MTCR input format
       # Handle species mapping
       # Extract temperatures and densities
   end

   # Apply MTCR results back to HallThruster
   function apply_mtcr_results(ht_solution, mtcr_results, cell_index)
       # Update species densities
       # Update energy/temperature
       # Maintain conservation laws
   end
   ```

3. **Species Mapping**
   - Map between HallThruster.jl and MTCR species conventions
   - Handle different propellants (Xe, Kr, Ar, N2)
   - Validate species conservation

4. **Documentation and Examples**
   - User guide for MTCR integration
   - Example coupled simulations
   - API reference documentation

**Deliverables**:
- Coupling interface functions
- Species mapping utilities
- Integration examples
- Documentation

## Technical Challenges and Solutions

### 1. Memory Management
**Challenge**: Julia and Fortran have different memory models.
**Solution**:
- Use `Ref` types for scalar passing
- Ensure contiguous array layout with `Array{Float64}`
- Proper cleanup in finalization functions

### 2. Unit Systems
**Challenge**: HallThruster.jl uses SI, MTCR uses CGS.
**Solution**:
- Comprehensive conversion utilities in `data_conversion.jl`
- Automatic unit conversion in interface functions
- Clear documentation of unit expectations

### 3. Species Mapping
**Challenge**: Different species naming conventions.
**Solution**:
- Species mapping dictionary in `create_species_mapping()`
- Validation functions to ensure consistency
- Flexible mapping for different propellants

### 4. Error Handling
**Challenge**: Fortran errors need proper Julia exception handling.
**Solution**:
- Wrapper functions catch and convert Fortran errors
- Comprehensive validation before Fortran calls
- Graceful degradation on failure

### 5. Performance
**Challenge**: Interface overhead could impact performance.
**Solution**:
- Minimize data copying with in-place operations
- Batch operations when possible
- Profile and optimize critical paths

## API Design

### High-Level Interface
```julia
# Configuration
config = MTCRConfig(
    species = ["N", "N2", "N+", "N2+", "E-"],
    mole_fractions = [1e-20, 0.9998, 1e-20, 0.0001, 0.0001],
    total_number_density = 1.0e13,
    temperatures = TemperatureConfig(750.0, 750.0, 115000.0),
    time_params = TimeIntegrationConfig(0.5e-5, 5.0, 1e3)
)

# Simulation
results = solve_mtcr_0d(config)

# Results access
final_densities = results.species_densities[:, end]
temperature_evolution = results.temperatures.tt
```

### Integration Interface
```julia
# For coupling with HallThruster.jl
function couple_with_mtcr(ht_solution, mtcr_config)
    # Extract conditions from HallThruster solution
    # Run MTCR for each cell/region
    # Update HallThruster solution with MTCR results
end
```

## Testing Strategy

### Unit Tests
- [ ] Data conversion functions
- [ ] Configuration validation
- [ ] Species mapping
- [ ] Unit conversions

### Integration Tests
- [ ] Fortran library loading
- [ ] API function calls
- [ ] Memory management
- [ ] Error handling

### Validation Tests
- [ ] 0D Nitrogen example replication
- [ ] Comparison with MTCR reference results
- [ ] Physical consistency checks
- [ ] Conservation law verification

### Performance Tests
- [ ] Interface overhead measurement
- [ ] Memory usage profiling
- [ ] Scaling with problem size
- [ ] Comparison with standalone MTCR

## Dependencies and Requirements

### Julia Dependencies
- Base Julia (≥1.10)
- DocStringExtensions.jl (for documentation)
- No additional Julia packages required

### System Dependencies
- Fortran compiler (gfortran or Intel Fortran)
- MPI library (OpenMPI or MPICH)
- MTCR source code and databases

### Build Requirements
- Make or CMake for building shared library
- Compatible C/Fortran compiler toolchain
- Proper library path configuration

## Deployment and Distribution

### Library Distribution
- MTCR shared library needs to be built on target system
- Consider providing pre-built binaries for common platforms
- Document build process for different operating systems

### Integration with HallThruster.jl
- MTCR module included as part of HallThruster.jl
- Optional dependency (graceful degradation if not available)
- Clear documentation for enabling MTCR features

## Future Enhancements

### Short Term
- [ ] Support for additional propellants (Ar, Kr, Xe)
- [ ] 1D MTCR integration
- [ ] Advanced time integration methods
- [ ] Parallel processing support

### Medium Term
- [ ] Adaptive chemistry (enable/disable based on conditions)
- [ ] Reduced chemistry models for performance
- [ ] Integration with other chemistry solvers
- [ ] GPU acceleration possibilities

### Long Term
- [ ] Multi-dimensional MTCR coupling
- [ ] Real-time chemistry adaptation
- [ ] Machine learning surrogate models
- [ ] Advanced visualization tools

## Status Summary

| Phase | Status | Completion |
|-------|--------|------------|
| Phase 1: Fortran Library | ⚠️ TODO | 0% |
| Phase 2: Julia-Fortran Interface | ✅ COMPLETE | 100% |
| Phase 3: High-Level Interface | ✅ COMPLETE | 100% |
| Phase 4: Test Case Implementation | ⚠️ TODO | 0% |
| Phase 5: HallThruster Integration | ⚠️ TODO | 0% |

**Overall Progress**: 40% complete

## Next Steps

1. **Immediate (Phase 1)**:
   - Modify MTCR Makefile to build shared library
   - Test library compilation on target platform
   - Verify API function exports

2. **Short Term (Phase 4)**:
   - Implement and test reference case
   - Validate results against MTCR output
   - Performance benchmarking

3. **Medium Term (Phase 5)**:
   - Develop coupling interface
   - Create integration examples
   - Documentation and testing

The foundation for MTCR integration is now in place with a complete Julia interface. The next critical step is building the Fortran shared library and validating the interface with the reference test case.
