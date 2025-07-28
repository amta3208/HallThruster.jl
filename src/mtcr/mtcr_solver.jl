"""
# MTCR Solver Module

This module provides the high-level interface for running MTCR simulations
from Julia, hiding the complexity of the Fortran interface and providing
a clean, Julia-native API.

## Key Functions

- `solve_mtcr_0d`: Main solver function for 0D simulations
- `initialize_mtcr`: Initialize MTCR system
- `finalize_mtcr`: Clean up MTCR system
- `nitrogen_10ev_example`: Run the 0D nitrogen example case

## Usage

```julia
# Create configuration
config = nitrogen_10ev_config()

# Run simulation
results = solve_mtcr_0d(config)

# Access results
println("Final species densities: ", results.species_densities[:, end])
println("Final temperature: ", results.temperatures.tt[end])
```
"""

using DocStringExtensions

# Global state tracking
const MTCR_INITIALIZED = Ref{Bool}(false)
const MTCR_NUM_SPECIES = Ref{Int32}(0)
const MTCR_NUM_DIMENSIONS = Ref{Int32}(0)

"""
$(SIGNATURES)

Initialize the MTCR system.

This function must be called before any MTCR calculations can be performed.
It sets up the Fortran library, initializes internal data structures, and
prepares the system for simulation.

# Arguments
- `lib_path::String`: Path to the MTCR shared library

# Returns
- `true` if initialization successful

# Throws
- `ErrorException` if initialization fails
"""
function initialize_mtcr(lib_path::String)
    if MTCR_INITIALIZED[]
        @warn "MTCR already initialized. Skipping initialization."
        return true
    end

    try
        # Set library path
        set_mtcr_lib_path!(lib_path)

        # Initialize the API
        num_species, num_dimensions = initialize_api_wrapper()

        # Store system information
        MTCR_NUM_SPECIES[] = num_species
        MTCR_NUM_DIMENSIONS[] = num_dimensions
        MTCR_INITIALIZED[] = true

        @info "MTCR initialized successfully" num_species num_dimensions
        return true

    catch e
        @error "Failed to initialize MTCR" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Finalize the MTCR system and clean up resources.

This function should be called when MTCR is no longer needed to properly
clean up memory and resources.
"""
function finalize_mtcr()
    if !MTCR_INITIALIZED[]
        @warn "MTCR not initialized. Nothing to finalize."
        return
    end

    try
        finalize_api_wrapper()
        MTCR_INITIALIZED[] = false
        MTCR_NUM_SPECIES[] = 0
        MTCR_NUM_DIMENSIONS[] = 0
        @info "MTCR finalized successfully"
    catch e
        @error "Error during MTCR finalization" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Check if MTCR is properly initialized.

# Returns
- `true` if MTCR is initialized and ready for use
"""
function is_mtcr_initialized()
    return MTCR_INITIALIZED[]
end

"""
$(SIGNATURES)

Solve a 0D MTCR simulation.

This is the main high-level interface for running MTCR simulations.
It handles all the complexity of data conversion, Fortran interfacing,
and result processing.

# Arguments
- `config::MTCRConfig`: Configuration for the simulation

# Returns
- `MTCRResults`: Results of the simulation

# Throws
- `ErrorException` if MTCR not initialized or simulation fails
"""
function solve_mtcr_0d(config::MTCRConfig)
    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_mtcr(lib_path) first.")
    end

    try
        @info "Starting MTCR 0D simulation" species=config.species

        # Convert configuration to initial conditions
        initial_state = config_to_initial_state(config)

        # Run the time integration
        results = integrate_0d_system(config, initial_state)

        @info "MTCR simulation completed successfully"
        return results

    catch e
        @error "MTCR simulation failed" exception=e
        return MTCRResults(
            Float64[], zeros(0, 0), (;), Float64[], nothing, false,
            "Simulation failed: $(string(e))",
        )
    end
end

"""
$(SIGNATURES)

Convert MTCRConfig to initial state vectors for MTCR.

# Arguments
- `config::MTCRConfig`: Configuration object

# Returns
- Named tuple with initial state vectors in CGS units
"""
function config_to_initial_state(config::MTCRConfig)
    # Get molecular weights
    molecular_weights = get_molecular_weights(config.species)

    # Convert mole fractions to mass densities (CGS units)
    number_density_cgs = convert_number_density_si_to_cgs(config.total_number_density)
    mass_densities_cgs = mole_fractions_to_mass_densities(
        config.mole_fractions, molecular_weights, number_density_cgs,
    )

    # Calculate initial total energy (simplified - assumes thermal equilibrium)
    # This is a placeholder - actual implementation would be more sophisticated
    rho_total = sum(mass_densities_cgs)

    # Estimate thermal energy (very simplified)
    # E = (3/2) * n * k * T for translational + electronic
    const BOLTZMANN_CGS = 1.380649e-16  # erg/K
    thermal_energy_per_volume = 1.5 * number_density_cgs * BOLTZMANN_CGS *
                                config.temperatures.Tt

    # Add electron thermal energy
    if "E-" in config.species
        electron_idx = findfirst(==(("E-")), config.species)
        if electron_idx !== nothing
            electron_number_density = config.mole_fractions[electron_idx] *
                                      number_density_cgs
            electron_thermal_energy = 1.5 * electron_number_density * BOLTZMANN_CGS *
                                      config.temperatures.Te
            thermal_energy_per_volume += electron_thermal_energy
        end
    end

    return (
        rho_sp = mass_densities_cgs,
        rho_etot = thermal_energy_per_volume,
        number_density = number_density_cgs,
        molecular_weights = molecular_weights,
    )
end

"""
$(SIGNATURES)

Integrate the 0D system over time.

# Arguments
- `config::MTCRConfig`: Configuration object
- `initial_state`: Initial state vectors

# Returns
- `MTCRResults`: Simulation results
"""
function integrate_0d_system(config::MTCRConfig, initial_state)
    # Time parameters
    dt = config.time_params.dt
    dtm = config.time_params.dtm
    tlim = config.time_params.tlim
    nstep = config.time_params.nstep

    # Initialize storage arrays
    n_output_steps = Int(ceil(tlim / dtm)) + 1
    n_species = length(config.species)

    time_points = Float64[]
    species_densities = zeros(n_species, 0)
    temperatures_tt = Float64[]
    temperatures_te = Float64[]
    temperatures_tv = Float64[]
    total_energies = Float64[]

    # Current state
    current_time = 0.0
    current_state = deepcopy(initial_state)
    output_time = 0.0

    # Store initial conditions
    push!(time_points, current_time)
    species_densities = hcat(species_densities, current_state.rho_sp)

    # Calculate initial temperatures
    temps = calculate_temperatures_wrapper(
        current_state.rho_sp, current_state.rho_etot,
    )
    push!(temperatures_tt, temps.tt)
    push!(temperatures_te, temps.teex)
    push!(temperatures_tv, temps.tvib)
    push!(total_energies, current_state.rho_etot)

    step = 0
    while current_time < tlim && step < nstep
        # Calculate source terms
        sources = calculate_sources_wrapper(
            current_state.rho_sp, current_state.rho_etot,
        )

        # Simple forward Euler integration (placeholder)
        # Real implementation would use the specified integration method
        current_state = (
            rho_sp = current_state.rho_sp .+ dt .* sources.drho_sp,
            rho_etot = current_state.rho_etot + dt * sources.drho_etot,
            number_density = current_state.number_density,
            molecular_weights = current_state.molecular_weights,
        )

        # Ensure non-negative densities
        current_state = (
            rho_sp = max.(current_state.rho_sp, 1e-30),
            rho_etot = current_state.rho_etot,
            number_density = current_state.number_density,
            molecular_weights = current_state.molecular_weights,
        )

        current_time += dt
        step += 1

        # Check if it's time to output
        if current_time >= output_time + dtm || current_time >= tlim
            push!(time_points, current_time)
            species_densities = hcat(species_densities, current_state.rho_sp)

            # Calculate temperatures
            temps = calculate_temperatures_wrapper(
                current_state.rho_sp, current_state.rho_etot,
            )
            push!(temperatures_tt, temps.tt)
            push!(temperatures_te, temps.teex)
            push!(temperatures_tv, temps.tvib)
            push!(total_energies, current_state.rho_etot)

            output_time = current_time

            # Progress reporting
            if step % 1000 == 0
                @info "Integration progress" time=current_time step=step
            end
        end
    end

    # Convert results back to SI units
    species_densities_si = convert_density_cgs_to_si.(eachcol(species_densities))
    species_densities_si_matrix = hcat(species_densities_si...)
    total_energies_si = [convert_energy_density_cgs_to_si(e) for e in total_energies]

    # Create results structure
    temperatures = (
        tt = temperatures_tt,
        te = temperatures_te,
        tv = temperatures_tv,
    )

    success = current_time >= tlim || step >= nstep
    message = success ? "Simulation completed successfully" : "Simulation terminated early"

    return MTCRResults(
        time_points,
        species_densities_si_matrix,
        temperatures,
        total_energies_si,
        nothing,  # source_terms - could be added later
        success,
        message,
    )
end

"""
$(SIGNATURES)

Run the 0D Nitrogen Te=10eV example case.

This function provides a convenient way to run the reference test case
that matches the MTCR example in `/mtcr/examples/0D_Nitrogen_Te_10eV`.

# Arguments
- `lib_path::String`: Path to MTCR shared library

# Returns
- `MTCRResults`: Results of the simulation

# Example
```julia
results = nitrogen_10ev_example("/path/to/libmtcr.dylib")
```
"""
function nitrogen_10ev_example(lib_path::String)
    # Initialize MTCR if not already done
    if !is_mtcr_initialized()
        initialize_mtcr(lib_path)
    end

    # Create configuration for the example case
    config = nitrogen_10ev_config()

    @info "Running 0D Nitrogen Te=10eV example case"
    @info "Configuration" species=config.species mole_fractions=config.mole_fractions
    @info "Temperatures" Tt=config.temperatures.Tt Te=config.temperatures.Te
    @info "Time parameters" dt=config.time_params.dt tlim=config.time_params.tlim

    # Run simulation
    results = solve_mtcr_0d(config)

    if results.success
        @info "Example simulation completed successfully"
        @info "Final conditions" time=results.time[end]

        # Print final species densities
        for (i, species) in enumerate(config.species)
            final_density = results.species_densities[i, end]
            @info "Final density" species=species density=final_density unit="kg/m³"
        end

        @info "Final temperatures" Tt=results.temperatures.tt[end] Te=results.temperatures.te[end]
    else
        @error "Example simulation failed" message=results.message
    end

    return results
end

"""
$(SIGNATURES)

Validate simulation results for physical consistency.

# Arguments
- `results::MTCRResults`: Simulation results to validate

# Returns
- `true` if results pass validation, `false` otherwise
"""
function validate_results(results::MTCRResults)
    if !results.success
        @warn "Simulation was not successful"
        return false
    end

    # Check for negative densities
    if any(results.species_densities .< 0)
        @warn "Negative species densities found"
        return false
    end

    # Check for NaN or Inf values
    if any(isnan.(results.species_densities)) || any(isinf.(results.species_densities))
        @warn "NaN or Inf values found in species densities"
        return false
    end

    if any(isnan.(results.temperatures.tt)) || any(isinf.(results.temperatures.tt))
        @warn "NaN or Inf values found in temperatures"
        return false
    end

    # Check temperature ranges (should be positive and reasonable)
    if any(results.temperatures.tt .<= 0) || any(results.temperatures.tt .> 1e6)
        @warn "Unreasonable translational temperatures found"
        return false
    end

    if any(results.temperatures.te .<= 0) || any(results.temperatures.te .> 1e6)
        @warn "Unreasonable electron temperatures found"
        return false
    end

    @info "Results validation passed"
    return true
end

"""
$(SIGNATURES)

Save MTCR results to file.

# Arguments
- `results::MTCRResults`: Results to save
- `filename::String`: Output filename (CSV format)

# Returns
- `true` if save successful
"""
function save_results(results::MTCRResults, filename::String)
    try
        # Prepare data for CSV output
        n_times = length(results.time)
        n_species = size(results.species_densities, 1)

        # Create header
        header = ["time", "total_energy", "T_trans", "T_electron", "T_vib"]
        for i in 1:n_species
            push!(header, "species_$(i)_density")
        end

        # Create data matrix
        data = zeros(n_times, length(header))
        data[:, 1] = results.time
        data[:, 2] = results.total_energy
        data[:, 3] = results.temperatures.tt
        data[:, 4] = results.temperatures.te
        data[:, 5] = results.temperatures.tv

        for i in 1:n_species
            data[:, 5 + i] = results.species_densities[i, :]
        end

        # Write to file
        open(filename, "w") do io
            # Write header
            println(io, join(header, ","))

            # Write data
            for row in eachrow(data)
                println(io, join(row, ","))
            end
        end

        @info "Results saved successfully" filename=filename
        return true

    catch e
        @error "Failed to save results" filename=filename exception=e
        return false
    end
end
