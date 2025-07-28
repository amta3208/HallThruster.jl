"""
# Fortran Wrapper Module

This module provides low-level Julia bindings to the MTCR Fortran library using `ccall`.
It handles direct interfacing with the Fortran API functions defined in `interface.f90`.

## Key Functions

- `initialize_api_wrapper`: Initialize MTCR system
- `finalize_api_wrapper`: Clean up MTCR system
- `calculate_sources_wrapper`: Calculate thermochemical source terms
- `calculate_temperatures_wrapper`: Calculate temperatures from state
- `calculate_total_energy_wrapper`: Calculate total energy

## Memory Management

This module handles the complexities of passing data between Julia and Fortran,
including proper array allocation, memory management, and type conversions.
"""

# MTCR library state
const MTCR_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const MTCR_LIB_PATH = Ref{String}("")

"""
$(SIGNATURES)

Set the path to the MTCR shared library and load it.

# Arguments
- `path::String`: Path to the MTCR shared library file

# Throws
- `ErrorException`: If the library cannot be loaded
"""
function set_mtcr_lib_path!(path::String)
    # Close existing handle if open
    if MTCR_HANDLE[] != C_NULL
        Libdl.dlclose(MTCR_HANDLE[])
        MTCR_HANDLE[] = C_NULL
    end

    # Validate path exists
    if !isfile(path)
        error("MTCR library file not found: $path")
    end

    # Open new library
    try
        MTCR_HANDLE[] = Libdl.dlopen(path)
        MTCR_LIB_PATH[] = path  # Store the path for ccall usage
    catch e
        error("Failed to load MTCR library from $path: $(e.msg)")
    end
end

"""
$(SIGNATURES)

Get the handle to the loaded MTCR library.

# Returns
- `Ptr{Cvoid}`: Handle to the loaded library

# Throws
- `ErrorException`: If no library is loaded
"""
function get_mtcr_handle()
    if MTCR_HANDLE[] == C_NULL
        error("MTCR library not loaded. Call set_mtcr_lib_path!(path) first.")
    end
    return MTCR_HANDLE[]
end

"""
$(SIGNATURES)

Get the path to the loaded MTCR library.

# Returns
- `String`: Path to the loaded library

# Throws
- `ErrorException`: If no library is loaded
"""
function get_mtcr_lib_path()
    if MTCR_HANDLE[] == C_NULL
        error("MTCR library not loaded. Call set_mtcr_lib_path!(path) first.")
    end
    return MTCR_LIB_PATH[]
end

"""
$(SIGNATURES)

Check if the MTCR library is currently loaded.

# Returns
- `Bool`: True if library is loaded, false otherwise
"""
function is_mtcr_loaded()
    return MTCR_HANDLE[] != C_NULL
end

"""
$(SIGNATURES)

Close the MTCR library and free resources.
"""
function close_mtcr_library()
    if MTCR_HANDLE[] != C_NULL
        Libdl.dlclose(MTCR_HANDLE[])
        MTCR_HANDLE[] = C_NULL
        MTCR_LIB_PATH[] = ""  # Clear the path
    end
end

"""
$(SIGNATURES)

Initialize the MTCR API system.

# Arguments
- `num_species::Int32`: Number of species in the system
- `num_dimensions::Int32`: Number of spatial dimensions (default: 0)
- `case_path::String`: Path to directory containing input/ subdirectory (default: current directory)

# Throws
- `ErrorException`: If case_path doesn't exist, input file is missing, or Fortran call fails
"""
function initialize_api_wrapper(
        num_species::Int32; num_dimensions::Int32 = Int32(0), case_path::String = pwd(),)
    # Validate inputs
    if !isdir(case_path)
        error("Case path does not exist: $case_path")
    end

    input_file = joinpath(case_path, "input", "prob_setup.inp")
    if !isfile(input_file)
        error("Required input file not found: $input_file")
    end

    # Ensure output directory structure exists
    output_dir = joinpath(case_path, "output")
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Ensure required output subdirectories exist
    sources_dir = joinpath(output_dir, "sources")
    states_dir = joinpath(output_dir, "states")

    if !isdir(sources_dir)
        mkpath(sources_dir)
    end

    if !isdir(states_dir)
        mkpath(states_dir)
    end

    # Store current directory and change to case path
    original_dir = pwd()

    try
        cd(case_path)

        # Create references for the Fortran call
        # Note: Fortran interface expects (num_dimensions, num_species) order
        num_dimensions_ref = Ref{Int32}(num_dimensions)
        num_species_ref = Ref{Int32}(num_species)

        # Call Fortran function
        ccall((:initialize_api, get_mtcr_lib_path()), Cvoid,
            (Ref{Int32}, Ref{Int32}),
            num_dimensions_ref, num_species_ref,)

        # If we get here, the call succeeded
        return nothing

    catch e
        # Re-throw with more context about the failure
        if isa(e, Base.SystemError) || isa(e, ErrorException)
            error("MTCR initialization failed in directory $case_path: $(e.msg)")
        else
            rethrow(e)
        end
    finally
        # Always restore original directory
        cd(original_dir)
    end
end

"""
$(SIGNATURES)

Finalize the MTCR API system and clean up resources.
"""
function finalize_api_wrapper()
    handle = get_mtcr_handle()

    # TODO: Uncomment when Fortran library is available (Phase 1 complete)
    # ccall((:finalize_api, handle), Cvoid, ())
end

"""
$(SIGNATURES)

Get the maximum number of species supported by MTCR.
"""
function get_max_number_of_species_wrapper()
    handle = get_mtcr_handle()

    # TODO: Uncomment when Fortran library is available (Phase 1 complete)
    # return ccall((:get_max_number_of_species, handle), Int32, ())

    return Int32(10)  # Placeholder: reasonable maximum for testing
end

"""
$(SIGNATURES)

Calculate nonequilibrium source terms.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_etot::Float64`: Total energy density
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `rho_vx::Array{Float64,3}`: Vibrational state densities (optional)
- `rho_erot::Float64`: Rotational energy density (optional)
- `rho_eeex::Float64`: Electron-electronic energy density (optional)
- `rho_evib::Float64`: Vibrational energy density (optional)

# Returns
- Tuple of derivative arrays corresponding to input arrays
"""
function calculate_sources_wrapper(rho_sp::Vector{Float64},
        rho_etot::Float64;
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        rho_u::Union{Float64, Nothing} = nothing,
        rho_v::Union{Float64, Nothing} = nothing,
        rho_w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing,)
    handle = get_mtcr_handle()

    # Prepare output arrays
    drho_sp = similar(rho_sp)
    drho_etot = Ref{Float64}(0.0)

    # Handle optional arguments
    drho_ex = rho_ex !== nothing ? similar(rho_ex) : nothing
    drho_vx = rho_vx !== nothing ? similar(rho_vx) : nothing
    drho_erot = rho_erot !== nothing ? Ref{Float64}(0.0) : nothing
    drho_eeex = rho_eeex !== nothing ? Ref{Float64}(0.0) : nothing
    drho_evib = rho_evib !== nothing ? Ref{Float64}(0.0) : nothing

    # TODO: Uncomment when Fortran library is available (Phase 1 complete)
    # Call Fortran subroutine
    # Note: This is a simplified interface - the actual implementation will need
    # to handle the complex optional argument passing to Fortran
    # ccall((:calculate_nonequilibrium_sources, lib_path), Cvoid,
    #     (Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}),
    #     rho_sp, rho_etot, drho_sp, drho_etot,)

    # Placeholder: zero source terms for testing
    fill!(drho_sp, 0.0)
    drho_etot[] = 0.0

    return (drho_sp = drho_sp,
        drho_etot = drho_etot[],
        drho_ex = drho_ex,
        drho_vx = drho_vx,
        drho_erot = drho_erot !== nothing ? drho_erot[] : nothing,
        drho_eeex = drho_eeex !== nothing ? drho_eeex[] : nothing,
        drho_evib = drho_evib !== nothing ? drho_evib[] : nothing,)
end

"""
$(SIGNATURES)

Calculate temperatures from thermodynamic state.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities
- `rho_etot::Float64`: Total energy density
- Additional optional energy components

# Returns
- Named tuple with temperatures (tt, trot, teex, tvib, tex, tvx)
"""
function calculate_temperatures_wrapper(rho_sp::Vector{Float64},
        rho_etot::Float64;
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        rho_u::Union{Float64, Nothing} = nothing,
        rho_v::Union{Float64, Nothing} = nothing,
        rho_w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing,)
    handle = get_mtcr_handle()

    # Prepare output variables
    tt = Ref{Float64}(0.0)
    trot = Ref{Float64}(0.0)
    teex = Ref{Float64}(0.0)
    tvib = Ref{Float64}(0.0)

    # Prepare arrays for species-specific temperatures
    max_species = get_max_number_of_species_wrapper()
    tex = zeros(Float64, max_species)
    tvx = zeros(Float64, 10, max_species)  # Assuming max 10 electronic states

    # TODO: Uncomment when Fortran library is available (Phase 1 complete)
    # Call Fortran subroutine
    # Note: Simplified interface - actual implementation needs proper argument handling
    # ccall((:calculate_temperatures, lib_path), Cvoid,
    #     (Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    #         Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},),
    #     rho_sp, rho_etot, tt, trot, teex, tvib, tex, tvx,)

    # Placeholder: reasonable temperature values for testing
    tt[] = 300.0  # Translational temperature (K)
    trot[] = 300.0  # Rotational temperature (K)
    teex[] = 1000.0  # Electronic excitation temperature (K)
    tvib[] = 500.0  # Vibrational temperature (K)
    fill!(tex, 1000.0)  # Electronic temperatures
    fill!(tvx, 500.0)  # Vibrational temperatures

    return (tt = tt[], trot = trot[], teex = teex[], tvib = tvib[],
        tex = tex, tvx = tvx,)
end

"""
$(SIGNATURES)

Calculate total energy from state variables.
"""
function calculate_total_energy_wrapper(tt::Float64,
        rho_sp::Vector{Float64};
        kwargs...,)
    handle = get_mtcr_handle()

    rho_etot = Ref{Float64}(0.0)

    # TODO: Uncomment when Fortran library is available (Phase 1 complete)
    # Simplified call - actual implementation needs full argument handling
    # ccall((:calculate_total_energy, lib_path), Cvoid,
    #     (Ref{Float64}, Ptr{Float64}, Ref{Float64}),
    #     tt, rho_sp, rho_etot,)

    # Placeholder: calculate approximate total energy for testing
    rho_etot[] = sum(rho_sp) * tt * 1000.0  # Rough approximation
    return rho_etot[]
end

"""
$(SIGNATURES)

Get species names from MTCR.
"""
function get_species_names_wrapper()
    handle = get_mtcr_handle()

    max_species = get_max_number_of_species_wrapper()
    name_length = 32  # Assuming max name length

    # Allocate buffer for species names
    names_buffer = zeros(UInt8, name_length * max_species)

    # TODO: Uncomment when Fortran library is available (Phase 1 complete)
    # ccall((:get_species_names, lib_path), Cvoid,
    #     (Ptr{UInt8},), names_buffer,)

    # Placeholder: return typical nitrogen plasma species for testing
    species_names = ["N2", "N", "N+", "N2+", "e-"]

    return species_names
end
