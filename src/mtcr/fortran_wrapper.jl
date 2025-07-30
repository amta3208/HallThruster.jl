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
const MTCR_INITIALIZED = Ref{Bool}(false)

"""
$(SIGNATURES)

Set the path to the MTCR shared library and load it.

# Arguments
- `path::String`: Path to the MTCR shared library file

# Throws
- `ErrorException`: If the library cannot be loaded
"""
function load_mtcr_library!(path::String)
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
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
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
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
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
- `num_species::Integer`: Number of species in the system
- `num_dimensions::Integer`: Number of spatial dimensions (default: 0)
- `case_path::String`: Path to directory containing input/ subdirectory (default: current directory)

# Throws
- `ErrorException`: If case_path doesn't exist, input file is missing, or Fortran call fails
"""
function initialize_api_wrapper(
        num_species::Integer; num_dimensions::Integer = Int32(0), case_path::String = pwd(),)
    if MTCR_INITIALIZED[]
        @warn "MTCR already initialized in this Julia session - skipping"
        return nothing
    end

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
        MTCR_INITIALIZED[] = true
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
    if !MTCR_INITIALIZED[]
        @warn "MTCR not initialized - nothing to finalize"
        return nothing
    end

    ccall((:finalize_api, get_mtcr_lib_path()), Cvoid, ())

    MTCR_INITIALIZED[] = false
    return nothing
end

"""
$(SIGNATURES)

Get the maximum number of species supported by MTCR.

# Returns
- `Int32`: Maximum number of species

# Throws
- `ErrorException`: If MTCR library is not loaded
"""
function get_max_number_of_species_wrapper()
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    return ccall((:get_max_number_of_species, get_mtcr_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get species names from MTCR.

# Returns
- `Vector{String}`: Array of species names

# Throws
- `ErrorException`: If MTCR library is not loaded
"""
function get_species_names_wrapper()
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    # Get max species count first
    max_species = get_max_number_of_species_wrapper()

    # Fortran character length (from parameters module: nmlen)
    name_length = 32  # This should match nmlen in Fortran

    # Allocate buffer for species names (Fortran returns character array)
    names_buffer = zeros(UInt8, name_length * max_species)

    # Call Fortran subroutine
    ccall((:get_species_names, get_mtcr_lib_path()), Cvoid,
        (Ptr{UInt8},), names_buffer,)

    # Convert buffer to Julia strings
    species_names = String[]

    # The Fortran subroutine packs names sequentially with null terminators
    # Look for null-terminated strings
    i = 1
    while i <= length(names_buffer)
        # Find the next null terminator
        null_idx = findfirst(==(0), names_buffer[i:end])
        if null_idx === nothing
            break  # No more null terminators
        end

        # Extract the name bytes (excluding null terminator)
        name_end = i + null_idx - 2
        if name_end >= i
            name_bytes = names_buffer[i:name_end]
            name = String(name_bytes) |> strip
            if !isempty(name)
                push!(species_names, name)
            end
        end

        i += null_idx
        # Skip any additional null padding to get to next name
        while i <= length(names_buffer) && names_buffer[i] == 0
            i += 1
        end
    end

    return species_names
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

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
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
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    # Prepare output arrays
    drho_sp = similar(rho_sp)
    drho_etot = Ref{Float64}(0.0)

    # Handle optional output arguments
    drho_ex = rho_ex !== nothing ? similar(rho_ex) : nothing
    drho_vx = rho_vx !== nothing ? similar(rho_vx) : nothing
    drho_erot = rho_erot !== nothing ? Ref{Float64}(0.0) : nothing
    drho_eeex = rho_eeex !== nothing ? Ref{Float64}(0.0) : nothing
    drho_evib = rho_evib !== nothing ? Ref{Float64}(0.0) : nothing

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_nonequilibrium_sources, get_mtcr_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_sp
                Ptr{Float64},                                    # rho_ex (optional)
                Ptr{Float64},                                    # rho_vx (optional)
                Ptr{Float64},                                    # rho_u (optional)
                Ptr{Float64},                                    # rho_v (optional)
                Ptr{Float64},                                    # rho_w (optional)
                Ref{Float64},                                    # rho_etot
                Ptr{Float64},                                    # rho_erot (optional)
                Ptr{Float64},                                    # rho_eeex (optional)
                Ptr{Float64},                                    # rho_evib (optional)
                Ptr{Float64},                                    # drho_sp
                Ptr{Float64},                                    # drho_ex (optional)
                Ptr{Float64},                                    # drho_vx (optional)
                Ref{Float64},                                    # drho_etot
                Ptr{Float64},                                    # drho_erot (optional)
                Ptr{Float64},                                    # drho_eeex (optional)
                Ptr{Float64},),                                   # drho_evib (optional)
            rho_sp,
            rho_ex !== nothing ? rho_ex : C_NULL,
            rho_vx !== nothing ? rho_vx : C_NULL,
            rho_u !== nothing ? Ref{Float64}(rho_u) : C_NULL,
            rho_v !== nothing ? Ref{Float64}(rho_v) : C_NULL,
            rho_w !== nothing ? Ref{Float64}(rho_w) : C_NULL,
            rho_etot,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL,
            drho_sp,
            drho_ex !== nothing ? drho_ex : C_NULL,
            drho_vx !== nothing ? drho_vx : C_NULL,
            drho_etot,
            drho_erot !== nothing ? drho_erot : C_NULL,
            drho_eeex !== nothing ? drho_eeex : C_NULL,
            drho_evib !== nothing ? drho_evib : C_NULL,)
    catch e
        error("Failed to calculate source terms: $(e)")
    end

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

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
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
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    # Prepare output variables
    tt = Ref{Float64}(0.0)
    trot = Ref{Float64}(0.0)
    teex = Ref{Float64}(0.0)
    tvib = Ref{Float64}(0.0)

    # Prepare arrays for species-specific temperatures
    max_species = get_max_number_of_species_wrapper()
    tex = zeros(Float64, max_species)
    # Get max electronic states - using 10 as reasonable default for now
    max_electronic_states = 10  # TODO: Get from MTCR parameter function
    tvx = zeros(Float64, max_electronic_states, max_species)

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_temperatures, get_mtcr_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_sp
                Ptr{Float64},                                    # rho_ex (optional)
                Ptr{Float64},                                    # rho_vx (optional)
                Ptr{Float64},                                    # rho_u (optional)
                Ptr{Float64},                                    # rho_v (optional)
                Ptr{Float64},                                    # rho_w (optional)
                Ref{Float64},                                    # rho_etot
                Ptr{Float64},                                    # rho_erot (optional)
                Ptr{Float64},                                    # rho_eeex (optional)
                Ptr{Float64},                                    # rho_evib (optional)
                Ref{Float64},                                    # tt
                Ref{Float64},                                    # trot
                Ref{Float64},                                    # teex
                Ref{Float64},                                    # tvib
                Ptr{Float64},                                    # tex
                Ptr{Float64},),                                   # tvx
            rho_sp,
            rho_ex !== nothing ? rho_ex : C_NULL,
            rho_vx !== nothing ? rho_vx : C_NULL,
            rho_u !== nothing ? Ref{Float64}(rho_u) : C_NULL,
            rho_v !== nothing ? Ref{Float64}(rho_v) : C_NULL,
            rho_w !== nothing ? Ref{Float64}(rho_w) : C_NULL,
            rho_etot,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL,
            tt,
            trot,
            teex,
            tvib,
            tex,
            tvx,)
    catch e
        error("Failed to calculate temperatures: $(e)")
    end

    return (tt = tt[], trot = trot[], teex = teex[], tvib = tvib[],
        tex = tex, tvx = tvx,)
end

"""
$(SIGNATURES)

Calculate total energy from state variables.

# Arguments
- `tt::Float64`: Translational temperature (K)
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `rho_vx::Array{Float64,3}`: Vibrational state densities (optional)
- `u::Float64`: x-velocity component (optional)
- `v::Float64`: y-velocity component (optional)
- `w::Float64`: z-velocity component (optional)
- `rho_erot::Float64`: Rotational energy density (optional)
- `rho_eeex::Float64`: Electron-electronic energy density (optional)
- `rho_evib::Float64`: Vibrational energy density (optional)

# Returns
- `Float64`: Total energy density

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
"""
function calculate_total_energy_wrapper(tt::Float64,
        rho_sp::Vector{Float64};
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        u::Union{Float64, Nothing} = nothing,
        v::Union{Float64, Nothing} = nothing,
        w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing,)
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    rho_etot = Ref{Float64}(0.0)

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_total_energy, get_mtcr_lib_path()), Cvoid,
            (Ref{Float64},                                    # rho_etot (output)
                Ref{Float64},                                    # tt
                Ptr{Float64},                                    # rho_sp
                Ptr{Float64},                                    # rho_ex (optional)
                Ptr{Float64},                                    # rho_vx (optional)
                Ptr{Float64},                                    # u (optional)
                Ptr{Float64},                                    # v (optional)
                Ptr{Float64},                                    # w (optional)
                Ptr{Float64},                                    # rho_erot (optional)
                Ptr{Float64},                                    # rho_eeex (optional)
                Ptr{Float64},),                                   # rho_evib (optional)
            rho_etot,
            tt,
            rho_sp,
            rho_ex !== nothing ? rho_ex : C_NULL,
            rho_vx !== nothing ? rho_vx : C_NULL,
            u !== nothing ? Ref{Float64}(u) : C_NULL,
            v !== nothing ? Ref{Float64}(v) : C_NULL,
            w !== nothing ? Ref{Float64}(w) : C_NULL,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL,)
    catch e
        error("Failed to calculate total energy: $(e)")
    end

    return rho_etot[]
end
