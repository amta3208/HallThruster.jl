"""
# MTCR Configuration Module

This module handles configuration management for MTCR simulations,
including parameter validation, default values, and input file generation.

## Key Types

- `MTCRConfig`: Main configuration struct
- `MTCRResults`: Results container
- `TemperatureConfig`: Temperature configuration
- `TimeIntegrationConfig`: Time integration parameters

## Key Functions

- `MTCRConfig`: Constructor with validation
- `validate_config`: Validate configuration parameters
- `generate_input_files`: Generate MTCR input files from configuration
"""

using DocStringExtensions

"""
Temperature configuration for MTCR simulation.

# Fields
- `Tt::Float64`: Translational temperature (K)
- `Tv::Float64`: Vibrational temperature (K)
- `Te::Float64`: Electron temperature (K)
- `Tee::Float64`: Electron-electronic temperature (K)
"""
struct TemperatureConfig
    Tt::Float64
    Tv::Float64
    Te::Float64
    Tee::Float64

    function TemperatureConfig(Tt, Tv, Te, Tee = Tt)
        if any([Tt, Tv, Te, Tee] .<= 0)
            error("All temperatures must be positive")
        end
        new(Float64(Tt), Float64(Tv), Float64(Te), Float64(Tee))
    end
end

"""
Time integration configuration for MTCR simulation.

# Fields
- `dt::Float64`: Time step (s)
- `dtm::Float64`: Output time step (s)
- `tlim::Float64`: Final time (s)
- `nstep::Int`: Maximum number of time steps
- `method::Int`: Integration method (0=forward Euler, 1=high order explicit, 2=implicit)
"""
struct TimeIntegrationConfig
    dt::Float64
    dtm::Float64
    tlim::Float64
    nstep::Int
    method::Int

    function TimeIntegrationConfig(dt, dtm, tlim, nstep = 500000, method = 2)
        if dt <= 0 || dtm <= 0 || tlim <= 0
            error("Time parameters must be positive")
        end
        if nstep <= 0
            error("Number of steps must be positive")
        end
        if !(method in [0, 1, 2])
            error("Integration method must be 0, 1, or 2")
        end
        new(Float64(dt), Float64(dtm), Float64(tlim), Int(nstep), Int(method))
    end
end

"""
Physics modeling configuration for MTCR simulation.

# Fields
- `bbh_model::Int`: Bound-bound heavy particle model
- `esc_model::Int`: Escape model
- `ar_et_model::Int`: Ar-ET model
- `eex_noneq::Int`: Electron-electronic nonequilibrium flag
- `ev_relax_set::Int`: Electron-vibrational relaxation set
- `et_relax_set::Int`: Electron-translational relaxation set
"""
struct PhysicsConfig
    bbh_model::Int
    esc_model::Int
    ar_et_model::Int
    eex_noneq::Int
    ev_relax_set::Int
    et_relax_set::Int

    function PhysicsConfig(;
            bbh_model = 4,
            esc_model = 1,
            ar_et_model = 1,
            eex_noneq = 1,
            ev_relax_set = 1,
            et_relax_set = 1,
    )
        new(bbh_model, esc_model, ar_et_model, eex_noneq, ev_relax_set, et_relax_set)
    end
end

"""
Process flags configuration for MTCR simulation.

# Fields
- `consider_elec_bbe::Int`: Consider electron bound-bound excitation
- `consider_elec_bfe::Int`: Consider electron bound-free excitation
- `consider_elec_bbh::Int`: Consider electron bound-bound heavy
- `consider_elec_bfh::Int`: Consider electron bound-free heavy
- `consider_rad::Int`: Consider radiation
- `consider_rdr::Int`: Consider RDR
- `consider_chem::Int`: Consider chemistry
"""
struct ProcessConfig
    consider_elec_bbe::Int
    consider_elec_bfe::Int
    consider_elec_bbh::Int
    consider_elec_bfh::Int
    consider_rad::Int
    consider_rdr::Int
    consider_chem::Int

    function ProcessConfig(;
            consider_elec_bbe = 1,
            consider_elec_bfe = 1,
            consider_elec_bbh = 1,
            consider_elec_bfh = 1,
            consider_rad = 0,
            consider_rdr = 0,
            consider_chem = 1,
    )
        new(consider_elec_bbe, consider_elec_bfe, consider_elec_bbh,
            consider_elec_bfh, consider_rad, consider_rdr, consider_chem,)
    end
end

"""
Main configuration struct for MTCR simulations.

# Fields
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Initial mole fractions
- `total_number_density::Float64`: Total number density (1/cm³)
- `temperatures::TemperatureConfig`: Temperature configuration
- `time_params::TimeIntegrationConfig`: Time integration parameters
- `physics::PhysicsConfig`: Physics modeling options
- `processes::ProcessConfig`: Process flags
- `database_path::String`: Path to chemistry database
- `radiation_length::Float64`: Radiation length scale (cm)
- `print_source_terms::Bool`: Print source terms flag
- `get_electron_density_by_charge_balance::Bool`: Electron density by charge balance
- `min_sts_frac::Float64`: Minimum state-to-state fraction
- `is_isothermal_teex::Bool`: Isothermal electron-electronic flag
"""
struct MTCRConfig
    species::Vector{String}
    mole_fractions::Vector{Float64}
    total_number_density::Float64
    temperatures::TemperatureConfig
    time_params::TimeIntegrationConfig
    physics::PhysicsConfig
    processes::ProcessConfig
    database_path::String
    radiation_length::Float64
    print_source_terms::Bool
    get_electron_density_by_charge_balance::Bool
    min_sts_frac::Float64
    is_isothermal_teex::Bool

    function MTCRConfig(;
            species::Vector{String},
            mole_fractions::Vector{Float64},
            total_number_density::Float64,
            temperatures::TemperatureConfig,
            time_params::TimeIntegrationConfig,
            physics::PhysicsConfig = PhysicsConfig(),
            processes::ProcessConfig = ProcessConfig(),
            database_path::String = "../../databases/n2/elec_sts_expanded_electron_fits_ground",
            radiation_length::Float64 = 1.0,
            print_source_terms::Bool = true,
            get_electron_density_by_charge_balance::Bool = true,
            min_sts_frac::Float64 = 1e-30,
            is_isothermal_teex::Bool = true,
    )

        # Validate inputs
        validate_config(
            species, mole_fractions, total_number_density, temperatures, time_params,)

        new(species, mole_fractions, total_number_density, temperatures, time_params,
            physics, processes, database_path, radiation_length, print_source_terms,
            get_electron_density_by_charge_balance, min_sts_frac, is_isothermal_teex,)
    end
end

"""
Results container for MTCR simulations.

# Fields
- `time::Vector{Float64}`: Time points
- `species_densities::Matrix{Float64}`: Species densities over time
- `temperatures::NamedTuple`: Temperature evolution
- `total_energy::Vector{Float64}`: Total energy evolution
- `source_terms::Union{NamedTuple, Nothing}`: Source terms (if requested)
- `success::Bool`: Simulation success flag
- `message::String`: Status message
"""
struct MTCRResults
    time::Vector{Float64}
    species_densities::Matrix{Float64}
    temperatures::NamedTuple
    total_energy::Vector{Float64}
    source_terms::Union{NamedTuple, Nothing}
    success::Bool
    message::String
end

"""
$(SIGNATURES)

Validate MTCR configuration parameters.

# Arguments
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Initial mole fractions
- `total_number_density::Float64`: Total number density
- `temperatures::TemperatureConfig`: Temperature configuration
- `time_params::TimeIntegrationConfig`: Time integration parameters

# Throws
- `ArgumentError` if validation fails
"""
function validate_config(species::Vector{String},
        mole_fractions::Vector{Float64},
        total_number_density::Float64,
        temperatures::TemperatureConfig,
        time_params::TimeIntegrationConfig,)

    # Check species and mole fractions consistency
    if length(species) != length(mole_fractions)
        throw(ArgumentError("Species and mole_fractions arrays must have same length"))
    end

    if isempty(species)
        throw(ArgumentError("At least one species must be specified"))
    end

    # Check mole fractions sum to 1
    if abs(sum(mole_fractions) - 1.0) > 1e-10
        throw(ArgumentError("Mole fractions must sum to 1.0, got $(sum(mole_fractions))"))
    end

    # Check for negative mole fractions
    if any(mole_fractions .< 0)
        throw(ArgumentError("Mole fractions must be non-negative"))
    end

    # Check total number density
    if total_number_density <= 0
        throw(ArgumentError("Total number density must be positive"))
    end

    # Check for duplicate species
    if length(unique(species)) != length(species)
        throw(ArgumentError("Duplicate species names found"))
    end

    # Validate species names format
    for species_name in species
        if isempty(strip(species_name))
            throw(ArgumentError("Species names cannot be empty"))
        end
    end

    return true
end

"""
$(SIGNATURES)

Generate MTCR input files from configuration.

# Arguments
- `config::MTCRConfig`: MTCR configuration
- `output_dir::String`: Output directory for input files

# Returns
- `true` if files generated successfully
"""
function generate_input_files(config::MTCRConfig, output_dir::String)

    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Generate prob_setup.inp
    generate_prob_setup_file(config, joinpath(output_dir, "prob_setup.inp"))

    # Generate sources_setup.inp
    generate_sources_setup_file(config, joinpath(output_dir, "sources_setup.inp"))

    # Generate tau_scaling.inp (empty for now)
    generate_tau_scaling_file(config, joinpath(output_dir, "tau_scaling.inp"))

    return true
end

"""
$(SIGNATURES)

Generate prob_setup.inp file from configuration.
"""
function generate_prob_setup_file(config::MTCRConfig, filepath::String)
    open(filepath, "w") do io
        println(io, "####################################################")
        println(io, "# Location of database and output folders")
        println(io, "####################################################")
        println(io, "DATABASE_PATH=$(config.database_path)")
        println(io, "CHEM_FILE_NAME=chemistry.dat")
        println(io)
        println(io, "--- Turn on source term printouts")
        println(io, "PRINT_SOURCE_TERMS=$(config.print_source_terms ? 1 : 0)")
        println(io,
            "GET_ELECTRON_DENSITY_BY_CHARGE_BALANCE=$(config.get_electron_density_by_charge_balance ? 1 : 0)",)
        println(io, "MIN_STS_FRAC=$(config.min_sts_frac)")
        println(io)
        println(io, "####################################################")
        println(io, "# Freestream condition")
        println(io, "####################################################")
        println(io)
        println(io, "---  Number of species, must match chemistry.dat")
        println(io, "NSP=$(length(config.species))")
        println(io)
        println(io, "--- Species mole fractions ($(join(config.species, ", ")))")
        for (i, frac) in enumerate(config.mole_fractions)
            println(io, "X$(i)=$(frac)")
        end
        println(io)
        println(io, "--- Total number density (1/cm³)")
        println(io, "TOTAL_NUMBER_DENSITY=$(config.total_number_density)")
        println(io)
        println(io, "--- Temperatures (K)")
        println(io, "TT=$(config.temperatures.Tt)")
        println(io, "TV=$(config.temperatures.Tv)")
        println(io, "TE=$(config.temperatures.Te)")
        println(io, "TEE=$(config.temperatures.Tee)")
        println(io)
        println(io, "--- Radiation length scale (cm)")
        println(io, "RAD_LEN=$(config.radiation_length)")
        println(io)
        println(io, "####################################################")
        println(io, "# Physical modeling variables")
        println(io, "####################################################")
        println(io, "--- Physics options")
        println(io, "BBHMODEL=$(config.physics.bbh_model)")
        println(io, "ESC_MODEL=$(config.physics.esc_model)")
        println(io, "AR_ET_MODEL=$(config.physics.ar_et_model)")
        println(io, "EEX_NONEQ=$(config.physics.eex_noneq)")
        println(io, "EV_RELAX_SET=$(config.physics.ev_relax_set)")
        println(io, "ET_RELAX_SET=$(config.physics.et_relax_set)")
        println(io)
        println(io, "--- Process flags")
        println(io, "CONSIDER_ELEC_BBE=$(config.processes.consider_elec_bbe)")
        println(io, "CONSIDER_ELEC_BFE=$(config.processes.consider_elec_bfe)")
        println(io, "CONSIDER_ELEC_BBH=$(config.processes.consider_elec_bbh)")
        println(io, "CONSIDER_ELEC_BFH=$(config.processes.consider_elec_bfh)")
        println(io, "CONSIDER_RAD=$(config.processes.consider_rad)")
        println(io, "CONSIDER_RDR=$(config.processes.consider_rdr)")
        println(io, "CONSIDER_CHEM=$(config.processes.consider_chem)")
        println(io)
        println(io, "####################################################")
        println(io, "# Computational setup")
        println(io, "####################################################")
        println(io,
            "--- Time integration method: forward euler == 0, high order explicit == 1, numerical implicit == 2",)
        println(io, "TIME_METHOD=$(config.time_params.method)")
        println(io)
        println(io, "IS_ISOTHERMAL_TEEX=$(config.is_isothermal_teex ? 1 : 0)")
        println(io)
        println(io, "--- Number of dimensions, 0 or 1")
        println(io, "ND=0")
        println(io)
        println(io, "--- 0D time integration setup (time in microseconds)")
        println(io, "DT=$(config.time_params.dt)")
        println(io, "DTM=$(config.time_params.dtm)")
        println(io, "TLIM=$(config.time_params.tlim)")
        println(io)
        println(io, "-- Max number of iterations")
        println(io, "NSTEP=$(config.time_params.nstep)")
    end
end

"""
$(SIGNATURES)

Generate sources_setup.inp file from configuration.
"""
function generate_sources_setup_file(config::MTCRConfig, filepath::String)
    open(filepath, "w") do io
        println(io, "BEGIN SPECIES SOURCES")
        for species in config.species
            println(io, species)
        end
        println(io, "END SPECIES SOURCES")
        println(io)
        println(io, "BEGIN EXCITED STATE SOURCES")
        # For now, include common excited states for nitrogen
        # This should be made more general based on the species
        if "N" in config.species
            for i in 1:5
                println(io, "N($(i))")
            end
        end
        if "N2" in config.species
            for i in 1:9
                println(io, "N2($(i))")
            end
        end
        if "N2+" in config.species
            for i in 1:4
                println(io, "N2+($(i))")
            end
        end
        println(io, "END EXCITED STATE SOURCES")
    end
end

"""
$(SIGNATURES)

Generate tau_scaling.inp file from configuration.
"""
function generate_tau_scaling_file(config::MTCRConfig, filepath::String)
    # For now, create an empty file
    # This can be expanded later if tau scaling is needed
    open(filepath, "w") do io
        println(io, "# Tau scaling file - currently empty")
    end
end

"""
$(SIGNATURES)

Create a default configuration for the 0D Nitrogen Te=10eV example.

# Returns
- `MTCRConfig`: Configuration matching the example case
"""
function nitrogen_10ev_config()
    species = ["N", "N2", "N+", "N2+", "E-"]
    mole_fractions = [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001]
    total_number_density = 1.0e13  # 1/cm³

    temperatures = TemperatureConfig(750.0, 750.0, 115000.0, 750.0)
    time_params = TimeIntegrationConfig(0.5e-5, 5.0, 1e3, 500000, 2)

    return MTCRConfig(
        species = species,
        mole_fractions = mole_fractions,
        total_number_density = total_number_density,
        temperatures = temperatures,
        time_params = time_params,
    )
end

"""
$(SIGNATURES)

Get molecular weights for common species (g/mol).

# Arguments
- `species::Vector{String}`: Species names

# Returns
- `Vector{Float64}`: Molecular weights in g/mol
"""
function get_molecular_weights(species::Vector{String})
    # Common molecular weights (g/mol)
    molecular_weight_db = Dict(
        "N" => 14.007,
        "N2" => 28.014,
        "N+" => 14.007,
        "N2+" => 28.014,
        "E-" => 5.485799e-4,  # electron mass
        "Ar" => 39.948,
        "Ar+" => 39.948,
        "Xe" => 131.293,
        "Xe+" => 131.293,
        "Kr" => 83.798,
        "Kr+" => 83.798,
        "O" => 15.999,
        "O2" => 31.998,
        "O+" => 15.999,
        "O2+" => 31.998,
    )

    weights = Float64[]
    for sp in species
        if haskey(molecular_weight_db, sp)
            push!(weights, molecular_weight_db[sp])
        else
            error("Molecular weight not found for species: $(sp)")
        end
    end

    return weights
end
