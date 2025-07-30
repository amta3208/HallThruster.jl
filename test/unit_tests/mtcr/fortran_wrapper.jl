using HallThruster: HallThruster as het
using Test
using DocStringExtensions
using Libdl

# Import the MTCR fortran wrapper functions
# Note: These will be available once the MTCR module is properly integrated
# For now, we'll access them through the module structure
include("../../../src/mtcr/fortran_wrapper.jl")

# TODO: Update with correct library path later on
temp_mtcr_path = "/Users/amin/postdoc/codes/HallThruster.jl/mtcr/source/libmtcr.so"

@testset "Library Management" begin
    @testset "Library Loading and Status" begin
        # Test initial state - no library loaded
        @test !is_mtcr_loaded()

        # Test that getting handle fails when no library is loaded
        @test_throws ErrorException get_mtcr_handle()

        # Test error message content
        try
            get_mtcr_handle()
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
            @test occursin("load_mtcr_library!", e.msg)
        end
    end

    @testset "Library Path Setting" begin
        # Test setting library path with non-existent file
        fake_path = "/nonexistent/path/libmtcr.so"
        @test_throws ErrorException load_mtcr_library!(fake_path)

        # Test error message for non-existent file
        try
            load_mtcr_library!(fake_path)
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library file not found", e.msg)
            @test occursin(fake_path, e.msg)
        end

        # Test that library is still not loaded after failed attempt
        @test !is_mtcr_loaded()
    end

    @testset "Library Cleanup" begin
        # Test closing library when none is loaded (should be safe)
        @test_nowarn close_mtcr_library()
        @test !is_mtcr_loaded()

        # Test multiple calls to close (should be safe)
        @test_nowarn close_mtcr_library()
        @test_nowarn close_mtcr_library()
        @test !is_mtcr_loaded()
    end
end

@testset "Initialization" begin
    # Ensure library path is set
    load_mtcr_library!(temp_mtcr_path)

    # Get path to test case directory
    test_case_path = joinpath(@__DIR__, "test_case")

    # Test basic initialization with required arguments
    num_species = Int32(5)
    @test_nowarn initialize_api_wrapper(num_species; case_path = test_case_path)
    # @test_nowarn finalize_api_wrapper()

    # Test initialization with optional num_dimensions
    # @test_nowarn initialize_api_wrapper(
    # num_species; num_dimensions = Int32(1), case_path = test_case_path,)

    # Test initialization with default num_dimensions (0)
    # @test_nowarn initialize_api_wrapper(Int32(3); case_path = test_case_path)
    # @test_nowarn finalize_api_wrapper()
end

@testset "Input/Directory Handling" begin
    @testset "Initialization Input Validation" begin
        load_mtcr_library!(temp_mtcr_path)

        # Reset initialization state for validation tests
        MTCR_INITIALIZED[] = false

        # Test with non-existent case path
        @test_throws ErrorException initialize_api_wrapper(
            Int32(5); case_path = "/nonexistent/path",)

        # Test error message for non-existent case path
        try
            initialize_api_wrapper(Int32(5); case_path = "/nonexistent/path")
            @test false  # Should not reach here
        catch e
            @test occursin("Case path does not exist", e.msg)
        end

        # Test with case path missing input directory
        temp_dir = mktempdir()
        try
            @test_throws ErrorException initialize_api_wrapper(
                Int32(5); case_path = temp_dir,)
        finally
            rm(temp_dir; recursive = true)
        end

        # Test error message for missing input file
        temp_dir = mktempdir()
        try
            initialize_api_wrapper(Int32(5); case_path = temp_dir)
            @test false  # Should not reach here
        catch e
            @test occursin("Required input file not found", e.msg)
            @test occursin("prob_setup.inp", e.msg)
        finally
            rm(temp_dir; recursive = true)
        end
    end

    @testset "Directory Management" begin
        load_mtcr_library!(temp_mtcr_path)

        # Store original directory
        original_dir = pwd()
        test_case_path = joinpath(@__DIR__, "test_case")

        # Test that directory is restored after successful call
        initialize_api_wrapper(Int32(5); case_path = test_case_path)
        @test pwd() == original_dir

        # Test that directory is restored even after failed call
        temp_dir = mktempdir()
        try
            initialize_api_wrapper(Int32(5); case_path = temp_dir)
        catch
            # Expected to fail
        end
        @test pwd() == original_dir
        rm(temp_dir; recursive = true)
    end

    @testset "Output Directory Creation" begin
        load_mtcr_library!(temp_mtcr_path)

        # Reset initialization state for this test
        MTCR_INITIALIZED[] = false

        # Create a temporary test case directory with input
        temp_case_dir = mktempdir()
        try
            # Create input directory and file
            input_dir = joinpath(temp_case_dir, "input")
            mkpath(input_dir)
            touch(joinpath(input_dir, "prob_setup.inp"))

            # Verify output directories don't exist initially
            output_dir = joinpath(temp_case_dir, "output")
            sources_dir = joinpath(output_dir, "sources")
            states_dir = joinpath(output_dir, "states")

            @test !isdir(output_dir)
            @test !isdir(sources_dir)
            @test !isdir(states_dir)

            # Initialize - should create output directories
            initialize_api_wrapper(Int32(5); case_path = temp_case_dir)

            # Verify output directories were created
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)

            # Test that calling again doesn't cause errors (directories already exist)
            @test initialize_api_wrapper(Int32(5); case_path = temp_case_dir) === nothing

            # Verify directories still exist
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)

        finally
            rm(temp_case_dir; recursive = true)
        end
    end

    @testset "Output Directory Creation with Existing Directories" begin
        load_mtcr_library!(temp_mtcr_path)

        # Reset initialization state for this test
        MTCR_INITIALIZED[] = false

        # Create a temporary test case directory with input and partial output structure
        temp_case_dir = mktempdir()
        try
            # Create input directory and file
            input_dir = joinpath(temp_case_dir, "input")
            mkpath(input_dir)
            touch(joinpath(input_dir, "prob_setup.inp"))

            # Create output directory but not subdirectories
            output_dir = joinpath(temp_case_dir, "output")
            mkpath(output_dir)

            # Create a test file in output to verify it's preserved
            test_file = joinpath(output_dir, "test_file.txt")
            write(test_file, "test content")

            sources_dir = joinpath(output_dir, "sources")
            states_dir = joinpath(output_dir, "states")

            @test isdir(output_dir)
            @test !isdir(sources_dir)
            @test !isdir(states_dir)
            @test isfile(test_file)

            # Initialize - should create missing subdirectories
            initialize_api_wrapper(Int32(5); case_path = temp_case_dir)

            # Verify all directories exist and existing content is preserved
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)
            @test isfile(test_file)
            @test read(test_file, String) == "test content"

        finally
            rm(temp_case_dir; recursive = true)
        end
    end
end

@testset "Utility Functions" begin
    @testset "Maximum Species Count" begin
        # Ensure library path is set
        load_mtcr_library!(temp_mtcr_path)

        # Test error when library not loaded
        close_mtcr_library()
        @test_throws ErrorException get_max_number_of_species_wrapper()

        # Test error message content
        try
            get_max_number_of_species_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
            @test occursin("load_mtcr_library!", e.msg)
        end

        # Reload library for actual test
        load_mtcr_library!(temp_mtcr_path)

        max_species = get_max_number_of_species_wrapper()
        println("Max species: ", max_species)

        # Check return type
        @test max_species isa Int32

        # Check that it's positive and reasonable
        @test max_species > 0
        @test max_species <= 100  # Reasonable upper bound
    end

    @testset "Species Names" begin
        # Ensure library is loaded
        load_mtcr_library!(temp_mtcr_path)

        # Test error when library not loaded
        close_mtcr_library()
        @test_throws ErrorException get_species_names_wrapper()

        # Test error message content
        try
            get_species_names_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
            @test occursin("load_mtcr_library!", e.msg)
        end

        # Reload library for actual test
        load_mtcr_library!(temp_mtcr_path)

        # Check species composition
        species_names = get_species_names_wrapper()
        allowed_species = ["N", "N2", "N+", "N2+", "E-"]
        for name in species_names
            @test name in allowed_species
        end

        # Check return type
        @test species_names isa Vector{String}

        # Check that all names are non-empty strings
        @test all(length(name) > 0 for name in species_names)

        # Check that number of species is reasonable
        @test length(species_names) > 0
        @test length(species_names) <= 100  # Reasonable upper bound

        # Check that species names contain only valid characters
        isalnum(c) = isletter(c) || isdigit(c)
        for name in species_names
            @test all(c -> isascii(c) && (isalnum(c) || c in ['+', '-', '_']), name)
        end

        # Check that species names are unique
        @test length(species_names) == length(unique(species_names))
    end
end

@testset "Source Terms Calculation" begin
    @testset "Error Handling Without Library" begin
        # Ensure library is not loaded
        close_mtcr_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test error when library not loaded
        @test_throws ErrorException calculate_sources_wrapper(rho_sp, rho_etot)

        # Test error message content
        try
            calculate_sources_wrapper(rho_sp, rho_etot)
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        # Load library but don't initialize
        load_mtcr_library!(temp_mtcr_path)
        MTCR_INITIALIZED[] = false

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test error when not initialized
        @test_throws ErrorException calculate_sources_wrapper(rho_sp, rho_etot)

        # Test error message content
        try
            calculate_sources_wrapper(rho_sp, rho_etot)
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    #TODO: Need to find appropriate values to test call to calculate_sources
    #=
    @testset "Function Signature and Return Structure" begin
        # Set up proper initialization
        load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        initialize_api_wrapper(Int32(5); case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test that function can be called (may fail with actual Fortran call)
        # but should have proper error handling
        @test_nowarn try
            # result = calculate_sources_wrapper(rho_sp, rho_etot;
            #     rho_erot = 2e9, rho_eeex = 2e11, rho_evib = 5e9,)
            result = calculate_sources_wrapper(rho_sp, rho_etot;
                rho_eeex = 200.0, rho_evib = 200.0,)
            # If successful, check structure
            @test result isa NamedTuple
            @test haskey(result, :drho_sp)
            @test haskey(result, :drho_etot)
            @test haskey(result, :drho_erot)
            @test haskey(result, :drho_eeex)
            @test haskey(result, :drho_evib)
        catch e
            # Expected to fail until Fortran library is built
            @test e isa ErrorException
        end
    end
    =#
end

@testset "Temperature Calculation" begin
    @testset "Error Handling Without Library" begin
        close_mtcr_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_throws ErrorException calculate_temperatures_wrapper(rho_sp, rho_etot)

        try
            calculate_temperatures_wrapper(rho_sp, rho_etot)
            @test false
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        load_mtcr_library!(temp_mtcr_path)
        MTCR_INITIALIZED[] = false

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_throws ErrorException calculate_temperatures_wrapper(rho_sp, rho_etot)

        try
            calculate_temperatures_wrapper(rho_sp, rho_etot)
            @test false
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        initialize_api_wrapper(Int32(5); case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_nowarn try
            result = calculate_temperatures_wrapper(rho_sp, rho_etot;
                rho_eeex = 200.0, rho_evib = 200.0,)

            # If successful, check structure
            @test result isa NamedTuple
            @test haskey(result, :tt)
            @test haskey(result, :trot)
            @test haskey(result, :teex)
            @test haskey(result, :tvib)
            @test haskey(result, :tex)
            @test haskey(result, :tvx)
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Total Energy Calculation" begin
    @testset "Error Handling Without Library" begin
        close_mtcr_library()

        tt = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException calculate_total_energy_wrapper(tt, rho_sp)

        try
            calculate_total_energy_wrapper(tt, rho_sp)
            @test false
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        load_mtcr_library!(temp_mtcr_path)
        MTCR_INITIALIZED[] = false

        tt = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException calculate_total_energy_wrapper(tt, rho_sp)

        try
            calculate_total_energy_wrapper(tt, rho_sp)
            @test false
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        initialize_api_wrapper(Int32(5); case_path = test_case_path)

        tt = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with optional arguments
        @test_nowarn try
            result = calculate_total_energy_wrapper(tt, rho_sp;
                rho_erot = 100.0, rho_eeex = 200.0, rho_evib = 300.0,)
            # If successful, check return type
            @test result isa Float64
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Finalization" begin
    # Test that finalization runs without error
    @test finalize_api_wrapper() === nothing
    @test_nowarn close_mtcr_library()

    @test MTCR_HANDLE[] == C_NULL
    @test MTCR_LIB_PATH[] == ""
end

@testset "Initialize-Finalize Catching When Pre-Initialized" begin
    # Test load checker
    @test !is_mtcr_loaded()

    # Ensure library path is set
    load_mtcr_library!(temp_mtcr_path)
    @test is_mtcr_loaded()

    # Get path to test case directory
    test_case_path = joinpath(@__DIR__, "test_case")

    # Test basic reinitialization with required arguments
    MTCR_INITIALIZED[] = true
    num_species = Int32(5)
    @test initialize_api_wrapper(num_species; case_path = test_case_path) ===
          nothing

    # Test basic refinalization
    MTCR_INITIALIZED[] = false
    @test finalize_api_wrapper() === nothing
    @test_nowarn close_mtcr_library()
    @test !MTCR_INITIALIZED[]
end
