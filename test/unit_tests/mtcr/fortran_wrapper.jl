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
            @test occursin("set_mtcr_lib_path!", e.msg)
        end
    end

    @testset "Library Path Setting" begin
        # Test setting library path with non-existent file
        fake_path = "/nonexistent/path/libmtcr.so"
        @test_throws ErrorException set_mtcr_lib_path!(fake_path)

        # Test error message for non-existent file
        try
            set_mtcr_lib_path!(fake_path)
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

@testset "API Lifecycle Management" begin
    @testset "Initialization" begin
        # Ensure library path is set
        set_mtcr_lib_path!(temp_mtcr_path)

        # Get path to test case directory
        test_case_path = joinpath(@__DIR__, "test_case")

        # Test basic initialization with required arguments
        num_species = Int32(5)
        @test_nowarn initialize_api_wrapper(num_species; case_path = test_case_path)

        # Test initialization with optional num_dimensions
        # @test_nowarn initialize_api_wrapper(
        # Int32(5); num_dimensions = Int32(1), case_path = test_case_path,)

        # Test initialization with default num_dimensions (0)
        # @test_nowarn initialize_api_wrapper(Int32(3); case_path = test_case_path)
    end

    @testset "Initialization Input Validation" begin
        set_mtcr_lib_path!(temp_mtcr_path)

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
        set_mtcr_lib_path!(temp_mtcr_path)

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
        set_mtcr_lib_path!(temp_mtcr_path)

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
            @test_nowarn initialize_api_wrapper(Int32(5); case_path = temp_case_dir)

            # Verify directories still exist
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)

        finally
            rm(temp_case_dir; recursive = true)
        end
    end

    @testset "Output Directory Creation with Existing Directories" begin
        set_mtcr_lib_path!(temp_mtcr_path)

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

    @testset "Finalization" begin
        # Test that finalization runs without error
        @test_nowarn finalize_api_wrapper()

        # Test multiple calls to finalization (should be safe)
        @test_nowarn finalize_api_wrapper()
        @test_nowarn finalize_api_wrapper()
    end

    @testset "Initialization-Finalization Sequence" begin
        set_mtcr_lib_path!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")

        # Test proper sequence
        @test_nowarn begin
            initialize_api_wrapper(Int32(5); case_path = test_case_path)
            finalize_api_wrapper()
        end

        # Test multiple initialization-finalization cycles
        for i in 1:3
            @test_nowarn initialize_api_wrapper(Int32(5); case_path = test_case_path)
            @test_nowarn finalize_api_wrapper()
        end
    end
end

@testset "Utility Functions" begin
    @testset "Maximum Species Count" begin
        # Ensure library path is set
        set_mtcr_lib_path!(temp_mtcr_path)

        max_species = get_max_number_of_species_wrapper()

        # Check return type
        @test max_species isa Int32

        # Check placeholder value
        @test max_species == 10

        # Check that it's positive and reasonable
        @test max_species > 0
        @test max_species <= 100  # Reasonable upper bound
    end

    @testset "Species Names" begin
        species_names = get_species_names_wrapper()

        # Check return type
        @test species_names isa Vector{String}

        # Check placeholder values (nitrogen plasma species)
        expected_species = ["N2", "N", "N+", "N2+", "e-"]
        @test species_names == expected_species

        # Check that all names are non-empty strings
        @test all(length(name) > 0 for name in species_names)

        # Check that species names are reasonable
        @test "N2" in species_names
        @test "N" in species_names
        @test "e-" in species_names

        # Check number of species matches expected
        @test length(species_names) == 5
    end
end

@testset "Source Terms Calculation" begin
    @testset "Basic Functionality" begin
        # Set up test data based on nitrogen plasma example
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]  # g/cm³ (N2, N, N+, N2+, e-)
        rho_etot = 1e4  # erg/cm³

        # Test basic call
        result = calculate_sources_wrapper(rho_sp, rho_etot)

        # Check return type and structure
        @test result isa NamedTuple
        @test haskey(result, :drho_sp)
        @test haskey(result, :drho_etot)
        @test haskey(result, :drho_ex)
        @test haskey(result, :drho_vx)
        @test haskey(result, :drho_erot)
        @test haskey(result, :drho_eeex)
        @test haskey(result, :drho_evib)

        # Check output types and sizes
        @test result.drho_sp isa Vector{Float64}
        @test length(result.drho_sp) == length(rho_sp)
        @test result.drho_etot isa Float64

        # Check placeholder values (should be zero)
        @test all(result.drho_sp .== 0.0)
        @test result.drho_etot == 0.0

        # Check optional outputs are nothing when not provided
        @test result.drho_ex === nothing
        @test result.drho_vx === nothing
        @test result.drho_erot === nothing
        @test result.drho_eeex === nothing
        @test result.drho_evib === nothing
    end

    @testset "With Optional Arguments" begin
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test with electronic states
        rho_ex = [1e-5 2e-5; 3e-5 4e-5; 5e-5 6e-5; 7e-5 8e-5; 9e-5 1e-4]  # 5 species, 2 electronic states
        result_ex = calculate_sources_wrapper(rho_sp, rho_etot; rho_ex = rho_ex)

        @test result_ex.drho_ex isa Matrix{Float64}
        @test size(result_ex.drho_ex) == size(rho_ex)

        # Test with vibrational states
        rho_vx = zeros(Float64, 5, 2, 3)  # 5 species, 2 electronic states, 3 vibrational levels
        result_vx = calculate_sources_wrapper(rho_sp, rho_etot; rho_vx = rho_vx)

        @test result_vx.drho_vx isa Array{Float64, 3}
        @test size(result_vx.drho_vx) == size(rho_vx)

        # Test with energy components
        result_energy = calculate_sources_wrapper(rho_sp, rho_etot;
            rho_erot = 100.0,
            rho_eeex = 200.0,
            rho_evib = 300.0,)

        @test result_energy.drho_erot isa Float64
        @test result_energy.drho_eeex isa Float64
        @test result_energy.drho_evib isa Float64
        @test result_energy.drho_erot == 0.0  # Placeholder
        @test result_energy.drho_eeex == 0.0  # Placeholder
        @test result_energy.drho_evib == 0.0  # Placeholder

        # Test with velocity components
        result_velocity = calculate_sources_wrapper(rho_sp, rho_etot;
            rho_u = 1000.0,
            rho_v = 2000.0,
            rho_w = 3000.0,)

        # Should not affect output structure (velocities are input only)
        @test haskey(result_velocity, :drho_sp)
        @test haskey(result_velocity, :drho_etot)
    end

    @testset "Input Validation and Edge Cases" begin
        # Test with empty species array
        @test_nowarn calculate_sources_wrapper(Float64[], 0.0)

        # Test with single species
        result_single = calculate_sources_wrapper([1e-3], 1e4)
        @test length(result_single.drho_sp) == 1
        @test result_single.drho_sp[1] == 0.0

        # Test with zero densities
        zero_rho = zeros(Float64, 5)
        result_zero = calculate_sources_wrapper(zero_rho, 0.0)
        @test all(result_zero.drho_sp .== 0.0)
        @test result_zero.drho_etot == 0.0

        # Test with very small values
        small_rho = [1e-20, 1e-25, 1e-30, 1e-35, 1e-40]
        result_small = calculate_sources_wrapper(small_rho, 1e-10)
        @test length(result_small.drho_sp) == 5
        @test all(isfinite.(result_small.drho_sp))

        # Test with large values
        large_rho = [1e10, 1e8, 1e6, 1e4, 1e2]
        result_large = calculate_sources_wrapper(large_rho, 1e15)
        @test length(result_large.drho_sp) == 5
        @test all(isfinite.(result_large.drho_sp))

        # Test with mixed positive/negative values (negative densities shouldn't occur physically)
        mixed_rho = [1e-3, -1e-6, 1e-7, -1e-8, 1e-9]
        @test_nowarn calculate_sources_wrapper(mixed_rho, 1e4)
    end

    @testset "Memory Management" begin
        # Test that output arrays are properly allocated
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        result = calculate_sources_wrapper(rho_sp, rho_etot)

        # Check that output array is different object from input
        @test result.drho_sp !== rho_sp

        # Check that modifying input doesn't affect output
        original_output = copy(result.drho_sp)
        rho_sp[1] = 999.0
        @test result.drho_sp == original_output

        # Test with large arrays to check memory allocation
        large_rho = zeros(Float64, 1000)
        large_result = calculate_sources_wrapper(large_rho, 1e4)
        @test length(large_result.drho_sp) == 1000
        @test all(large_result.drho_sp .== 0.0)
    end
end

@testset "Temperature Calculation" begin
    @testset "Basic Functionality" begin
        # Set up test data
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]  # g/cm³
        rho_etot = 1e4  # erg/cm³

        result = calculate_temperatures_wrapper(rho_sp, rho_etot)

        # Check return type and structure
        @test result isa NamedTuple
        @test haskey(result, :tt)
        @test haskey(result, :trot)
        @test haskey(result, :teex)
        @test haskey(result, :tvib)
        @test haskey(result, :tex)
        @test haskey(result, :tvx)

        # Check output types
        @test result.tt isa Float64
        @test result.trot isa Float64
        @test result.teex isa Float64
        @test result.tvib isa Float64
        @test result.tex isa Vector{Float64}
        @test result.tvx isa Matrix{Float64}

        # Check placeholder values (from current implementation)
        @test result.tt == 300.0  # Translational temperature
        @test result.trot == 300.0  # Rotational temperature
        @test result.teex == 1000.0  # Electronic excitation temperature
        @test result.tvib == 500.0  # Vibrational temperature

        # Check array sizes
        max_species = get_max_number_of_species_wrapper()
        @test length(result.tex) == max_species
        @test size(result.tvx) == (10, max_species)  # 10 electronic states

        # Check that all temperatures are positive
        @test result.tt > 0
        @test result.trot > 0
        @test result.teex > 0
        @test result.tvib > 0
        @test all(result.tex .> 0)
        @test all(result.tvx .> 0)

        # Check that temperatures are physically reasonable
        @test result.tt < 1e6  # Not too hot
        @test result.trot < 1e6
        @test result.teex < 1e6
        @test result.tvib < 1e6
    end

    @testset "With Optional Arguments" begin
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test with electronic states
        rho_ex = [1e-5 2e-5; 3e-5 4e-5; 5e-5 6e-5; 7e-5 8e-5; 9e-5 1e-4]
        result_ex = calculate_temperatures_wrapper(rho_sp, rho_etot; rho_ex = rho_ex)

        # Should still return same structure
        @test haskey(result_ex, :tt)
        @test haskey(result_ex, :tex)
        @test haskey(result_ex, :tvx)

        # Test with vibrational states
        rho_vx = zeros(Float64, 5, 2, 3)
        result_vx = calculate_temperatures_wrapper(rho_sp, rho_etot; rho_vx = rho_vx)

        @test haskey(result_vx, :tvib)
        @test haskey(result_vx, :tvx)

        # Test with energy components
        result_energy = calculate_temperatures_wrapper(rho_sp, rho_etot;
            rho_erot = 100.0,
            rho_eeex = 200.0,
            rho_evib = 300.0,)

        @test haskey(result_energy, :trot)
        @test haskey(result_energy, :teex)
        @test haskey(result_energy, :tvib)

        # Test with velocity components
        result_velocity = calculate_temperatures_wrapper(rho_sp, rho_etot;
            rho_u = 1000.0,
            rho_v = 2000.0,
            rho_w = 3000.0,)

        # Should not affect temperature calculation
        @test result_velocity.tt == 300.0
    end

    @testset "Input Validation and Edge Cases" begin
        # Test with empty species array
        result_empty = calculate_temperatures_wrapper(Float64[], 0.0)
        @test result_empty isa NamedTuple
        @test all(isfinite.([
            result_empty.tt, result_empty.trot, result_empty.teex, result_empty.tvib,]))

        # Test with single species
        result_single = calculate_temperatures_wrapper([1e-3], 1e4)
        @test result_single.tt == 300.0
        @test result_single.trot == 300.0

        # Test with zero energy
        result_zero = calculate_temperatures_wrapper([1e-3, 1e-6], 0.0)
        @test all(isfinite.([
            result_zero.tt, result_zero.trot, result_zero.teex, result_zero.tvib,]))

        # Test with very small densities
        small_rho = [1e-30, 1e-35, 1e-40]
        result_small = calculate_temperatures_wrapper(small_rho, 1e-10)
        @test all(isfinite.([
            result_small.tt, result_small.trot, result_small.teex, result_small.tvib,]))

        # Test with large values
        large_rho = [1e10, 1e8, 1e6]
        result_large = calculate_temperatures_wrapper(large_rho, 1e15)
        @test all(isfinite.([
            result_large.tt, result_large.trot, result_large.teex, result_large.tvib,]))
    end

    @testset "Temperature Consistency" begin
        # Test that temperatures are consistent across calls
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        result1 = calculate_temperatures_wrapper(rho_sp, rho_etot)
        result2 = calculate_temperatures_wrapper(rho_sp, rho_etot)

        @test result1.tt == result2.tt
        @test result1.trot == result2.trot
        @test result1.teex == result2.teex
        @test result1.tvib == result2.tvib
        @test result1.tex == result2.tex
        @test result1.tvx == result2.tvx

        # Test that different inputs give potentially different results
        # (though with placeholder implementation, they'll be the same)
        different_rho = [2e-3, 2e-6, 2e-7, 2e-7, 2e-10]
        result_diff = calculate_temperatures_wrapper(different_rho, 2e4)

        # With placeholder implementation, should be same
        @test result_diff.tt == result1.tt
        @test result_diff.trot == result1.trot
    end
end

@testset "Total Energy Calculation" begin
    @testset "Basic Functionality" begin
        # Set up test data
        tt = 300.0  # K
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]  # g/cm³

        result = calculate_total_energy_wrapper(tt, rho_sp)

        # Check return type
        @test result isa Float64

        # Check that result is positive (energy should be positive)
        @test result > 0

        # Check that result is finite
        @test isfinite(result)

        # Check placeholder calculation (rough approximation)
        expected_approx = sum(rho_sp) * tt * 1000.0
        @test result ≈ expected_approx
    end

    @testset "With Different Temperatures" begin
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with different temperatures
        temperatures = [100.0, 300.0, 1000.0, 5000.0, 10000.0]

        for temp in temperatures
            result = calculate_total_energy_wrapper(temp, rho_sp)
            @test result > 0
            @test isfinite(result)

            # Energy should scale with temperature
            expected = sum(rho_sp) * temp * 1000.0
            @test result ≈ expected
        end

        # Test that higher temperature gives higher energy
        result_low = calculate_total_energy_wrapper(300.0, rho_sp)
        result_high = calculate_total_energy_wrapper(3000.0, rho_sp)
        @test result_high > result_low
    end

    @testset "With Different Densities" begin
        tt = 1000.0

        # Test with different density profiles
        density_sets = [
            [1e-3, 1e-6, 1e-7, 1e-7, 1e-10],
            [2e-3, 2e-6, 2e-7, 2e-7, 2e-10],
            [1e-4, 1e-7, 1e-8, 1e-8, 1e-11],
            [1e-2, 1e-5, 1e-6, 1e-6, 1e-9],
        ]

        for rho_sp in density_sets
            result = calculate_total_energy_wrapper(tt, rho_sp)
            @test result > 0
            @test isfinite(result)

            # Check placeholder calculation
            expected = sum(rho_sp) * tt * 1000.0
            @test result ≈ expected
        end

        # Test that higher total density gives higher energy
        low_density = [1e-4, 1e-7, 1e-8, 1e-8, 1e-11]
        high_density = [1e-2, 1e-5, 1e-6, 1e-6, 1e-9]

        result_low = calculate_total_energy_wrapper(tt, low_density)
        result_high = calculate_total_energy_wrapper(tt, high_density)
        @test result_high > result_low
    end

    @testset "Edge Cases" begin
        # Test with zero temperature
        result_zero_temp = calculate_total_energy_wrapper(0.0, [1e-3, 1e-6])
        @test result_zero_temp == 0.0

        # Test with zero densities
        result_zero_density = calculate_total_energy_wrapper(1000.0, [0.0, 0.0])
        @test result_zero_density == 0.0

        # Test with single species
        result_single = calculate_total_energy_wrapper(500.0, [1e-3])
        @test result_single > 0
        @test result_single ≈ 1e-3 * 500.0 * 1000.0

        # Test with empty species array
        result_empty = calculate_total_energy_wrapper(1000.0, Float64[])
        @test result_empty == 0.0

        # Test with very small values
        result_small = calculate_total_energy_wrapper(1e-10, [1e-30])
        @test result_small ≥ 0
        @test isfinite(result_small)

        # Test with large values
        result_large = calculate_total_energy_wrapper(1e6, [1e6])
        @test result_large > 0
        @test isfinite(result_large)
    end

    @testset "With Optional Arguments" begin
        tt = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7]

        # Test with various optional arguments (should not affect result in current implementation)
        result_basic = calculate_total_energy_wrapper(tt, rho_sp)

        result_with_kwargs = calculate_total_energy_wrapper(tt, rho_sp;
            rho_ex = [1e-5 2e-5; 3e-5 4e-5; 5e-5 6e-5],
            rho_vx = zeros(Float64, 3, 2, 3),
            rho_erot = 100.0,
            rho_eeex = 200.0,
            rho_evib = 300.0,)

        # With current placeholder implementation, should be same
        @test result_with_kwargs == result_basic
    end
end

@testset "Integration Readiness" begin
    @testset "Function Signatures" begin
        # Test that all functions have expected signatures and can be called

        # Library management functions
        @test_nowarn set_mtcr_lib_path!("/test/path")
        @test_nowarn is_mtcr_loaded()
        @test_nowarn close_mtcr_library()

        # API lifecycle (with proper arguments)
        test_case_path = joinpath(@__DIR__, "test_case")
        @test_nowarn initialize_api_wrapper(Int32(5); case_path = test_case_path)
        @test_nowarn finalize_api_wrapper()

        # Utility functions
        @test_nowarn get_max_number_of_species_wrapper()
        @test_nowarn get_species_names_wrapper()

        # Core computation functions
        test_rho_sp = [1e-3, 1e-6, 1e-7]
        test_rho_etot = 1e4
        test_tt = 1000.0

        @test_nowarn calculate_sources_wrapper(test_rho_sp, test_rho_etot)
        @test_nowarn calculate_temperatures_wrapper(test_rho_sp, test_rho_etot)
        @test_nowarn calculate_total_energy_wrapper(test_tt, test_rho_sp)
    end

    @testset "Return Value Consistency" begin
        # Test that return values have consistent structure for interface compatibility

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Sources calculation
        sources = calculate_sources_wrapper(rho_sp, rho_etot)
        @test sources isa NamedTuple
        @test length(sources) == 7  # All expected fields

        # Temperature calculation
        temps = calculate_temperatures_wrapper(rho_sp, rho_etot)
        @test temps isa NamedTuple
        @test length(temps) == 6  # All expected temperature fields

        # Energy calculation
        energy = calculate_total_energy_wrapper(1000.0, rho_sp)
        @test energy isa Float64
    end

    @testset "Realistic Usage Patterns" begin
        # Test patterns that would be used in actual MTCR integration

        # Initialize system
        set_mtcr_lib_path!("/path/to/libmtcr.so")
        test_case_path = joinpath(@__DIR__, "test_case")
        initialize_api_wrapper(Int32(5); case_path = test_case_path)

        # Get system information
        max_species = get_max_number_of_species_wrapper()
        species_names = get_species_names_wrapper()

        @test max_species > 0
        @test length(species_names) > 0

        # Typical calculation sequence
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Calculate temperatures from state
        temps = calculate_temperatures_wrapper(rho_sp, rho_etot)
        @test temps.tt > 0

        # Calculate source terms
        sources = calculate_sources_wrapper(rho_sp, rho_etot)
        @test length(sources.drho_sp) == length(rho_sp)

        # Calculate total energy
        energy = calculate_total_energy_wrapper(temps.tt, rho_sp)
        @test energy > 0

        # Finalize system
        finalize_api_wrapper()
    end

    @testset "Error Handling Readiness" begin
        # Test that functions handle errors gracefully (important for future Fortran integration)

        # Test with library not loaded
        close_mtcr_library()
        @test_throws ErrorException get_mtcr_handle()

        # Reset library for other tests
        set_mtcr_lib_path!("/tmp/libmtcr.so")

        # Test that functions don't crash with unusual inputs
        @test_nowarn calculate_sources_wrapper([NaN], NaN)
        @test_nowarn calculate_temperatures_wrapper([Inf], Inf)
        @test_nowarn calculate_total_energy_wrapper(-1000.0, [1e-3])

        # Test with very large arrays (memory stress test)
        huge_array = zeros(Float64, 10000)
        @test_nowarn calculate_sources_wrapper(huge_array, 1e4)
    end

    @testset "Future Fortran Integration Compatibility" begin
        # Test that the interface is ready for actual Fortran library integration
        # These tests verify the structure is correct for when ccall is uncommented

        # Test that all required functions exist and have correct signatures
        @test isdefined(Main, :set_mtcr_lib_path!)
        @test isdefined(Main, :get_mtcr_handle)
        @test isdefined(Main, :get_mtcr_lib_path)
        @test isdefined(Main, :is_mtcr_loaded)
        @test isdefined(Main, :close_mtcr_library)
        @test isdefined(Main, :initialize_api_wrapper)
        @test isdefined(Main, :finalize_api_wrapper)
        @test isdefined(Main, :get_max_number_of_species_wrapper)
        @test isdefined(Main, :get_species_names_wrapper)
        @test isdefined(Main, :calculate_sources_wrapper)
        @test isdefined(Main, :calculate_temperatures_wrapper)
        @test isdefined(Main, :calculate_total_energy_wrapper)

        # Test that functions return expected types (critical for ccall compatibility)
        set_mtcr_lib_path!("/tmp/libmtcr.so")
        test_case_path = joinpath(@__DIR__, "test_case")

        # Test initialization (returns nothing now)
        @test_nowarn initialize_api_wrapper(Int32(5); case_path = test_case_path)

        max_species = get_max_number_of_species_wrapper()
        @test max_species isa Int32

        species_names = get_species_names_wrapper()
        @test species_names isa Vector{String}

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        sources = calculate_sources_wrapper(rho_sp, rho_etot)
        @test sources isa NamedTuple
        @test sources.drho_sp isa Vector{Float64}
        @test sources.drho_etot isa Float64

        temps = calculate_temperatures_wrapper(rho_sp, rho_etot)
        @test temps isa NamedTuple
        @test temps.tt isa Float64
        @test temps.tex isa Vector{Float64}
        @test temps.tvx isa Matrix{Float64}

        energy = calculate_total_energy_wrapper(1000.0, rho_sp)
        @test energy isa Float64
    end
end
