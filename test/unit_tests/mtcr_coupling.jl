using Test

@testset verbose=true "MTCR Module Tests" begin
    @testset "Data Conversion" begin
        include("mtcr/data_conversion.jl")
    end

    @testset "Fortran Wrapper" begin
        include("mtcr/fortran_wrapper.jl")
    end

    @testset "MTCR Configuration" begin
        include("mtcr/mtcr_config.jl")
    end

    # @testset "MTCR Solver" begin
    #     include("mtcr/mtcr_solver.jl")
    # end
end
