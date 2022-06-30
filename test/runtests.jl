using Test
import DifferentialProximity as DP
using StaticArrays
using LinearAlgebra
import ForwardDiff as FD
import FiniteDiff
using BenchmarkTools
using Printf
using Random: seed!
seed!(1234)

@testset "MRP" begin
    include("mrp_tests.jl")
end

@testset "qp solvers" begin
    include("qp_solver_tests.jl")
end

@testset "proximity" begin
    include("proximity_tests.jl")
end

@testset "polytopes" begin
    include("polygon_tests.jl")
end

@testset "poly deriv" begin
    include("polygon_derivs_tests.jl")
end

@testset "caps deriv" begin
    include("capsule_derivs_test.jl")
end

@testset "caps v polyg" begin
    include("capsule_v_polygon_tests.jl")
end
