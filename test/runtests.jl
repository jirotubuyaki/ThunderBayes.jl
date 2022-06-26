using ThunderBayes
using Test
using Distributions, Random
using Statistics, StatsBase
using LinearAlgebra
using DataFrames
using Plots,ColorSchemes

@testset "ThunderBayes.jl" begin
mu_1 = [5, 5]
mu_2 = [2, 2]
mu_3 = [6, -6]
sigma_1 = Matrix{Float64}(1 * I, 2, 2)
sigma_2 = Matrix{Float64}(2 * I, 2, 2)
sigma_3 = Matrix{Float64}(1 * I, 2, 2)
p_1 =  MvNormal(mu_1, sigma_1)
p_2 = MvNormal(mu_2, sigma_2)
p_3 = MvNormal(mu_3, sigma_3)
data = Matrix{BigFloat}(undef, 600, 2)
for i = 1 : 200
    data[i, :] = rand(p_1, 1)
end
for i = 201 : 400
    data[i, :] = rand(p_2, 1)
end
for i = 401 : 600
    data[i, :] = rand(p_3, 1)
end
ok = data_check(data)
@test ok == true
crp_visualize([1:2],1)
end
