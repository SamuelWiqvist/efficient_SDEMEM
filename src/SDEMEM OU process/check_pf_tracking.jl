# script to run inference for the OU SDEMEM model
using Pkg
using PyPlot
using LinearAlgebra
using KernelDensity
using DataFrames

# TODO: Which parameters can I log-transform?!?

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 40

θ_1 = 0.5
σ_ϵ = 0.2
θ_3 = 0.1

y,x,t_vec,dt,η,κ,σ_ϵ,ϕ,prior_η,prior_κ,prior_σ_ϵ = set_up(σ_ϵ = σ_ϵ, θ_3 = θ_3,θ_1=θ_1, M=M_subjects,N=N_time)

PyPlot.figure()
PyPlot.plot(Array(x'))

PyPlot.figure()
PyPlot.plot(Array(y'))

yt = Matrix(CSV.read("data/SDEMEM OU/data_ctms_small_t3.csv"; allowmissing=:auto))
y = Array(yt[:,2]')
PyPlot.figure()
PyPlot.plot(Array(y[1,:]))

@time loglik, w,x = pf(y, κ, σ_ϵ, ϕ, dt, 1000, true)
@profiler loglik, w,x = pf(y, κ, σ_ϵ, ϕ, dt, 1000, true)

sum(loglik)

PyPlot.figure()
PyPlot.plot(Array(x[:,end,:]'))

PyPlot.figure()
PyPlot.plot(Array(y[1,:]), "b")
PyPlot.plot(Array(x[1,end,:]), "r--")


loglik_m = zeros(M_subjects,500)

@time begin
for i in 1:500
    loglik_m[:,i] = pf(y, κ, σ_ϵ, ϕ, dt, 10)
end
end

std(loglik_m,dims=2)
std(sum(loglik_m, dims = 1))
