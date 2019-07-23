# script to run inference for the OU SDEMEM model
using Pkg
using PyPlot
using LinearAlgebra
using KernelDensity
using DataFrames

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 40

σ_ϵ = 0.2

y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ = set_up(σ_ϵ = σ_ϵ,M=M_subjects,N=N_time)


PyPlot.figure()
PyPlot.plot(t_vec, Array(x'))

PyPlot.figure()
PyPlot.plot(t_vec, Array(y'))

#=
yt = Matrix(CSV.read("data/SDEMEM OU/data_ctms_small_t3.csv"; allowmissing=:auto))
y = Array(yt[:,2]')
PyPlot.figure()
PyPlot.plot(Array(y[1,:]))
=#

# test kalman filter

# TODO: check kalman filter the path tracking is quite bad...
a = rand(2)


ϕ

σ_ϵ


@time loglik, y_hat, x_hat = kalman(y, σ_ϵ, ϕ, dt, true)

data = zeros(N_time,2)
data[:,1] = t_vec
data[:,2] = y

CSV.write("data/SDEMEM OU/data_julia_large_theta_3.csv", DataFrame(data))


sum(loglik)


sum(loglik)

PyPlot.figure()
PyPlot.plot(Array(x_hat'))

PyPlot.figure()
PyPlot.plot(Array(y_hat'))

PyPlot.figure()
PyPlot.plot(Array(x'))

PyPlot.figure()
PyPlot.plot(t_vec, Array(y[10,:]), "b")
PyPlot.plot(t_vec, Array(y_hat[10,:]), "r--")
