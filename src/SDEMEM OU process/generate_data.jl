# script to run inference for the OU SDEMEM model
using Pkg
using LinearAlgebra
using DataFrames

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
#include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

#M_subjects = parse(Int,ARGS[1])
#N_time = parse(Int,ARGS[2])

N_time = 200
M_subjects = 40
seed = 100

y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M=M_subjects,N=N_time,seed=seed)


ϕ




exp.(ϕ)




using PyPlot

PyPlot.figure()
h = PyPlot.plt[:hist](ϕ[:,1],10)

PyPlot.figure()
h = PyPlot.plt[:hist](ϕ[:,2],10)

PyPlot.figure()
h = PyPlot.plt[:hist](ϕ[:,3],10)

PyPlot.figure()
h = PyPlot.plt[:hist](exp.(ϕ[:,1]),10)

PyPlot.figure()
h = PyPlot.plt[:hist](exp.(ϕ[:,2]),10)

PyPlot.figure()
h = PyPlot.plt[:hist](exp.(ϕ[:,3]),10)
h = PyPlot.plot(σ_ϵ, 0.12, "k*")


PyPlot.figure()
PyPlot.plot(Array(x'))

PyPlot.figure()
PyPlot.plot(Array(y'))

# prior predictive simulations

using Distributions

prior_dist_σ_ϵ = Gamma(prior_parameters_σ_ϵ[1],1/prior_parameters_σ_ϵ[2])

prior_dist_τ_1 = Gamma(prior_parameters_η[1,3],1/prior_parameters_η[1,4])
prior_dist_τ_2 = Gamma(prior_parameters_η[2,3],1/prior_parameters_η[2,4])
prior_dist_τ_3 = Gamma(prior_parameters_η[3,3],1/prior_parameters_η[3,4])


prior_dist_τ_2 = Gamma(2,2)
tau_2_samples = rand(prior_dist_τ_2,500)
PyPlot.figure()
h = PyPlot.plt[:hist](tau_2_samples,50)

# check loglik est val

for i in 1:10
    Random.seed!(i)

    σ_ϵ_new = rand(prior_dist_σ_ϵ)

    τ_1_new  = rand(prior_dist_τ_1)
    τ_2_new  = rand(prior_dist_τ_2)
    τ_3_new  = rand(prior_dist_τ_3)

    μ_1_new  = rand(Normal(prior_parameters_η[1,1], sqrt(1/(prior_parameters_η[1,2]*τ_1_new))))
    μ_2_new  = rand(Normal(prior_parameters_η[2,1], sqrt(1/(prior_parameters_η[2,2]*τ_2_new))))
    μ_3_new  = rand(Normal(prior_parameters_η[3,1], sqrt(1/(prior_parameters_η[3,2]*τ_3_new))))

    log_θ_1_new = μ_1_new .+ sqrt(1/τ_1_new)*randn(M_subjects)
    log_θ_2_new = μ_1_new .+ sqrt(1/τ_2_new)*randn(M_subjects)
    log_θ_3_new = μ_1_new .+ sqrt(1/τ_3_new)*randn(M_subjects)

    η_new = [μ_1_new; μ_2_new; μ_3_new; τ_1_new; τ_2_new; τ_3_new]
    ϕ_new = [log_θ_1_new log_θ_2_new log_θ_3_new]


    y,x,t_vec = generate_data(N_time, M_subjects, η_new, ϕ_new, σ_ϵ_new,dt)

    PyPlot.figure()
    PyPlot.plot(Array(y'))

end
