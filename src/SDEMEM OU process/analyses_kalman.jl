# load packages
using Pkg
using PyPlot
using LinearAlgebra
using KernelDensity
using CSV
using Distributions

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")

N_time = 200
M_mixtures = 100
θ_3 = 0.1

y,x,t_vec,dt,η,κ,σ_ϵ,ϕ,prior_η,prior_κ,prior_σ_ϵ = set_up(θ_3 = θ_3, M=M_mixtures,N=N_time)

job = string(M_mixtures)*"_"*string(N_time)

chain_1 = Matrix(CSV.read("data/SDEMEM OU/kalman/chain_1_"*job*".csv"; allowmissing=:auto))'
chain_2 = Matrix(CSV.read("data/SDEMEM OU/kalman/chain_2_"*job*".csv"; allowmissing=:auto))'
chain_3 = Matrix(CSV.read("data/SDEMEM OU/kalman/chain_3_"*job*".csv"; allowmissing=:auto))'

loglik2_loglik_y = Matrix(CSV.read("data/SDEMEM OU/kalman/loglik2_loglik_y_"*job*".csv"; allowmissing=:auto))
loglik2_loglik_ϕ = Matrix(CSV.read("data/SDEMEM OU/kalman/loglik2_loglik_ϕ_"*job*".csv"; allowmissing=:auto))

#loglik1 = Matrix(CSV.read("data/SDEMEM OU/kalman/loglik1_"*job*".csv"; allowmissing=:auto))
#loglik3 = Matrix(CSV.read("data/SDEMEM OU/kalman/loglik3_"*job*".csv"; allowmissing=:auto))

accept_vec = Matrix(CSV.read("data/SDEMEM OU/kalman/accept_vec_"*job*".csv"; allowmissing=:auto))'

burn_in = 10000
R = size(chain_1,2)

# plot results
PyPlot.figure()
PyPlot.plot(chain_1[1,:])
PyPlot.plot(σ_ϵ*ones(R), "k--")

PyPlot.figure()
PyPlot.plot(chain_1[2,:])
PyPlot.plot(κ[1]*ones(R), "k--")

PyPlot.figure()
PyPlot.plot(chain_1[3,:])
PyPlot.plot(κ[2]*ones(R), "k--")

for i in 1:M_mixtures
    PyPlot.figure()
    PyPlot.plot(chain_2[i,:])
    PyPlot.plot(ϕ[i]*ones(R), "k--")
end

PyPlot.figure()
PyPlot.plot(chain_3[1,:])
PyPlot.plot(η[1]*ones(R), "k--")

PyPlot.figure()
PyPlot.plot(chain_3[2,:])
PyPlot.plot(η[2]*ones(R), "k--")

# plot results after burn-in
PyPlot.figure()
PyPlot.plot(chain_1[1,burn_in:end])
PyPlot.plot(σ_ϵ*ones(R), "k--")

PyPlot.figure()
PyPlot.plot(chain_1[2,burn_in:end])
PyPlot.plot(κ[1]*ones(R), "k--")

PyPlot.figure()
PyPlot.plot(chain_1[3,burn_in:end])
PyPlot.plot(κ[2]*ones(R), "k--")

for i in 1:M_mixtures
    PyPlot.figure()
    PyPlot.plot(chain_2[i,burn_in:end])
    PyPlot.plot(ϕ[i]*ones(R), "k--")
end

PyPlot.figure()
PyPlot.plot(chain_3[1,burn_in:end])
PyPlot.plot(η[1]*ones(R), "k--")

PyPlot.figure()
PyPlot.plot(chain_3[2,burn_in:end])
PyPlot.plot(η[2]*ones(R), "k--")

#plot loglik
PyPlot.figure()
PyPlot.plot(loglik1)

PyPlot.figure()
PyPlot.plot(Array(loglik2_loglik_y))

PyPlot.figure()
PyPlot.plot(Array(loglik2_loglik_ϕ))

PyPlot.figure()
PyPlot.plot(loglik3)

# acceptance rate
sum(accept_vec[1,:]/R*100)

sum(accept_vec[2,:]/(M_mixtures*R)*100)

sum(accept_vec[3,:]/R*100)

# posterior means
mean(chain_1[1,burn_in:end])
mean(chain_1[2,burn_in:end])
mean(chain_2[1,burn_in:end])
mean(chain_1[3,burn_in:end])


# plot posterior

# posterior σ_ϵ
prior_grid = LinRange(0,1,100)
prior = pdf.(Gamma(prior_σ_ϵ[1].parameters[1],prior_σ_ϵ[1].parameters[2]), prior_grid)
posterior = kde(chain_1[1,burn_in:end])

PyPlot.figure()
PyPlot.plot(prior_grid,prior, "g")
PyPlot.plot(posterior.x,posterior.density, "b")
PyPlot.plot((σ_ϵ, σ_ϵ), (0, maximum(posterior.density)), "k")

# posterior θ_1
prior_grid = LinRange(0,2,100)
prior = pdf.(Gamma(prior_κ[1].parameters[1],prior_κ[1].parameters[2]), prior_grid)
posterior = kde(chain_1[2,burn_in:end])

PyPlot.figure()
PyPlot.plot(prior_grid,prior, "g")
PyPlot.plot(posterior.x,posterior.density, "b")
PyPlot.plot((κ[1], κ[1]), (0, maximum(posterior.density)), "k")

# posterior θ_3
prior_grid = LinRange(0,2,100)
prior = pdf.(Gamma(prior_κ[2].parameters[1],prior_κ[2].parameters[2]), prior_grid)
posterior = kde(chain_1[3,burn_in:end])

PyPlot.figure()
PyPlot.plot(prior_grid,prior, "g")
PyPlot.plot(posterior.x,posterior.density, "b")
PyPlot.plot((κ[2], κ[2]), (0, maximum(posterior.density)), "k")

# posterior ϕ

for i = 1:M_mixtures

    posterior = kde(chain_2[i,burn_in:end])

    PyPlot.figure()
    PyPlot.plot(posterior.x,posterior.density, "b")
    PyPlot.plot((ϕ[i], ϕ[i]), (0, maximum(posterior.density)), "k")
end

# posterior μ_θ_2
prior_grid = LinRange(-10,10,100)
prior = pdf.(Normal(prior_η[1].parameters[1],prior_η[1].parameters[2]), prior_grid)
posterior = kde(chain_3[1,burn_in:end])

PyPlot.figure()
PyPlot.plot(prior_grid,prior, "g")
PyPlot.plot(posterior.x,posterior.density, "b")
PyPlot.plot((η[1], η[1]), (0, maximum(posterior.density)), "k")

# posterior σ_θ_2
prior_grid = LinRange(0,5,100)
prior = pdf.(Gamma(prior_η[2].parameters[1],prior_η[2].parameters[2]), prior_grid)
posterior = kde(chain_3[2,burn_in:end])

PyPlot.figure()
PyPlot.plot(prior_grid,prior, "g")
PyPlot.plot(posterior.x,posterior.density, "b")
PyPlot.plot((η[2], η[2]), (0, maximum(posterior.density)), "k")
