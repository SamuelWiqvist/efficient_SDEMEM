using Pkg
using LinearAlgebra
using DataFrames
import Statistics.mean
import Statistics.std
import Statistics.cor
using Printf
using CSV


include(pwd()*"/src/SDEMEM OU neuron data/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU neuron data/mcmc.jl")

seed = 100

M_subjects = 1

y, prior_parameters_η, prior_parameters_σ_ϵ = set_up(M_subjects, seed)

startval_ϕ = ones(M_subjects,3)

for j = 1:3

    μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

    startval_ϕ[:,j,1] = μ_0_j*ones(M_subjects)

end

startval_log_σ_ϵ = log(0.5)

# set time step
dt = diff(y[1].mSec)[1]

# posterior mean ests from kalman filter
log_σ_ϵ = -3.60203
ϕ  = [-3.26986 -0.888105 -0.780767]
nbr_particles_cor = 10

# estimate correlation
nbr_pf_eval_corr = 100

ll1 = zeros(nbr_pf_eval_corr)
ll2 = zeros(nbr_pf_eval_corr)

ρ = 0.99

for i = 1:nbr_pf_eval_corr

    nbr_particles = nbr_particles_cor*ones(Int64, M_subjects)

    # set number of obs for each data set
    T_vec = zeros(Int64, M_subjects)
    for i = 1:M_subjects; T_vec[i] = length(y[i].mV); end


    u_prop_old = [randn(nbr_particles[1], T_vec[1]+1)]
    u_resample_old = [randn(T_vec[1], 2)]

    for i in 2:M_subjects;
        append!(u_prop_old, [randn(nbr_particles[i],T_vec[i]+1)])
        append!(u_resample_old, [randn(T_vec[i],2)])
    end

    u_prop_new = [randn(nbr_particles[1], T_vec[1]+1)]
    u_resample_new = [randn(T_vec[1], 2)]

    for i in 2:M_subjects;
        append!(u_prop_new, [randn(nbr_particles[i],T_vec[i]+1)])
        append!(u_resample_new, [randn(T_vec[i],2)])
    end

    # correlat with old standard normal numbers
    u_prop_tilde = ρ*u_prop_old + sqrt(1-ρ^2)*u_prop_new
    u_resample_tilde = ρ*u_resample_old + sqrt(1-ρ^2)*u_resample_new

    ll1[i] = sum(bridge_pf(y, exp(log_σ_ϵ), ϕ, dt, u_prop_old, u_resample_old, nbr_particles,true))
    ll2[i] = sum(bridge_pf(y, exp(log_σ_ϵ), ϕ, dt, u_prop_tilde, u_resample_tilde, nbr_particles,true))

end

corr_ll = cor(ll1,ll2)

var_loglik_target = 2.16^2/(1-corr_ll^2)

#nbr_pf_eval = 100
#nbr_particles = 50

function run_bridge_pf(particels)

    nbr_particles = particels*ones(Int64, M_subjects)
    ρ = 0.0

    # set number of obs for each data set
    T_vec = zeros(Int64, M_subjects)
    for i = 1:M_subjects; T_vec[i] = length(y[i].mV); end

    u_prop_old = [randn(nbr_particles[1], T_vec[1]+1)]
    u_resample_old = [randn(T_vec[1], 2)]

    for i in 2:M_subjects;
        append!(u_prop_old, [randn(nbr_particles[i],T_vec[i]+1)])
        append!(u_resample_old, [randn(T_vec[i],2)])
    end

    run_time = @elapsed ll1 =  bridge_pf(y, exp(log_σ_ϵ), ϕ, dt, u_prop_old, u_resample_old, nbr_particles, true)

    return ll1, run_time

end

function loglik_est(nbr_particles_pf, nbr_pf_eval)

    ll = zeros(nbr_pf_eval)

    for i in 1:nbr_pf_eval
        ll_part, ~ = run_bridge_pf(nbr_particles_pf)
        ll[i] = sum(ll_part)
    end

    return ll

end


nbr_pf_eval = 100
nbr_particles_pf = [1, 2, 5, 10, 20, 50, 100]


loglik_matrix = zeros(nbr_pf_eval, length(nbr_particles_pf))
run_times = zeros(length(nbr_particles_pf))

i = 1
run_times[i]  = @elapsed loglik_matrix[:,i] = loglik_est(nbr_particles_pf[i], nbr_pf_eval)

Random.seed!(1)

for i in 1:length(nbr_particles_pf)
    @printf "Current nbr of particels: %.2f\n" nbr_particles_pf[i]
    run_times[i]  = @elapsed loglik_matrix[:,i]  = loglik_est(nbr_particles_pf[i], nbr_pf_eval)
end

using PyPlot


PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles_pf, mean(loglik_matrix, dims = 1)[:], "--*")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Avg. loglik est. val.")
#PyPlot.savefig("figures/avg_loglik_val_vs_nbr_particels_cpf_0999.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles_pf, var(loglik_matrix, dims = 1)[:], "--*")
PyPlot.plot(nbr_particles_pf, var_loglik_target*ones(length(nbr_particles_pf),1), "k")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Var of loglik est.")
#PyPlot.savefig("figures/loglik_var_vs_nbr_particles_cpf_0999.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles_pf, run_times/100, "--*")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Run time (sec)")
#PyPlot.savefig("figures/run_time_vs_nbr_particls_cpf_0999.png")
