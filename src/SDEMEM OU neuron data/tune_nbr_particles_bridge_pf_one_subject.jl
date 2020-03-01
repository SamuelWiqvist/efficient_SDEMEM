using Pkg
using LinearAlgebra
using DataFrames
import Statistics.mean
import Statistics.std
using Printf
using CSV

# load functions
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


nbr_particles = 100*ones(Int64, M_subjects)
ρ = 0.0

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

    run_time = @elapsed ll1 =  bridge_pf(y, exp(log_σ_ϵ), ϕ, dt, u_prop_old, u_resample_old, nbr_particles, false)

    return ll1, run_time

end

################################################################################
# run at one specific # particles
################################################################################

#=

nbr_particels = 1000

ll1, run_time = @time run_bridge_pf(nbr_particels)

@printf "Runtime Bridge filter: %.4f\n" run_time

nbr_pf_eval = 50
ll_vec = zeros(nbr_pf_eval)
run_time_vec = zeros(nbr_pf_eval)
ll_all = zeros(M_subjects, nbr_pf_eval)

for i = 1:nbr_pf_eval

    println(i)

    ll1, run_time = run_bridge_pf(nbr_particels)
    ll_all[:,i] = ll1
    ll_vec[i] = sum(ll1)
    run_time_vec[i] = run_time

end

@printf "----------------\n"
@printf "Bridge filter\n"
@printf "Nbr particles: %.2f\n" nbr_particels
@printf "Runtime: %.4f\n" mean(run_time_vec)
@printf "Loglik (mean): %.4f\n" mean(ll_vec)
@printf "Loglik (var): %.4f\n" var(ll_vec)
@printf "Loglik (std): %.4f\n" std(ll_vec)
@printf "Loglik (var) ind ll est:\n"
println(var(ll_all, dims = 2))

=#

################################################################################
# run over multiple number of # particles
################################################################################

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
#PyPlot.savefig("figures/neuronal_avg_loglik_val_vs_nbr_particels_pf.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles_pf, var(loglik_matrix, dims = 1)[:], "--*")
PyPlot.plot(nbr_particles_pf, 2*ones(length(nbr_particles_pf),1), "k")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Var of loglik est.")
#PyPlot.savefig("figures/neuronal_loglik_var_vs_nbr_particles_pf.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles_pf, run_times/100, "--*")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Run time (sec)")
#PyPlot.savefig("figures/neuronal_run_time_vs_nbr_particls_pf.png")
