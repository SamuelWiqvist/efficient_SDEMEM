using Pkg
using LinearAlgebra
using DataFrames
import Statistics.mean
import Statistics.std
using Printf
using CSV

println("start run script for cpmmh")

# load functions
include(pwd()*"/src/SDEMEM OU neuron data/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU neuron data/mcmc.jl")

seed = 100

M_subjects = 100

y, prior_parameters_η, prior_parameters_σ_ϵ = set_up(M_subjects, seed)

startval_ϕ = ones(M_subjects,3)

for j = 1:3

    μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

    startval_ϕ[:,j,1] = μ_0_j*ones(M_subjects)

end

startval_log_σ_ϵ = log(0.5)

# set time step
dt = diff(y[1].mSec)[1]

posterior_mean_kalman = CSV.read("data/SDEMEM OU neuron data/posterior_mean_sigma_epsilon_and_phi.csv")

log_σ_ϵ  = posterior_mean_kalman[1,1]
ϕ  = Matrix(posterior_mean_kalman[2:end,:])

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


function loglik_est(nbr_particles_pf, nbr_pf_eval)

    ll = zeros(nbr_pf_eval)

    for i in 1:nbr_pf_eval
        ll_part, ~ = run_bridge_pf(nbr_particles_pf)
        ll[i] = sum(ll_part)
    end

    return ll

end


nbr_pf_eval = 50
nbr_particles_pf = [20, 50, 100, 200, 250, 500, 1000, 2000, 3000, 4000]

loglik_matrix = zeros(nbr_pf_eval, length(nbr_particles_pf)+1)
run_times = zeros(length(nbr_particles_pf))

#i = 1
#run_times[i]  = @elapsed loglik_matrix[:,i] = loglik_est(nbr_particles_pf[i], nbr_pf_eval)

Random.seed!(1)

for i in 1:length(nbr_particles_pf)
    @printf "Current nbr of particels: %.2f\n" nbr_particles_pf[i]
    run_times[i]  = @elapsed loglik_matrix[:,i+1]  = loglik_est(nbr_particles_pf[i], nbr_pf_eval)
end


loglik_matrix[1:length(nbr_particles_pf),1] = run_times


CSV.write("data/SDEMEM OU neuron data/tune_nbr_particles_bridge_pf.csv", DataFrame(loglik_matrix))
