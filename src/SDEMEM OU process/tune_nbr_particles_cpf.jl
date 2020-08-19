using Pkg
using LinearAlgebra
using DataFrames
using Printf
using Random

println("test corr-pf")

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 40
nbr_particles_cor = 50

y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M=M_subjects,N=N_time,seed=100)


# estimate correlation
nbr_pf_eval_corr = 100

ll1 = zeros(nbr_pf_eval_corr)
ll2 = zeros(nbr_pf_eval_corr)

ρ = 0.999
corr_ll_vec = zeros(25)

Random.seed!(1)

for j in 1:25

    for i = 1:nbr_pf_eval_corr

        u_prop_old_block_1 = randn(nbr_particles_cor,N_time + 1,M_subjects)
        u_resample_old_block_1 = randn(N_time,2,M_subjects)

        u_prop_new = randn(nbr_particles_cor,N_time + 1,M_subjects)
        u_resample_new = randn(N_time,2,M_subjects)

        # correlat with old standard normal numbers
        u_prop_tilde = ρ*u_prop_old_block_1 + sqrt(1-ρ^2)*u_prop_new
        u_resample_tilde = ρ*u_resample_old_block_1 + sqrt(1-ρ^2)*u_resample_new

        ll1[i] = sum(cpf(y, σ_ϵ, ϕ, dt, u_prop_old_block_1, u_resample_old_block_1, nbr_particles_cor,true))
        ll2[i] = sum(cpf(y, σ_ϵ, ϕ, dt, u_prop_tilde, u_resample_tilde, nbr_particles_cor,true))

    end

    corr_ll_vec[j] = cor(ll1,ll2)

end

corr_ll = mean(corr_ll_vec)

var_loglik_target = 2.16^2/(1-corr_ll^2)

#nbr_pf_eval = 100
#nbr_particles = 50


function loglik_est(nbr_particles, nbr_pf_eval=100)

    ll = zeros(nbr_pf_eval)

    for i = 1:nbr_pf_eval
        ll[i] = sum(cpf(y, σ_ϵ, ϕ, dt, randn(nbr_particles,N_time + 1,M_subjects), randn(N_time,2,M_subjects), nbr_particles, true))
    end

    return ll

end


nbr_pf_eval = 100
nbr_particles = [20,50,100,150,200,250,500]

loglik_matrix = zeros(nbr_pf_eval, length(nbr_particles))
run_times = zeros(length(nbr_particles))

run_times[1]  = @elapsed loglik_matrix[:,1] = loglik_est(nbr_particles[1], nbr_pf_eval)

for i in 1:length(nbr_particles)
    @printf "Current nbr of particels: %.2f\n" nbr_particles[i]
    run_times[i]  = @elapsed loglik_matrix[:,i]  = loglik_est(nbr_particles[i], nbr_pf_eval)
end

using PyPlot

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles, mean(loglik_matrix, dims = 1)[:], "--*")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Avg. loglik est. val.")
PyPlot.savefig("figures/avg_loglik_val_vs_nbr_particels_cpf_099.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles, var(loglik_matrix, dims = 1)[:], "--*")
PyPlot.plot(nbr_particles, var_loglik_target*ones(length(nbr_particles),1), "k")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Var of loglik est.")
PyPlot.savefig("figures/loglik_var_vs_nbr_particles_cpf_099.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles, run_times/100, "--*")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Run time (sec)")
PyPlot.savefig("figures/run_time_vs_nbr_particls_cpf_099.png")
