using Pkg
using LinearAlgebra
using DataFrames
using Printf
using Random

println("test pf")

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 40

nbr_particles_pf = 3000
ρ = 0.0
run_sort = false

y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M=M_subjects,N=N_time,seed=100)

################################################################################
# Run pf at one # of particles
################################################################################

u_prop_old_block_1 = randn(nbr_particles_pf,N_time + 1,M_subjects)
u_resample_old_block_1 = randn(N_time,2,M_subjects)

u_prop_new = randn(nbr_particles_pf,N_time + 1,M_subjects)
u_resample_new = randn(N_time,2,M_subjects)

# correlat with old standard normal numbers
u_prop_tilde = ρ*u_prop_old_block_1 + sqrt(1-ρ^2)*u_prop_new
u_resample_tilde = ρ*u_resample_old_block_1 + sqrt(1-ρ^2)*u_resample_new

@time  cpf(y, σ_ϵ, ϕ, dt, u_prop_old_block_1, u_resample_old_block_1, nbr_particles_pf, run_sort)


function time_pf()

    for i = 1:100


        u_prop_old_block_1 = randn(nbr_particles_pf,N_time + 1,M_subjects)
        u_resample_old_block_1 = randn(N_time,2,M_subjects)

        u_prop_new = randn(nbr_particles_pf,N_time + 1,M_subjects)
        u_resample_new = randn(N_time,2,M_subjects)

        # correlat with old standard normal numbers
        u_prop_tilde = ρ*u_prop_old_block_1 + sqrt(1-ρ^2)*u_prop_new
        u_resample_tilde = ρ*u_resample_old_block_1 + sqrt(1-ρ^2)*u_resample_new

        ll =  cpf(y, σ_ϵ, ϕ, dt, u_prop_old_block_1, u_resample_old_block_1, nbr_particles_pf, run_sort)

    end

end

runtime_pf = @elapsed time_pf()

runtime_pf = runtime_pf/100

################################################################################
# Find correct nbr of particels
################################################################################


function loglik_est(nbr_particles_pf, nbr_pf_eval=100)

    ll = zeros(nbr_pf_eval)

    for i = 1:nbr_pf_eval
        u_prop_old_block_1 = randn(nbr_particles_pf,N_time + 1,M_subjects)
        u_resample_old_block_1 = randn(N_time,2,M_subjects)

        u_prop_new = randn(nbr_particles_pf,N_time + 1,M_subjects)
        u_resample_new = randn(N_time,2,M_subjects)

        # correlat with old standard normal numbers
        u_prop_tilde = ρ*u_prop_old_block_1 + sqrt(1-ρ^2)*u_prop_new
        u_resample_tilde = ρ*u_resample_old_block_1 + sqrt(1-ρ^2)*u_resample_new

        ll[i] = sum(cpf(y, σ_ϵ, ϕ, dt, u_prop_old_block_1, u_resample_old_block_1, nbr_particles_pf, run_sort))

    end

    return ll

end


nbr_pf_eval = 100
nbr_particles_pf = [100, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]

loglik_matrix = zeros(nbr_pf_eval, length(nbr_particles_pf))
run_times = zeros(length(nbr_particles_pf))

run_times[1]  = @elapsed loglik_matrix[:,1] = loglik_est(nbr_particles_pf[1], nbr_pf_eval)

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
PyPlot.savefig("figures/avg_loglik_val_vs_nbr_particels_pf.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles_pf, var(loglik_matrix, dims = 1)[:], "--*")
PyPlot.plot(nbr_particles_pf, 2*ones(length(nbr_particles_pf),1), "k")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Var of loglik est.")
PyPlot.savefig("figures/loglik_var_vs_nbr_particles_pf.png")

PyPlot.figure(figsize=(20,10))
PyPlot.plot(nbr_particles_pf, run_times/100, "--*")
PyPlot.xlabel("Nbr particles.")
PyPlot.ylabel("Run time (sec)")
PyPlot.savefig("figures/run_time_vs_nbr_particls_pf.png")
