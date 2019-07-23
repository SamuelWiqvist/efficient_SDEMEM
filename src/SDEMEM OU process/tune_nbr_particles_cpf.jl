using Pkg
using LinearAlgebra
using DataFrames

println("test corr-pf")

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 40
nbr_particles_cor = 3000

y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M=M_subjects,N=N_time,seed=100)


# estimate correlation
nbr_pf_eval_corr = 50

ll1 = zeros(nbr_pf_eval_corr)
ll2 = zeros(nbr_pf_eval_corr)

ρ = 0.999

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


corr_ll = cor(ll1,ll2)

var_loglik = 2.16^2/(1-corr_ll^2)

nbr_pf_eval = 100

nbr_particles = 50

ll = zeros(nbr_pf_eval)

for i = 1:nbr_pf_eval
    ll[i] = sum(cpf(y, σ_ϵ, ϕ, dt, randn(nbr_particles,N_time + 1,M_subjects), randn(N_time,2,M_subjects), nbr_particles, true))
end

mean(ll)
var(ll)
