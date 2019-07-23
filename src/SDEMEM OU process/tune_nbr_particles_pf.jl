using Pkg
using LinearAlgebra
using DataFrames

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


nbr_pf_eval = 100

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

mean(ll)
var(ll)
std(ll)
