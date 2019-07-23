using Pkg
using LinearAlgebra
using DataFrames

println("test corr-pf")

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 10
nbr_particles = 500

ρ = 0.9

θ_1 = 0.5
σ_ϵ = 0.2
θ_3 = 0.1

y,x,t_vec,dt,η,κ,σ_ϵ,ϕ,prior_η,prior_κ,prior_σ_ϵ = set_up(σ_ϵ = σ_ϵ, θ_3 = θ_3,θ_1=θ_1, M=M_subjects,N=N_time)

κ = [θ_1,θ_3]

function time_pf()

    for i = 1:1000
        ll = @time sum(pf(y, κ, σ_ϵ, ϕ, dt, nbr_particles))
    end

end

function time_corr_pf()

    for i = 1:1000
        ll = sum(cpf(y, κ, σ_ϵ, ϕ, dt, randn(nbr_particles,N_time + 1,M_subjects), randn(N_time,2,M_subjects), nbr_particles))
    end

end


runtime_pf = @elapsed time_pf()

runtime_pf = runtime_pf/1000

runtime_cpf = @elapsed time_corr_pf()

runtime_cpf = runtime_cpf/1000

runtime_pf/runtime_cpf
