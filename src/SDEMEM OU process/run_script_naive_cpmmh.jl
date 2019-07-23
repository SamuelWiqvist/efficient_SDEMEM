# script to run inference for the OU SDEMEM model
using Pkg
using LinearAlgebra
using DataFrames

println("start run script")

# load functions
include(pwd()*"/src/SDEMEM OU process/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU process/mcmc.jl")

N_time = 200
M_subjects = 40

nbr_particles = parse(Int,ARGS[1])
ρ = parse(Float64,ARGS[2])
seed = parse(Int,ARGS[3])


y,x,t_vec,dt,η,σ_ϵ,ϕ,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M=M_subjects,N=N_time,seed=seed)

# run MH-Gibbs

R = 15000
burn_in = 500

job = string(M_subjects)*"_"*string(N_time)*"_"*string(nbr_particles)*"_"*string(ρ)

cov_ϕ = [0.04;0.08;0.04]
cov_σ_ϵ = 0.02


# hard coded start values
#startval_ϕ = ϕ
#startval_σ_ϵ = σ_ϵ

startval_ϕ = ones(M_subjects,3)

for j = 1:3

    μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

    startval_ϕ[:,j,1] = μ_0_j*ones(M_subjects)

end

startval_σ_ϵ = 0.5


# estimate parameters using exact Gibbs sampling

# set random numbers
Random.seed!(seed)

run_time_cpmmh = @elapsed chain_ϕ_cpmmh, chain_σ_ϵ_cpmmh, chain_η_cpmmh, accept_vec_cpmmh = naive_gibbs_cpmmh(R,
                                                                                                              y,
                                                                                                              dt,
                                                                                                              cov_ϕ,
                                                                                                              cov_σ_ϵ,
                                                                                                              startval_ϕ,
                                                                                                              startval_σ_ϵ,
                                                                                                              prior_parameters_η,
                                                                                                              prior_parameters_σ_ϵ,
                                                                                                              nbr_particles,
                                                                                                              ρ);

println(run_time_cpmmh)
println(sum(accept_vec_cpmmh[1,:])/(M_subjects*R))
println(sum(accept_vec_cpmmh[2,:])/R)
println(sum(accept_vec_cpmmh[3,:])/(3*R))

sim_data = zeros(4,1)
sim_data[1] = run_time_cpmmh
sim_data[2] = sum(accept_vec_cpmmh[1,:])/(M_subjects*R)
sim_data[3] = sum(accept_vec_cpmmh[2,:])/R
sim_data[4] = sum(accept_vec_cpmmh[3,:])/(3*R)


# Save results

chain_ϕ_export = zeros(M_subjects*3, R)

idx = 0
for m = 1:M_subjects
    for j = 1:3
        global idx = idx + 1
        chain_ϕ_export[idx,:] = chain_ϕ_cpmmh[m,j,:]
    end
end

CSV.write("data/SDEMEM OU/naive_cpmmh/sim_data_"*string(seed)*"_"*job*".csv", DataFrame(sim_data))
CSV.write("data/SDEMEM OU/naive_cpmmh/chain_sigma_epsilon_"*string(seed)*"_"*job*".csv", DataFrame(chain_σ_ϵ_cpmmh'))
CSV.write("data/SDEMEM OU/naive_cpmmh/chain_eta_"*string(seed)*"_"*job*".csv", DataFrame(chain_η_cpmmh'))
CSV.write("data/SDEMEM OU/naive_cpmmh/chain_phi_"*string(seed)*"_"*job*".csv", DataFrame(chain_ϕ_export'))
