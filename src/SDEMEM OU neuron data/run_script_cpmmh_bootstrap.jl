# script to run inference for the OU SDEMEM model
using Pkg
using LinearAlgebra
using DataFrames

println("start run script for cpmmh")

# load functions
include(pwd()*"/src/SDEMEM OU neuron data/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU neuron data/mcmc.jl")

M_subjects = 100

ρ = parse(Float64,ARGS[1])
nbr_particles = parse(Int,ARGS[2])*ones(Int64, M_subjects)
seed = parse(Int,ARGS[3])


y, prior_parameters_η, prior_parameters_σ_ϵ = set_up(M_subjects, seed)


# set number of obs for each data set
T_vec = zeros(Int64, M_subjects)
for i = 1:M_subjects; T_vec[i] = length(y[i].mV); end

# run MH-Gibbs

R = 35000 #15000 #10000
burn_in = 2000

job = string(M_subjects)*"_"*string(ρ)*"_"*string(nbr_particles[1])*"_"*"bootstrap"

startval_ϕ = ones(M_subjects,3)

for j = 1:3

    μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

    startval_ϕ[:,j,1] = μ_0_j*ones(M_subjects)

end


startval_η = zeros(6)

for j = 1:3

    # get paramaters
    μ_0_j, M_0_j, α_j, β_j = prior_parameters_η[j,:]

    startval_η[j] = μ_0_j # set start value to mean
    startval_η[j+3] = (α_j - 1)/β_j

end


startval_log_σ_ϵ = log(0.2)

# set time step
dt = diff(y[1].mSec)[1]

u_prop_old = [randn(nbr_particles[1], T_vec[1]+1)]
u_resample_old = [randn(T_vec[1], 2)]


for i in 2:M_subjects;
    append!(u_prop_old, [randn(nbr_particles[i],T_vec[i]+1)])
    append!(u_resample_old, [randn(T_vec[i],2)])
end

ll = @time cpf(y, exp(startval_log_σ_ϵ), startval_ϕ, dt, u_prop_old, u_resample_old, nbr_particles, true)
sum(ll)


Σ_i_σ_ϵ = 1^2

Σ_i_ϕ = Matrix{Float64}[]
for i in 1:M_subjects; push!(Σ_i_ϕ, [0.004 0 0;0 0.04 0;0 0 0.004]); end


γ_ϕ_0 = 1.
γ_σ_ϵ_0 = 1.
μ_i_ϕ = startval_ϕ
μ_i_σ_ϵ = startval_log_σ_ϵ
α_star_ϕ = 0.25
α_star_σ_ϵ = 0.25
log_λ_i_ϕ = log.(2.4/sqrt(3)*ones(M_subjects))
log_λ_i_σ_ϵ = log(2.4)
update_interval = 1
α_power = 0.7
start_update = 100


# set random numbers
Random.seed!(seed)

run_time_cpmmh = @elapsed chain_ϕ_cpmmh, chain_log_σ_ϵ_cpmmh, chain_η_cpmmh, accept_vec_cpmmh, loglik_vec = gibbs_cpmmh(R,
                                                                                                                        y,
                                                                                                                        dt,
                                                                                                                        Σ_i_σ_ϵ,
                                                                                                                        Σ_i_ϕ,
                                                                                                                        γ_ϕ_0,
                                                                                                                        γ_σ_ϵ_0,
                                                                                                                        μ_i_ϕ,
                                                                                                                        μ_i_σ_ϵ,
                                                                                                                        α_star_ϕ,
                                                                                                                        α_star_σ_ϵ,
                                                                                                                        log_λ_i_ϕ,
                                                                                                                        log_λ_i_σ_ϵ,
                                                                                                                        update_interval,
                                                                                                                        start_update,
                                                                                                                        α_power,
                                                                                                                        startval_ϕ,
                                                                                                                        startval_log_σ_ϵ,
                                                                                                                        startval_η,
                                                                                                                        prior_parameters_η,
                                                                                                                        prior_parameters_σ_ϵ,
                                                                                                                        nbr_particles,
                                                                                                                        T_vec,
                                                                                                                        ρ,
                                                                                                                        cpf);

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

loglik_data = zeros(R,1)
loglik_data[:] = loglik_vec
DataFrame(loglik_data)

CSV.write("data/SDEMEM OU neuron data/cpmmh/sim_data_"*string(seed)*"_"*job*".csv", DataFrame(sim_data))
CSV.write("data/SDEMEM OU neuron data/cpmmh/chain_sigma_epsilon_"*string(seed)*"_"*job*".csv", DataFrame(chain_log_σ_ϵ_cpmmh'))
CSV.write("data/SDEMEM OU neuron data/cpmmh/chain_eta_"*string(seed)*"_"*job*".csv", DataFrame(chain_η_cpmmh'))
CSV.write("data/SDEMEM OU neuron data/cpmmh/chain_phi_"*string(seed)*"_"*job*".csv", DataFrame(chain_ϕ_export'))
CSV.write("data/SDEMEM OU neuron data/cpmmh/loglik_data_"*string(seed)*"_"*job*".csv", DataFrame(loglik_data))
