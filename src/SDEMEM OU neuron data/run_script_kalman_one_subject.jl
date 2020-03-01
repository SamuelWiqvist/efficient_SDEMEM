# script to run inference for the OU SDEMEM model

using Pkg
using LinearAlgebra
using DataFrames

println("start run script for kalman")

# load functions
include(pwd()*"/src/SDEMEM OU neuron data/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU neuron data/mcmc.jl")

M_subjects = 1

#seed = parse(Int,ARGS[1])

seed = 100

y,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M_subjects,seed)

job = string(M_subjects)

# run MH-Gibbs
R = 100000 #500000
burn_in = 20000 #250000

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


# we should transform this parameter, to have the "walk" in a better space
startval_log_σ_ϵ = log(0.2)

# set time step
dt = diff(y[1].mSec)[1]

# estimate parameters using exact Gibbs sampling

ll = @time kalman(y, exp(startval_log_σ_ϵ), startval_ϕ, dt)
sum(ll)



Σ_i_σ_ϵ = 1^2

Σ_i_ϕ = Matrix{Float64}[]
for i in 1:M_subjects; push!(Σ_i_ϕ, [0.04 0 0;0 0.4 0;0 0 0.04]); end


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

run_time_kalman = @elapsed chain_ϕ_kalman, chain_σ_ϵ_kalman, chain_η_kalman, accept_vec_kalman, loglik_vec = gibbs_exact(R,
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
                                                                                                                         prior_parameters_σ_ϵ);


println(run_time_kalman)
println(sum(accept_vec_kalman[1,:])/(M_subjects*R))
println(sum(accept_vec_kalman[2,:])/R)
println(sum(accept_vec_kalman[3,:])/(3*R))

sim_data = zeros(4,1)

sim_data[1] = run_time_kalman
sim_data[2] = sum(accept_vec_kalman[1,:])/(M_subjects*R)
sim_data[3] = sum(accept_vec_kalman[2,:])/R
sim_data[4] = sum(accept_vec_kalman[3,:])/(3*R)


chain_ϕ_export = zeros(M_subjects*3, R)

idx = 0
for m = 1:M_subjects
    for j = 1:3
        global idx = idx + 1
        chain_ϕ_export[idx,:] = chain_ϕ_kalman[m,j,:]
    end
end


loglik_data = zeros(R,1)
loglik_data[:] = loglik_vec
DataFrame(loglik_data)

CSV.write("data/SDEMEM OU neuron data/kalman/sim_data_"*string(seed)*"_"*job*".csv", DataFrame(sim_data))
CSV.write("data/SDEMEM OU neuron data/kalman/chain_sigma_epsilon_"*string(seed)*"_"*job*".csv", DataFrame(chain_σ_ϵ_kalman'))
CSV.write("data/SDEMEM OU neuron data/kalman/chain_eta_"*string(seed)*"_"*job*".csv", DataFrame(chain_η_kalman'))
CSV.write("data/SDEMEM OU neuron data/kalman/chain_phi_"*string(seed)*"_"*job*".csv", DataFrame(chain_ϕ_export'))
CSV.write("data/SDEMEM OU neuron data/kalman/loglik_data_"*string(seed)*"_"*job*".csv", DataFrame(loglik_data))


println(Threads.nthreads())

println("end run script")
