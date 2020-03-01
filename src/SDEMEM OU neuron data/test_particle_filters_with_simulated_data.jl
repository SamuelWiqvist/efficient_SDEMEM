# parameter settings: \nu, \lambda, \sigma, \sigma_e)=(0.001, 0.02, 0.1, 0.001,
# parameter settings: sigma_epsilon = 0.001
# start value: x_0 = 5


using Pkg
using DataFrames
import Statistics.mean
import Statistics.var
using LinearAlgebra
using CSV
using PyPlot


# load functions
include(pwd()*"/src/SDEMEM OU neuron data/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU neuron data/mcmc.jl")

data = CSV.read("data/SDEMEM OU neuron data/OUsim.dat"; header=0)

y = data[2:end,1]



dt = 1

time_grid = 1:dt:length(y)

PyPlot.figure()
PyPlot.plot(time_grid, y)


M_subjects = 1

ϕ = zeros(M_subjects,3) # lambda nu sigma

ϕ[1,:] = log.([0.001, 0.02, 0.1])

log_σ_ϵ = log(0.001)


nbr_particles = 1*ones(Int64, M_subjects)
ρ = 0.0

# set number of obs for each data set
T_vec = zeros(Int64, M_subjects)
for i = 1:M_subjects; T_vec[i] = length(y); end

u_prop_old = [randn(nbr_particles[1], T_vec[1]+1)]
u_resample_old = [randn(T_vec[1], 2)]

for i in 2:M_subjects;
    append!(u_prop_old, [randn(nbr_particles[i],T_vec[i]+1)])
    append!(u_resample_old, [randn(T_vec[i],2)])
end

m = 1


dt
# Plot paths

loglik_est, x_hat_kalman = @time kalman_filter_test(y, exp(log_σ_ϵ),  ϕ[m,:], dt, true)

PyPlot.figure()
PyPlot.plot(time_grid, y, "g")
PyPlot.plot(time_grid, x_hat_kalman, "r--")


nbr_particles_bootstrap = 10000*ones(Int64, M_subjects)

# set number of obs for each data set
T_vec = zeros(Int64, M_subjects)
for i = 1:M_subjects; T_vec[i] = length(y); end

u_prop_old_bootstrap = [randn(nbr_particles_bootstrap[1], T_vec[1]+1)]
u_resample_old_bootstrap = [randn(T_vec[1], 2)]

for i in 2:M_subjects;
    append!(u_prop_old_bootstrap, [randn(nbr_particles_bootstrap[i],T_vec[i]+1)])
    append!(u_resample_old_bootstrap, [randn(T_vec[i],2)])
end

loglik, x,nbr_resamlping = @time particlefilter_test(y, exp(log_σ_ϵ), ϕ[m,:], dt, u_prop_old_bootstrap[m], u_resample_old_bootstrap[m], nbr_particles_bootstrap[m], true, true)


PyPlot.figure()
PyPlot.plot(time_grid, x[1:100,:]', "r--")
PyPlot.plot(time_grid, y, "g")


loglik, x,nbr_resamlping = diffusionbridgeparticlefilter_test(y, exp(log_σ_ϵ), ϕ[m,:], dt, u_prop_old[m], u_resample_old[m], nbr_particles[m], true, true)

PyPlot.figure()
PyPlot.plot(time_grid, y, "g")
PyPlot.plot(time_grid, x[:], "r--")


# run multiple loglik calcs

runs = 100

loglik_matrix = zeros(3, runs)

for i in 1:runs

    T_vec = zeros(Int64, M_subjects)
    for i = 1:M_subjects; T_vec[i] = length(y); end

    u_prop_old = [randn(nbr_particles[1], T_vec[1]+1)]
    u_resample_old = [randn(T_vec[1], 2)]

    for i in 2:M_subjects;
        append!(u_prop_old, [randn(nbr_particles[i],T_vec[i]+1)])
        append!(u_resample_old, [randn(T_vec[i],2)])
    end


    nbr_particles_bootstrap = 10000*ones(Int64, M_subjects)

    # set number of obs for each data set
    T_vec = zeros(Int64, M_subjects)
    for i = 1:M_subjects; T_vec[i] = length(y); end

    u_prop_old_bootstrap = [randn(nbr_particles_bootstrap[1], T_vec[1]+1)]
    u_resample_old_bootstrap = [randn(T_vec[1], 2)]

    for i in 2:M_subjects;
        append!(u_prop_old_bootstrap, [randn(nbr_particles_bootstrap[i],T_vec[i]+1)])
        append!(u_resample_old_bootstrap, [randn(T_vec[i],2)])
    end

    m = 1

    loglik_matrix[1,i], ~ = kalman_filter_test(y, exp(log_σ_ϵ),  ϕ[m,:], dt, true)
    loglik_matrix[2,i], ~ = particlefilter_test(y, exp(log_σ_ϵ), ϕ[m,:], dt, u_prop_old_bootstrap[m], u_resample_old_bootstrap[m], nbr_particles_bootstrap[m], true, true)
    loglik_matrix[3,i], ~ = diffusionbridgeparticlefilter_test(y, exp(log_σ_ϵ), ϕ[m,:], dt, u_prop_old[m], u_resample_old[m], nbr_particles[m], true, true)

end



println("Mean:")
println(round.(mean(loglik_matrix, dims = 2),digits=2))

println("Std:")
println(round.(std(loglik_matrix, dims = 2),digits=2))
