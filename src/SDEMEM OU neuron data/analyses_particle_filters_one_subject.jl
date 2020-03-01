# script to run inference for the OU SDEMEM model
using Pkg
using LinearAlgebra
using DataFrames
import Statistics.mean
import Statistics.std
using Printf

# load functions
include(pwd()*"/src/SDEMEM OU neuron data/ou_sdemem.jl")
include(pwd()*"/src/SDEMEM OU neuron data/mcmc.jl")

M_subjects = 1

seed = 100 #parse(Int,ARGS[1])

y,prior_parameters_η,prior_parameters_σ_ϵ = set_up(M_subjects,seed)


# plot data
using PyPlot

# Plot data (this plot should be similar to Fig 3 in the paper)
PyPlot.figure()
for i in 1:M_subjects
    PyPlot.plot(y[i].mSec, y[i].mV, "k", alpha = 0.2)
end
PyPlot.xlabel("Time (msec)")
PyPlot.ylabel("depolarization mV")


nbr_obs = 0
for i in 1:M_subjects
    global nbr_obs = nbr_obs + length(y[i].mSec)
end

nbr_obs

# posterior mean ests from kalman filter
posterior_mean_kalman_log_σ_ϵ = -3.60203
posterior_mean_kalman_ϕ = [-3.26986 -0.888105 -0.780767]


# loglik eval:
log_σ_ϵ = posterior_mean_kalman_log_σ_ϵ
ϕ = posterior_mean_kalman_ϕ


# set time step
dt = diff(y[1].mSec)[1]


# run kalman filter + plot paths

ll1, path1 = @time kalman(y, exp(log_σ_ϵ), ϕ, dt, true)

loglik_kalman = sum(ll1)


# bootstrap filter

function run_bootstrap(particels)

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

    run_time = @elapsed ll1,path1, nbr_resamlping =  cpf(y, exp(log_σ_ϵ), ϕ, dt, u_prop_old, u_resample_old, nbr_particles, true, true)

    return ll1, run_time, path1,nbr_resamlping

end


ll1, run_time, path1, nbr_resamlping = run_bootstrap(10000)

sum(ll1)


# bridge pf

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

    run_time = @elapsed ll1,path1,nbr_resamlping =  bridge_pf(y, exp(log_σ_ϵ), ϕ, dt, u_prop_old, u_resample_old, nbr_particles, true, true)

    return ll1,run_time,path1,nbr_resamlping
end



ll1b,run_timeb,path1b,nbr_resamlping = run_bridge_pf(100)

sum(ll1b)

subject = rand(1:M_subjects)
PyPlot.figure(figsize=(8,5))
PyPlot.subplot(121)
PyPlot.plot(path1[subject][:,:]', "r")
PyPlot.plot(y[subject].mV, "g")
PyPlot.xlabel("Time (msec)")
PyPlot.ylabel("depolarization mV")
PyPlot.title("Bootstrap")
PyPlot.subplot(122)
PyPlot.plot(path1b[subject][:,:]', "r")
PyPlot.plot(y[subject].mV, "g")
PyPlot.title("Bridge")
PyPlot.xlabel("Time (msec)")
PyPlot.ylabel("depolarization mV")
#PyPlot.savefig("figures/boostrap_bridge_tracking.pdf", format="pdf", dpi=1000)
#run(`pdftops -eps figures/boostrap_bridge_tracking.pdf figures/boostrap_bridge_tracking.eps`)

# randomly selected subject is subject 60
subject = rand(1:M_subjects)
PyPlot.figure(figsize=(12,5))
PyPlot.subplot(131)
PyPlot.plot(y[subject].mSec, y[subject].mV, "k")
PyPlot.title("Data")
PyPlot.xlabel("Time (msec)")
PyPlot.ylabel("depolarization mV")
PyPlot.subplot(132)
PyPlot.plot(y[subject].mSec,path1[subject][:,:]', "k")
PyPlot.xlabel("Time (msec)")
PyPlot.ylabel("depolarization mV")
PyPlot.title("Bootstrap")
PyPlot.subplot(133)
PyPlot.plot(y[subject].mSec,path1b[subject][:,:]', "k")
PyPlot.title("Bridge")
PyPlot.xlabel("Time (msec)")
PyPlot.ylabel("depolarization mV")
#PyPlot.savefig("figures/boostrap_bridge_tracking.pdf", format="pdf", dpi=1000)
#run(`pdftops -eps figures/boostrap_bridge_tracking.pdf figures/boostrap_bridge_tracking.eps`)

# test variance of loglik est

N_filters = 100

data = zeros(4, N_filters)
run_time_kalman = zeros(N_filters)
ll_all_boostrap = zeros(M_subjects, N_filters)
ll_all_bridge = zeros(M_subjects, N_filters)

for i in 1:N_filters

    println(i)


    ll_boostrap, run_time_boostrap, path1 = run_bootstrap(10000)
    ll1_bridge, run_time_bridge, path1 = run_bridge_pf(10)

    ll_all_boostrap[:,i] = ll_boostrap
    ll_all_bridge[:,i] = ll1_bridge

    data[:,i] = [sum(ll_boostrap);run_time_boostrap; sum(ll1_bridge);run_time_bridge]

end

for i in 1:N_filters
    run_time_kalman[i] = @elapsed ll1, path1 = kalman(y, exp(log_σ_ϵ), ϕ, dt, true)
end

@printf "----------------\n"
@printf "Loglik kalman: %.4f\n" round(loglik_kalman, digits = 0)
@printf "Kalman - runtime (mean): %.4f\n" round(mean(run_time_kalman), digits = 4)
@printf "Bootstrap (mean): %.4f\n" round(mean(data[1,:]), digits = 0)
@printf "Bootstrap (var): %.4f\n" round(var(data[1,:]), digits = 0)
@printf "Bootstrap (std): %.4f\n" round(std(data[1,:]), digits = 0)
@printf "Bootstrap (var) ind ll est:\n"
println(var(ll_all_boostrap, dims = 2))
@printf "Bootstrap - runtime (mean): %.4f\n" round(mean(data[2,:]), digits = 4)
@printf "Bridge (mean): %.4f\n" round(mean(data[3,:]), digits = 0)
@printf "Bridge (var): %.4f\n" round(var(data[3,:]), digits = 0)
@printf "Bridge (std): %.4f\n" round(std(data[3,:]), digits = 0)
@printf "Bridge (var) ind ll est:\n"
println(var(ll_all_bridge, dims = 2))
@printf "Bridge - runtime (mean): %.4f\n" round(mean(data[4,:]), digits = 4)
